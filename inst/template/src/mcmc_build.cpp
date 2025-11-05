/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build.cpp
 *
 */

#include <chrono>

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/string.h"
#include "tools/read.h"
#include "tools/error.h"
#include "tools/status.h"
#include "tools/params.h"
#include "tools/rates.h"
#include "function.h"
#include "function_registry.h"
#include "variable.h"
#include "household.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : Impl : Setup : build  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Default seed for pseudo random number generator.
std::int32_t default_setup_seed()
{
    const auto time_raw = std::chrono::steady_clock::now().time_since_epoch().count();
    // ensure time_now is 64bit
    const std::uint64_t time_now = static_cast<std::uint64_t>(time_raw);
    // xor between upper and lower part of time stamp
    const std::uint64_t mixed = time_now ^ (time_now >> 32);
    // note: constrain to 32bit for compatibility with R (mask sign bit)
    return static_cast<std::int32_t>(mixed & 0x7FFFFFFF);
}

void MCMC::Impl::Setup::build(const u_dict<std::string> &deflist)
{
    // default seed
    seed = default_setup_seed();

    if (!deflist.empty())
    {
        // check entries
        const str_vec &headers = hsp::read::list_headers(deflist);
        hsp::str::assert_no_unused_entries(headers, file_headers);

        // Random Initialization of Parameters
        if (deflist.count("progress") && !hsp::na::is(deflist.at("progress")))
            progress = hsp::str::toul(deflist.at("progress"));

        // N. of Threads
        if (deflist.count("threads") && !hsp::na::is(deflist.at("threads")))
            n_threads = hsp::str::toi(deflist.at("threads"), "n. of threads");

        // PRNG Seed
        if (deflist.count("seed") && !hsp::na::is(deflist.at("seed")))
            seed = hsp::str::toi32(deflist.at("seed"), "random seed");

        // N. of Chains
        if (deflist.count("chains") && !hsp::na::is(deflist.at("chains")))
            n_chains = hsp::str::toul(deflist.at("chains"), "n. of chains");

        // N. of Steps
        if (deflist.count("steps") && !hsp::na::is(deflist.at("steps")))
            n_steps = hsp::str::toul(deflist.at("steps"), "n. of steps");

        // N. of Adaptive Steps
        if (deflist.count("adapt") && !hsp::na::is(deflist.at("adapt")))
            n_adapt = hsp::str::toul(deflist.at("adapt"), "n. of adaptive steps");

        // N. of Thinning Steps
        if (deflist.count("thin") && !hsp::na::is(deflist.at("thin")))
            n_thin = hsp::str::toul(deflist.at("thin"), "thin");

        // N. of Steps to "Burn-In"
        if (deflist.count("burnin") && !hsp::na::is(deflist.at("burnin")))
            n_burnin = hsp::str::toul(deflist.at("burnin"), "burn-in");

        // Individual ID for Augmented Trace
        if (deflist.count("id_draws") && !hsp::na::is(deflist.at("id_draws")))
            id_draws = hsp::read::parse_list(deflist.at("id_draws"));
    }

    // reset internal states
    chain = 1; // current chain
    rung = 0;  // current thermodynamic rung

    set_n_threads(); // define number of threads
    set_prng_seed(); // pseudo-random number generator seeding
}

/// Set the number of threads for parrell computations (default is 1).
void MCMC::Impl::Setup::set_n_threads()
{
    if (n_threads == 0)
        throw std::invalid_argument("N. of threads was not greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_num_procs();

    if (n_threads > max_threads)
    {
        hsp::msg::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available threads.");
        n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n_threads);

#else
    if (n_threads != 1)
        hsp::msg::warning(
            "Multi-threading was disabled. OpenMP is not supported on current machine.");

    n_threads = 1;
#endif
}

/// Assign seed for random generators.
void MCMC::Impl::Setup::set_prng_seed()
{
    // seed for the unique thread rand num generator
    prng_.seed(seed);
    get_prng_ = &Impl::Setup::do_prng_unique_thread_;

#ifdef OPENMP
    // seeds for the case of parallelization
    if (n_threads > 1)
    {
        // create different generators for different parallel threads
        // note: they are used only to parallelize data augmentantion
        get_prng_ = &Impl::Setup::do_prng_parallel_threads_;

        std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

        prng_parallel_.reserve(n_threads);
        for (unsigned n = 0; n < n_threads; ++n)
            prng_parallel_[n].seed(seed_generator(prng_));
    }
#endif
}

void MCMC::Impl::Setup::assert_mcmc_setup()
{
    if (inv_T < 0. || inv_T > 1.)
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " was out of range [0, 1].");

    if (n_chains < 1)
        throw std::out_of_range("Number of chains was not greater than 0.");

    // Interpret thinning equal to 0 has been interpreted as 1.
    if (n_thin == 0)
        n_thin = 1;

    if (n_steps < 1)
        throw std::out_of_range(
            "N. of steps was not greater than 0.");
    if (n_burnin >= n_steps)
        throw std::out_of_range(
            "Burn-in must was not smaller than n. of steps.");
    if (n_thin > n_steps)
        throw std::out_of_range(
            "Thinning was not be smaller than (or equal to) n. of steps. ");

    c_sample_size = (n_steps - n_burnin) / n_thin;

    if (c_sample_size == 0)
        throw std::out_of_range(
            "Choice of steps, burn-in and thin generated an empty sample.");

    // assign number of sweeps
    c_n_sweeps = n_adapt / n_thin;
    c_n_sweeps += n_burnin / n_thin;
    c_n_sweeps += c_sample_size;
}

/// The text to print in the Markov chian progress bar.
std::string MCMC::Impl::Setup::fmt_progtxt()
{
    if (is_power_post())
        return "power post " + std::to_string(rung) + "/" + std::to_string(n_rungs);
    if (n_chains > 1)
        return "sample " + std::to_string(chain) + "/" + std::to_string(n_chains);
    return "sample";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : setter  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Link an individual to the household.
void Household::init_member(size_t index,
                            const std::string &id,
                            int status,
                            const u_dict<double> &data)
{
    // emplace new individual
    members_.emplace_back(index, house_size_, id, status, data);

    // update status lists within the household
    if (members_.back().is_infected_status())
        infected_list_.push_back(house_size_);
    else if (members_.back().is_healthy_status())
        healthy_list_.push_back(house_size_);
    else
        missing_list_.push_back(house_size_);

    ++house_size_;
}

/// Change_status of individual i to new status (without controls).
///	NOTE: this introduces augmented variables -> it requires RJMCMC in case of inference.
void Household::set_status_(Individual &member, int new_status)
{
    if (new_status == member.status())
        return;

    const size_t i = static_cast<size_t>(member); // get within house index

    if (member.is_infected_status())
        infected_list_.erase(std::find(infected_list_.begin(), infected_list_.end(), i));
    else if (member.is_healthy_status())
        healthy_list_.erase(std::find(healthy_list_.begin(), healthy_list_.end(), i));
    else
        missing_list_.erase(std::find(missing_list_.begin(), missing_list_.end(), i));

    member.status_ = new_status; // change member status
    if (member.is_infected_status())
        infected_list_.push_back(i);
    else if (member.is_healthy_status())
        healthy_list_.push_back(i);
    else
        missing_list_.push_back(i);
}

/// Associate a pointer to parameter value for the contact from->to.
void Household::init_contact_rates(size_t from, size_t to, const Attribute &x)
{
    // initialize contact rate matrix if necessary
    if (contact_rates_.empty())
        contact_rates_.resize(
            house_size_,
            std::vector<Attribute>(house_size_, Attribute(hsp::rates::VALUE_REFERENCE)));

    contact_rates_[from][to] = x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Recover previous values of the log-likelihood and cached depedencies.
void Household::recover_buffer_values(bool is_cache_depen)
{
    obs_log_lik_ = obs_buffer_;
    if (is_cache_depen)
    {
        // swap buffer into  cache attributes
        cache_int_vec_.swap(buffer_cache_int_vec_);
        cache_size_t_vec_.swap(buffer_cache_size_t_vec_);
        cache_dbl_vec_.swap(buffer_cache_dbl_vec_);
        cache_dbl_mtx_.swap(buffer_cache_dbl_mtx_);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC::Impl : build data  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign data after pre-processing.
void MCMC::Impl::DataWide::build(const u_dict<str_vec> &input)
{
    if (input.empty())
        return;

    // drop all previous data and prepare for new
    build_helper_clear_data_();

    // extract headers
    for (const auto &x : input)
        headers_.push_back(x.first);
    // assess the presence of house IDs
    if (!(headers_ % "house_id"))
        hsp::msg::warning(
            "Ignored household structure. Missing data column 'house_id'.");

    // safe load of households
    report_indiv.total = input.at(headers_[0]).size();
    for (size_t i = 0; i < report_indiv.total; ++i)
    {
        std::string id = std::to_string(i + 1);                          // default id
        std::string house_id = std::to_string(house_vector_.size() + 1); // default house_id
        int status = hsp::status::MISSING;
        u_dict<double> data = {};

        // read individual data
        for (const std::string &key : headers_)
            if (key == "id" || key == "house_id")
            {
                if (hsp::na::is(input.at(key)[i]))
                    throw std::invalid_argument(
                        "Individual in row n." + std::to_string(i + 2) + " had missing " + key + ".");
            }
            else if (key == "status")
            {
                // note: NA status is a missing status
                if (!hsp::na::is(input.at("status")[i]))
                    status = hsp::str::toi(input.at("status")[i], "individual status");
            }
            else
            {
                data.emplace(key, hsp::str::tod(input.at(key)[i], "data '" + key + "'"));
            }

        // if household index doesn't exist, add new household
        if (!hhid_to_index_.count(house_id))
        {
            // map data index to household number
            hhid_to_index_.emplace(house_id, house_vector_.size());
            // create a new empty household
            house_vector_.emplace_back(house_id, house_vector_.size());
        }
        // extract the household reference
        Household &house = house_vector_.at(hhid_to_index_.at(house_id));

        // assign the individual the respective household
        house.init_member(i, id, status, data);
    }
    report_house.total = house_vector_.size();

    headers_ = headers_ - str_vec({"id", "house_id", "status"});

    // update report for data
    build_helper_fill_report_data_();
}

/// Use the loaded data to fill-in the respective report.
void MCMC::Impl::DataWide::build_helper_clear_data_()
{
    headers_.clear();
    house_vector_.clear();
    hhid_to_index_.clear();
    all_infected_status_list.clear();
    report_indiv = ReportIndividual();
    report_house = ReportHousehold();
    report_hh_size = ReportHouseSizes();
}

/// Use the loaded data to fill-in the respective report.
void MCMC::Impl::DataWide::build_helper_fill_report_data_()
{
    // set report counts
    for (const Household &house : house_vector_)
    {
        // list of unique infected statuses
        for (const Individual &member : house.infected())
            if (!(all_infected_status_list % member.status()))
                all_infected_status_list.push_back(member.status());

        // observed quantities in-house
        const unsigned short house_size = house.size();
        const unsigned short n_infec = house.n_infec_status();
        const unsigned short n_miss = house.n_missing_status();
        const unsigned short n_heal = house.n_healthy_status();

        report_indiv.infected += n_infec;
        report_indiv.missing += n_miss;
        report_indiv.healthy += n_heal;

        report_house.no_cases += (n_infec == 0);
        report_house.at_least_1_case += (n_infec > 0);
        report_house.at_least_2_cases += (n_infec > 1);

#define UPDATE_REPORT_SIZE(Report, Value, Symbol) \
    if (Report.first Symbol Value)                \
    {                                             \
        Report.first = Value;                     \
        Report.second = 1;                        \
    }                                             \
    else if (Report.first == Value)               \
        ++Report.second;

        // max number of house cases
        UPDATE_REPORT_SIZE(report_hh_size.max_cases, n_infec, <);
        // max number of missing status
        UPDATE_REPORT_SIZE(report_hh_size.max_miss, n_miss, <);
        // minimum house number of members
        UPDATE_REPORT_SIZE(report_hh_size.min_size, house_size, >);
        // maximum house number of members
        UPDATE_REPORT_SIZE(report_hh_size.max_cases, house_size, <);
    }
}

/// Extra integer vectors for all households.
void MCMC::Impl::DataWide::init_cache_integer_vector(const std::string &e)
{
    for (Household &house : restricted())
        house.init_cache_integer_vector(e);
}
/// Extra size_t vectors for all households.
void MCMC::Impl::DataWide::init_cache_size_t_vector(const std::string &e)
{
    for (Household &house : restricted())
        house.init_cache_size_t_vector(e);
}
/// Extra double vectors for all households.
void MCMC::Impl::DataWide::init_cache_double_vector(const std::string &e)
{
    for (Household &house : restricted())
        house.init_cache_double_vector(e);
}
/// Extra double square matrix for all households.
void MCMC::Impl::DataWide::init_cache_double_matrix(const std::string &e)
{
    for (Household &house : restricted())
        house.init_cache_double_matrix(e);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : build multiplicities  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build multiplicities.
// UPDATE I should find a way to deal internally with multiplicity.
void MCMC::Impl::Multiplicity::build(const u_dict<str_vec> &deftable)
{
    indexes.clear();
    size_t n_houses = pImpl_.data_wide.report_house.total;
    if (pImpl_.augm.is_any() || (deftable.size() == 0))
    {
        // note: multiplicities are not compatible with data augmentation
        if (!deftable.empty())
            hsp::msg::warning(
                "Ignored multiplicity optimization: incompatible with data augmentation.");

        // default choice for restricted household indexes
        indexes.resize(n_houses);
        for (size_t H = 0; H < n_houses; ++H)
        {
            pImpl_.data_wide(H).set_multiplicity(1, H);
            indexes[H] = H;
        }
    }
    else
    {
        // UPDATE: this error is purely internal
        if (deftable.at("house_id").size() != n_houses)
            throw std::out_of_range(
                "Multiplicity definitions had different dimensions than data.");

        for (size_t r = 0; r < n_houses; ++r)
        {
            Household &house = pImpl_.data_wide(deftable.at("house_id")[r]);
            size_t H0 = pImpl_.data_wide(deftable.at("class_id")[r]).house_index();
            unsigned long m = hsp::str::toul(deftable.at("multiplicity")[r]);

            house.set_multiplicity(m, H0);

            if (house.is_multiplicity_head_class())
                indexes.emplace_back(house.house_index());
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : DataLong : build  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign data after pre-processing.
void MCMC::Impl::DataLong::build(const u_dict<str_vec> &input)
{
    if (input.empty())
        return;

    // data sequence headers
    headers_.clear();
    for (const auto &x : input)
        headers_.push_back(x.first);

    // check for individual IDs
    if (!(headers_ % "id"))
        throw std::invalid_argument(
            "Individual identifier column 'id' was missing.");
    // remove "id"
    headers_ = headers_ - str_vec({"id"});

    // individual map for id -> index
    // note: this could be moved elsewhere if required
    std::unordered_map<std::string, size_t> id_to_index = {};
    size_t c = 0;
    for (const Household &house : pImpl_.data_wide)
        for (const Individual &member : house)
        {
            if (id_to_index.count(member.id()))
                throw std::invalid_argument(
                    "There were multiple individuals with id equal to '" + member.id() + "'.");
            id_to_index.emplace(member.id(), c++);
        }

    // data-long allocation
    indiv_datalist_.clear();
    indiv_datalist_.resize(pImpl_.data_wide.report_indiv.total);

    for (u_dict<dbl_vec> &xmap : indiv_datalist_)
        for (const std::string &key : headers_)
            xmap[key].clear();

    // load individual data list
    const size_t n_rows = input.at(headers_[0]).size();
    for (size_t r = 0; r < n_rows; ++r)
    {
        std::string id = input.at("id")[r];
        if (hsp::na::is(id))
            throw std::invalid_argument(
                "Invalid id found in row n." + std::to_string(r + 2) + ".");
        if (!id_to_index.count(id))
        {
            // UPDATE: these errors should be more self-explicative
            if (pImpl_.data_wide.is_empty())
                throw std::invalid_argument(
                    "Tried to define data long format before wide format.");
            throw std::invalid_argument(
                "Id '" + id + "' from data long format had no match in loaded data wide format.");
        }

        // append (since vector is initialized to empty)
        u_dict<dbl_vec> &xdict = indiv_datalist_.at(id_to_index.at(id));
        for (const std::string &key : headers_)
            xdict.at(key).emplace_back(
                hsp::str::tod(input.at(key).at(r), "data long format '" + key + "'"));
    }

    // assign pointers to relative data-long to all individuals
    c = 0;
    for (Household &house : pImpl_.data_wide)
        for (Individual &member : house)
            member.init_pointer_to_data_long(&indiv_datalist_.at(c++));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : GroupRates : build  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Update model report and set group rates associate to new individual.
void MCMC::Impl::GroupRates::build(const u_dict<str_vec> &deftable)
{
    if (deftable.empty())
        return;

    const str_vec data_headers(helper_assert_definitions_(deftable));

    // loop over all household members
    const size_t n_rows = deftable.at("group_rate").size();
    for (Household &house : pImpl_.data_wide)
        for (Individual &member : house)
        {
            u_dict<str_vec> these_rates = {};

            // obtain group rates for given member
            for (size_t r = 0; r < n_rows; ++r)
            {
                // check if all conditions are satisfied
                bool is_this_class = true;
                for (const std::string &d : data_headers)
                    if (!hsp::na::is(deftable.at(d)[r]))
                        is_this_class &= (member.data(d) == hsp::str::tod(deftable.at(d)[r]));

                // if it's the case, add the parameter name
                if (is_this_class)
                    these_rates[deftable.at("group_rate")[r]].push_back(deftable.at("param_name")[r]);
            }

            // assign rates to the individual
            helper_init_all_rates_(member, these_rates);
        }
}

/// Verify that group rates definitions are correct.
str_vec MCMC::Impl::GroupRates::helper_assert_definitions_(const u_dict<str_vec> &deftable)
{
    // check that tempalte headers are present
    str_vec data_headers = hsp::read::table_headers(deftable);
    hsp::str::assert_no_missing_entries(data_headers, file_headers_);

    // remove template headers to keep only data
    data_headers = data_headers - file_headers_;
    // check that required data columns are present
    hsp::str::assert_no_unused_entries(data_headers, pImpl_.data_wide.headers());

    // assert the group rate names are allowed
    headers_.clear();
    for (const std::string &g : deftable.at("group_rate"))
    {
        if (g == hsp::rates::CONTACT || g == hsp::rates::MISSING)
            throw std::invalid_argument(
                "Tried to use the reserved group rate name '" + g + "'.\nPlease, choose a different name.");
        if (hsp::na::is(g))
            throw std::invalid_argument(
                "Found a missing group rate name in the definitions.");
        headers_.push_back(g);
    }
    headers_ = hsp::str::unique(headers_);

    return data_headers;
}

/// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
void MCMC::Impl::GroupRates::helper_init_all_rates_(Individual &member,
                                                    const u_dict<str_vec> &these_rates)
{
    for (const auto &pair : these_rates)
    {
        const std::string &g = pair.first;

        // extract non-reference parameter pointers
        std::vector<const Attribute *> vec_par_ptr = {};
        for (const std::string &p : pair.second)
        {
            if (hsp::na::is(p) || p == hsp::rates::MISSING)
                throw std::invalid_argument(
                    "A group rate defintion '" + g + "' was missing for " + member.err_name() + ".");
            if (p != hsp::rates::REFERENCE)
                vec_par_ptr.emplace_back(pImpl_.param.at(p).get_ptr());
        }

        if (vec_par_ptr.empty()) // stick to reference
            member.init_group_rate_to_reference(g);

        else if (vec_par_ptr.size() == 1) // use simple pointer
            member.init_group_rate(g, Attribute(vec_par_ptr[0]));

        else // use variadic product
        {
            std::shared_ptr<Function> fun = compose_func_registries(
                "prod", "function", {&ArithmeticRegistry::inst()}, false);

            // initialize attributes to pointers
            fun->init_attr(hsp::params::X, Attribute(vec_par_ptr[0]));
            for (size_t i = 1; i < vec_par_ptr.size(); ++i)
                fun->init_attr(hsp::params::X + std::to_string(i + 1), Attribute(vec_par_ptr[i]));

            member.init_group_rate(g, Attribute(fun));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : Parameter : build  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign parameter definitions to dictionary.
void MCMC::Impl::Parameter::build(const u_dict<str_vec> &deftable)
{
    if (deftable.empty())
        return;

    // check entries
    const str_vec file_headers = {"name", "value", "const", "discr",
                                  "prior", hsp::params::MIN, hsp::params::MAX,
                                  "rule", "scale"};
    const str_vec headers = hsp::read::table_headers(deftable);
    hsp::str::assert_no_missing_entries(headers, file_headers);
    hsp::str::assert_no_unused_entries(headers, file_headers);

    // drop current list of parameters
    dictionary_.clear();
    list_const.clear();
    list_sample.clear();
    list_sample_standard.clear();
    list_adaptable_sample_standard.clear();
    list_non_adapt_sample_standard.clear();
    list_sample_cached.clear();
    list_adaptable_sample_cached.clear();
    list_non_adapt_sample_cached.clear();

    // assert that parameter is not already defined
    if (deftable.at("name").size() != hsp::str::unique(deftable.at("name")).size())
        throw std::invalid_argument(
            "Found multiple definitions for the same parameters.");

    // load parameters from the row with options
    for (size_t k = 0; k < deftable.at("name").size(); ++k)
        if (hsp::str::tobool(deftable.at("const")[k]))
            build_constant_parameter_(deftable, k);
        else
            build_parameter_to_sample_(deftable, k);

    // initialize all parameters to standard
    list_sample_standard = list_sample;
    for (const std::string &p : list_sample_standard)
        if (this->at(p).is_adaptable())
            list_adaptable_sample_standard.push_back(p);
        else
            list_non_adapt_sample_standard.push_back(p);
}

/// Helper to build constant parameter.
void MCMC::Impl::Parameter::build_constant_parameter_(const u_dict<str_vec> &deftable, size_t k)
{
    const std::string &p = deftable.at("name")[k]; // parameter name

    if (hsp::na::is(deftable.at("value")[k]))
        throw std::invalid_argument(
            "Constant parameter '" + p + "' had undefined value.");

    // assign constant parameter to the set of variables
    double value = pImpl_.read_numeric_(deftable.at("value")[k]).get();
    dictionary_.emplace(p, Variable(p, "parameter", value));
    // update the lists of constant parameters
    list_const.push_back(p);
}

/// Helper to build parameter to sample.
void MCMC::Impl::Parameter::build_parameter_to_sample_(const u_dict<str_vec> &deftable, size_t k)
{
    const std::string &p = deftable.at("name")[k]; // parameter name

    // assert parameter to sample has a valid value initialization
    if (hsp::na::is(deftable.at("value")[k]) && hsp::na::is(deftable.at("prior")[k]))
        throw std::invalid_argument(
            "Non-constant parameter '" + p + "' had undefined prior.");

    // load rule name
    std::string rule_name = "sample";
    if (!hsp::na::is(deftable.at("rule")[k]))
        rule_name = deftable.at("rule")[k];

    // load rule scale
    double rule_scale = 1.;
    if (!hsp::na::is(deftable.at("scale")[k]))
        rule_scale = hsp::str::tod(deftable.at("scale")[k]);

    // create a new parameter Variable
    dictionary_.emplace(
        p, Variable(p, "parameter", pImpl_.read_numeric_(deftable.at("value")[k]).get(), false));
    Variable &new_param = this->at(p);

    // assign a prior distribution to the parameter
    hsp::parse::Definition parser(deftable.at("prior")[k]);
    new_param.init_distr(parser);
    // load prior arguments
    for (const auto &child : parser.root.children)
        new_param.init_distr_attr(child.key, pImpl_.read_numeric_(child));
    // load parameter range
    new_param.init_distr_attr(hsp::params::MIN, pImpl_.read_numeric_(deftable.at(hsp::params::MIN)[k]));
    new_param.init_distr_attr(hsp::params::MAX, pImpl_.read_numeric_(deftable.at(hsp::params::MAX)[k]));
    new_param.init_distr_reshaping(hsp::str::tobool(deftable.at("discr")[k]));

    // assign the rule for the proposer to the parameter
    new_param.init_proposer(rule_name, rule_scale);
    // update the lists of parameters to sample
    list_sample.push_back(p);
}

/// Adjust lists of parameters to be inferred making parameter as dependent on cache.
void MCMC::Impl::Parameter::make_cache_dependent(const std::string &p)
{
    // remove a parameter from standard if it was defined as such
    if (list_sample_standard % p)
    {
        list_sample_standard = list_sample_standard - str_vec({p});
        list_adaptable_sample_standard = list_adaptable_sample_standard - str_vec({p});
        list_non_adapt_sample_standard = list_non_adapt_sample_standard - str_vec({p});
    }

    // add the parameter as cache if not alread present
    if (!(list_sample_cached % p))
    {
        list_sample_cached.push_back(p);
        if (at(p).is_adaptable())
            list_adaptable_sample_cached.push_back(p);
        else
            list_non_adapt_sample_cached.push_back(p);
    }
}

/// Randomly initialize values of the parameters to be inferred.
void MCMC::Impl::Parameter::user_reset()
{
    for (const std::string &p : list_sample)
        this->at(p).user_reset(pImpl_.setup.prng());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : constructor  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC()
    : pImpl_(std::make_unique<Impl>(
          [this]()
          { this->init_model_for_inference(); },
          [this]()
          { this->init_model_for_simulation(); },
          [this](Household &h) -> double
          { return this->compute_household_obs_log_lik(h); },
          [this](Household &h)
          { return this->simulate_single_household(h); },
          [this](Household &h, Individual &i)
          { this->update_cache_dependecies(h, i); })) {}
MCMC::~MCMC() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::uniform() { return pImpl_->setup.unif(); }
double MCMC::param(const std::string &p) const { return pImpl_->param.get(p); }
const ReportIndividual &MCMC::get_report_indiv() const { return pImpl_->data_wide.report_indiv; }
const ReportHousehold &MCMC::get_report_house() const { return pImpl_->data_wide.report_house; }
const ReportHouseSizes &MCMC::get_report_hh_size() const { return pImpl_->data_wide.report_hh_size; }
void MCMC::assert_no_missing_status() const
{
    if (pImpl_->data_wide.report_indiv.missing > 0)
        hsp::msg::error(
            "Individual with missing status are incompatible with current model.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : cache attributes  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::cache_house_integer_vector(const std::string &e) { pImpl_->data_wide.init_cache_integer_vector(e); }
void MCMC::cache_house_size_t_vector(const std::string &e) { pImpl_->data_wide.init_cache_size_t_vector(e); }
void MCMC::cache_house_double_vector(const std::string &e) { pImpl_->data_wide.init_cache_double_vector(e); }
void MCMC::cache_house_double_matrix(const std::string &e) { pImpl_->data_wide.init_cache_double_matrix(e); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : contact rates ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build all group rates and contact rates and set rate per group. // UPDATE
void MCMC::build_all_rates()
{
    try
    {
        for (Household &house : pImpl_->data_wide)
            init_contact_rates_(house);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                        std::string(e.what()));
    }
}

/// Set contact rate according to user assignement. // UPDATE
void MCMC::init_contact_rates_(Household &house)
{
    // loop over all couples
    for (const Individual &member_k : house)     // from : member_k
        for (const Individual &member_i : house) // -> to : member_i
            if (member_k != member_i)
            {
                _pair_contact_rate = hsp::na::STRING;
                init_household_contact_rates(house, member_k, member_i);

                if (!hsp::na::is(_pair_contact_rate))
                {
                    // default reference rate
                    Attribute y(hsp::rates::VALUE_REFERENCE);

                    if (_pair_contact_rate == hsp::rates::MISSING)
                        throw std::invalid_argument(
                            "Missing contact rate for individuals '" +
                            member_k.id() + "' and '" + member_i.id() +
                            "' (" + house.err_name() + ").");

                    else if (_pair_contact_rate != hsp::rates::REFERENCE)
                        y = pImpl_->param.at(_pair_contact_rate); // parameter rate

                    house.init_contact_rates(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
                }
            }
}
