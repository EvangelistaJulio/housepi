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

#include "tools/string.h"
#include "tools/read.h"
#include "tools/error.h"
#include "tools/status.h"
#include "tools/params.h"
#include "tools/rates.h"
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
        if (deflist % "progress" && !hsp::na::is(deflist.at("progress")))
            progress = hsp::str::toul(deflist.at("progress"));

        // N. of Threads
        if (deflist % "threads" && !hsp::na::is(deflist.at("threads")))
            n_threads = hsp::str::toi(deflist.at("threads"), "n. of threads");

        // PRNG Seed
        if (deflist % "seed" && !hsp::na::is(deflist.at("seed")))
            seed = hsp::str::toi32(deflist.at("seed"), "random seed");

        // N. of Chains
        if (deflist % "chains" && !hsp::na::is(deflist.at("chains")))
            n_chains = hsp::str::toul(deflist.at("chains"), "n. of chains");

        // N. of Steps
        if (deflist % "steps" && !hsp::na::is(deflist.at("steps")))
            n_steps = hsp::str::toul(deflist.at("steps"), "n. of steps");

        // N. of Adaptive Steps
        if (deflist % "adapt" && !hsp::na::is(deflist.at("adapt")))
            n_adapt = hsp::str::toul(deflist.at("adapt"), "n. of adaptive steps");

        // N. of Thinning Steps
        if (deflist % "thin" && !hsp::na::is(deflist.at("thin")))
            n_thin = hsp::str::toul(deflist.at("thin"), "thin");

        // N. of Steps to "Burn-In"
        if (deflist % "burnin" && !hsp::na::is(deflist.at("burnin")))
            n_burnin = hsp::str::toul(deflist.at("burnin"), "burn-in");

        // Individual ID for Augmented Trace
        if (deflist % "id_draws" && !hsp::na::is(deflist.at("id_draws")))
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

#if defined(_OPENMP)
    unsigned int max_threads = omp_get_num_procs();

    if (n_threads > max_threads)
    {
        hsp::msg::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available processors.");
        n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n_threads);

#else
    if (n_threads != 1)
        hsp::msg::warning("Parallelization with OpenMP is not supported.");

    n_threads = 1;
#endif
}

/// Assign seed for random generators.
void MCMC::Impl::Setup::set_prng_seed()
{
    // seed for the unique thread rand num generator
    prng_.seed(seed);
    get_prng_ = &Impl::Setup::do_prng_unique_thread_;

#if defined(_OPENMP)
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
    if ((inv_T < 0.) || (inv_T > 1.))
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
    else if (n_chains > 1)
        return "sample " + std::to_string(chain) + "/" + std::to_string(n_chains);
    else
        return "sample";
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC::Impl : build data  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign data after pre-processing.
void MCMC::Impl::DataWide::build(const u_dict<str_vec> &input)
{
    if (!input.empty())
    {
        // drop all previous data and prepare for new
        headers_.clear();
        for (const auto &x : input)
            headers_.push_back(x.first);

        n_individuals_ = input.at(headers_[0]).size();

        house_vector_.clear();
        restrict_indexes_.clear();
        hhid_to_index_.clear();

        // assess the presence of house IDs
        if (!(headers_ % "house_id"))
            hsp::msg::warning(
                "Ignored household structure. Missing data column \"house_id\".");

        // safe load of households
        for (size_t i = 0; i < n_individuals_; ++i)
        {
            // read individual data
            std::string id = std::to_string(i + 1);                // default id
            std::string house_id = std::to_string(n_houses() + 1); // default house_id
            int status = hsp::status::MISSING;
            u_dict<double> data = {};

            for (const std::string &key : headers_)
                if (key == "id")
                {
                    id = input.at("id")[i];
                    if (hsp::na::is(id))
                        throw std::invalid_argument(
                            "Individual in row n." + std::to_string(i + 2) + " had missing id.");
                }
                else if (key == "house_id")
                {
                    house_id = input.at("house_id")[i];
                    if (hsp::na::is(house_id))
                        throw std::invalid_argument(
                            "Individual in row n." + std::to_string(i + 2) + " had missing house_id.");
                }
                else if (key == "status")
                {
                    // note: NA status is a missing status
                    if (!hsp::na::is(input.at("status")[i]))
                        status = hsp::str::toi(input.at("status")[i], "individual status");
                }
                else
                {
                    data.emplace(key, hsp::str::tod(input.at(key)[i], "data \"" + key + "\""));
                }

            // if household index doesn't exist, add new household
            if (!(hhid_to_index_ % house_id))
            {
                // map data index to household number
                hhid_to_index_.emplace(house_id, n_houses());
                // default choice for restricted indeces
                restrict_indexes_.push_back(n_houses());
                // create a new empty household
                house_vector_.emplace_back(house_id, n_houses());
            }

            // extract the household reference
            Household &house = house_vector_.at(hhid_to_index_.at(house_id));

            // assign the individual the respective household
            house.init_member(i, id, status, data);
        }
        headers_ = headers_ - str_vec({"id", "house_id", "status"});

        // update report for data
        build_helper_fill_report_data_();
    }
}

/// Use the loaded data to fill-in the respective report.
void MCMC::Impl::DataWide::build_helper_fill_report_data_()
{
    all_infected_status_list.clear();

    report_indiv.total = n_individuals_;
    report_indiv.missing = 0;
    report_indiv.healthy = 0;
    report_indiv.infected = 0;

    report_house.total = n_houses();
    report_house.no_cases = 0;
    report_house.at_least_1_case = 0;
    report_house.at_least_2_cases = 0;

    report_hh_size.max_cases = {0, 0};
    report_hh_size.max_miss = {0, 0};
    report_hh_size.min_size = {1000, 0};
    report_hh_size.max_size = {0, 0};

    // set report counts
    for (const Household &house : house_vector_)
    {
        // list of unique infected statuses
        for (const Individual &member : house)
            if (member.is_infected_status())
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
void MCMC::Impl::DataWide::build_multiplicity(const u_dict<str_vec> &table)
{
    if (pImpl_.augm.is_any() || (table.size() == 0))
    {
        // note: multiplicities are not compatible with data augmentation
        if (!table.empty())
            hsp::msg::warning(
                "Ignored multiplicity optimization: incompatible with data augmentation.");

        // default choice for restricted household indexes
        restrict_indexes_.clear();
        restrict_indexes_.resize(n_houses());
        for (size_t H = 0; H < n_houses(); ++H)
        {
            house_vector_[H].set_multiplicity(1, H);
            restrict_indexes_[H] = H;
        }
    }
    else
    {
        // UPDATE: this error is purely internal
        if (n_houses() != table.at("house_id").size())
            throw std::out_of_range(
                "Multiplicity definitions had different dimensions than dataset.");

        restrict_indexes_.clear();
        for (size_t r = 0; r < n_houses(); ++r)
        {
            size_t H = hhid_to_index_.at(table.at("house_id")[r]);
            size_t H0 = hhid_to_index_.at(table.at("class_id")[r]);
            unsigned long m = hsp::str::toul(table.at("multiplicity")[r]);

            house_vector_[H].set_multiplicity(m, H0);

            if (house_vector_[H].is_multiplicity_head_class())
                restrict_indexes_.emplace_back(H);
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
    if (!input.empty())
    {
        // data sequence headers
        headers_.clear();
        for (const auto &x : input)
            headers_.push_back(x.first);

        // check for individual IDs
        if (!(headers_ % "id"))
            throw std::invalid_argument(
                "Individual identifier column \"id\" was missing.");
        // remove "id"
        headers_ = headers_ - str_vec({"id"});

        // individual map for id -> index
        // note: this could be moved elsewhere if required
        std::unordered_map<std::string, size_t> id_to_index = {};
        size_t c = 0;
        for (const Household &house : pImpl_.dataset)
            for (const Individual &member : house)
                if (id_to_index % member.id())
                    throw std::invalid_argument(
                        "There were multiple individuals with id equal to \"" + member.id() + "\".");
                else
                    id_to_index.emplace(member.id(), c++);

        // datalong allocation
        if (indiv_datalist_.empty())
            indiv_datalist_.resize(pImpl_.dataset.n_individuals());

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
            else if (!(id_to_index % id))
            {
                // UPDATE: these errors should be more self-explicative
                if (pImpl_.dataset.is_empty())
                    throw std::invalid_argument(
                        "Tried to define data long format before wide format.");
                else
                    throw std::invalid_argument(
                        "Id \"" + id + "\" from data long format had no match in loaded data wide format.");
            }

            // append (since vector is initialized to empty)
            u_dict<dbl_vec> &xdict = indiv_datalist_.at(id_to_index.at(id));
            for (const std::string &key : headers_)
                xdict.at(key).emplace_back(
                    hsp::str::tod(input.at(key).at(r), "datalong \"" + key + "\""));
        }

        // assign pointers to relative datalong to all individuals
        c = 0;
        for (Household &house : pImpl_.dataset)
            for (Individual &member : house)
                member.init_pointer_to_datalong(&indiv_datalist_.at(c++));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : GroupRates : build  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Update model report and set group rates associate to new individual.
void MCMC::Impl::GroupRates::build(const u_dict<str_vec> &deftable)
{
    if (!deftable.empty())
    {
        const str_vec data_headers(helper_assert_definitions_(deftable));

        // loop over all household members
        const size_t n_rows = deftable.at("group_rate").size();
        for (Household &house : pImpl_.dataset)
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
}

/// Verify that group rates definitions are correct.
str_vec MCMC::Impl::GroupRates::helper_assert_definitions_(const u_dict<str_vec> &deftable)
{
    // check that tempalte headers are present
    str_vec data_headers = hsp::read::table_headers(deftable);
    hsp::str::assert_no_missing_entries(data_headers, file_headers_);

    // remove template headers to keep only dataset
    data_headers = data_headers - file_headers_;
    // check that required dataset columns are present
    hsp::str::assert_no_unused_entries(data_headers, pImpl_.dataset.headers());

    // assert the group rate names are allowed
    headers_.clear();
    for (const std::string &g : deftable.at("group_rate"))
        if (g == hsp::rates::CONTACT || g == hsp::rates::MISSING)
            throw std::invalid_argument(
                "Tried to use the reserved group rate name \"" + g + "\".\nPlease, choose a different name.");
        else if (hsp::na::is(g))
            throw std::invalid_argument(
                "Found a missing group rate name in the definitions.");
        else
            headers_.push_back(g);
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
        std::vector<observ_dbl_ptr> vec_par_ptr = {};
        for (const std::string &p : pair.second)
            if (hsp::na::is(p) || p == hsp::rates::MISSING)
                throw std::invalid_argument(
                    "A group rate defintion \"" + g + "\" was missing for " + member.err_name() + ".");
            else if (p != hsp::rates::REFERENCE)
                vec_par_ptr.emplace_back(pImpl_.param.get_ptr(p));

        if (vec_par_ptr.empty()) // stick to reference
            member.init_group_rate_to_reference(g);

        else if (vec_par_ptr.size() == 1) // use simple pointer
            member.init_group_rate(g, Attribute(vec_par_ptr[0]));

        else // use variadic product
        {
            std::shared_ptr<Function> fun = compose_func_registries(
                "prod", "function", {&ArithmeticRegistry::inst()});

            // initialize attributes to pointers
            fun->init_attr(hsp::params::X, vec_par_ptr[0]);
            for (size_t i = 1; i < vec_par_ptr.size(); ++i)
                fun->init_attr(hsp::params::X + std::to_string(i + 1), vec_par_ptr[i]);

            member.init_group_rate(g, fun);
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
    if (!deftable.empty())
    {
        const str_vec file_headers = {"name", "value",
                                      "const", "discr",
                                      "prior",
                                      hsp::params::MIN, hsp::params::MAX,
                                      "rule", "scale"};

        // check entries
        str_vec headers = hsp::read::table_headers(deftable);
        hsp::str::assert_no_missing_entries(headers, file_headers);
        hsp::str::assert_no_unused_entries(headers, file_headers);

        // drop current list of parameters
        dictionary_.clear();
        list_const.clear();
        list_sample.clear();
        list_sample_standard.clear();
        list_sample_cache.clear();

        // assert that parameter is not already defined
        if (deftable.at("name").size() != hsp::str::unique(deftable.at("name")).size())
            throw std::invalid_argument(
                "Found multiple definitions for the same parameters.");

        // load parameters from the row with options
        for (size_t k = 0; k < deftable.at("name").size(); ++k)
            if (hsp::str::tobool(deftable.at("const")[k]))
                build_constant_parameter_(deftable, k);
            else
            {
                list_sample.push_back(deftable.at("name")[k]);
                build_parameter_to_sample_(deftable, k);
            }

        // initialize all parameters to standard
        list_sample_standard = list_sample;
    }
}

/// Helper to build constant parameter.
void MCMC::Impl::Parameter::build_constant_parameter_(const u_dict<str_vec> &deftable, size_t k)
{
    const std::string &p = deftable.at("name")[k]; // parameter name

    if (hsp::na::is(deftable.at("value")[k]))
        throw std::invalid_argument(
            "Constant parameter \"" + p + "\" had undefined value.");

    // assign constant parameter to the set of variables
    double value = pImpl_.read_numeric_(deftable.at("value")[k]);
    dictionary_.emplace(p, Variable(p, "parameter", value));
    // update the lists of constant parameters
    list_const.push_back(p);
}

/// Helper to build parameter to sample.
void MCMC::Impl::Parameter::build_parameter_to_sample_(const u_dict<str_vec> &deftable, size_t k)
{
    const std::string &p = deftable.at("name")[k]; // parameter name

    // assert parameter to sample has a valid value initialization
    if (hsp::na::is(deftable.at("value")[k]))
        if (hsp::na::is(deftable.at("prior")[k]))
            throw std::invalid_argument(
                "Non-constant parameter \"" + p + "\" had undefined prior.");

    // load boolean flag for discrete variables
    bool is_discrete = hsp::str::tobool(deftable.at("discr")[k]);

    // load rule name
    std::string rule_name = "sample";
    if (!hsp::na::is(deftable.at("rule")[k]))
        rule_name = deftable.at("rule")[k];

    // load rule scale
    double rule_scale = 1.;
    if (!hsp::na::is(deftable.at("scale")[k]))
        rule_scale = hsp::str::tod(deftable.at("scale")[k]);

    // create a new Parameter Variable
    dictionary_.emplace(
        p, Variable(p, "parameter", pImpl_.read_numeric_(deftable.at("value")[k]), false));
    Variable &new_param = this->at(p);

    // assign a prior distribution to the parameter
    DefinitionParser parser(deftable.at("prior")[k]);
    new_param.init_distr(parser);

    for (const auto &child : parser.root.children)
        new_param.init_distr_attr(child.key, pImpl_.read_numeric_(child));

    // load parameter range
    new_param.init_distr_attr(hsp::params::MIN,
                              pImpl_.read_numeric_(deftable.at(hsp::params::MIN)[k]));
    new_param.init_distr_attr(hsp::params::MAX,
                              pImpl_.read_numeric_(deftable.at(hsp::params::MAX)[k]));
    new_param.init_distr_reshaping(is_discrete);

    // assign the rule for the proposer to the parameter
    new_param.init_proposer(rule_name, rule_scale);
}

/// Adjust lists of parameters to be inferred making parameter as dependent on cache.
void MCMC::Impl::Parameter::make_cache_dependent(const std::string &p)
{
    // remove a parameter from standard if it was defined as such
    if (list_sample_standard % p)
        list_sample_standard = list_sample_standard - str_vec({p});

    // add the parameter as cache
    if (!(list_sample_cache % p))
        list_sample_cache.push_back(p);
}

/// Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::Parameter::user_reset()
{
    for (const std::string &p : list_sample)
        this->at(p).user_reset(pImpl_.setup.prng());
}
