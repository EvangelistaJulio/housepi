/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build.cpp
 *
 */

#include <chrono>

#include "mcmc_impl.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (setter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::Impl::ImplSetup::build(const dict_str_str &deflist)
{
    try
    {
        if (!deflist.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deflist(deflist);
            tools::check_for_unused_entries(headers, file_headers);

            // Random Initialization of Parameters
            if (deflist.count("progress") && (!tools::is_na(deflist.at("progress"))))
                progress = tools::stoul(deflist.at("progress"));

            // N. of Threads
            if (deflist.count("n_threads") && (!tools::is_na(deflist.at("n_threads"))))
                n_threads = tools::stoi(deflist.at("n_threads"), "n. of threads");

            // PRNG Seed
            if (deflist.count("seed") && (!tools::is_na(deflist.at("seed"))))
                seed = tools::stoi32(deflist.at("seed"), "random seed");

            // N. of Chains
            if (deflist.count("n_chains") && (!tools::is_na(deflist.at("n_chains"))))
                n_chains = tools::stoul(deflist.at("n_chains"), "n. of chains");

            // Single Chain Sample Size
            if (deflist.count("sample_size") && (!tools::is_na(deflist.at("sample_size"))))
                c_sample_size = tools::stoul(deflist.at("sample_size"), "sample size");

            // N. of Steps
            if (deflist.count("n_steps") && (!tools::is_na(deflist.at("n_steps"))))
                n_steps = tools::stoul(deflist.at("n_steps"), "n. of steps");

            // N. of Adaptive Steps
            if (deflist.count("n_adapt") && (!tools::is_na(deflist.at("n_adapt"))))
                n_adapt = tools::stoul(deflist.at("n_adapt"), "n. of adaptive steps");

            // N. of Thinning Steps
            if (deflist.count("thin") && (!tools::is_na(deflist.at("thin"))))
                thin = tools::stoul(deflist.at("thin"), "thin");

            // N. of Steps to "Burn-In"
            if (deflist.count("n_burnin") && (!tools::is_na(deflist.at("n_burnin"))))
                n_burnin = tools::stoul(deflist.at("n_burnin"), "burn-in");

            // Random Initialization of Parameters
            if (deflist.count("rand_init") && (!tools::is_na(deflist.at("rand_init"))))
                rand_init_params = tools::stobool(deflist.at("rand_init"));

            // Individual ID for Augmented Trace
            if (deflist.count("id_trace") && (!tools::is_na(deflist.at("id_trace"))))
                id_trace = tools::parse_list_args(deflist.at("id_trace"));
        }

        // reset internal states
        chain = 1;   // current chain
        rung = 0;    // current thermodynamic rung
        n_rungs = 0; // number of rungs

        set_n_threads(); // define number of threads
        set_prng_seed(); // pseudo-random number generator seeding
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing setup definitions:\n" +
                     std::string(e.what()));
    }
}

/* Set the number of threads for parrell computations (default is 1). */
void MCMC::Impl::ImplSetup::set_n_threads()
{
    if (n_threads == 0)
        throw std::invalid_argument("N. of threads must be greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_num_procs();

    if (n_threads > max_threads)
    {
        tools::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available processors.");
        n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n_threads);

#else
    if (n_threads != 1)
        tools::warning("Parallelization with OpenMP not supported.");

    n_threads = 1;
#endif
}

/* Assign seed for random generators. */
void MCMC::Impl::ImplSetup::set_prng_seed()
{
    // seed for the unique thread rand num generator
    prng_.seed(seed);
    get_prng_ = &Impl::ImplSetup::do_prng_unique_thread_;

#ifdef _OPENMP
    // seeds for the case of parallelization
    if (n_threads > 1)
    {
        // create different generators for different parallel threads
        // note: they are used only to parallelize data augmentantion
        get_prng_ = &Impl::ImplSetup::do_prng_parallel_threads_;

        std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

        prng_parallel_.reserve(n_threads);
        for (unsigned n = 0; n < n_threads; ++n)
            prng_parallel_[n].seed(seed_generator(prng_));
    }
#endif
}

void MCMC::Impl::ImplSetup::assert_mcmc_setup()
{
    if ((inv_T < 0.) || (inv_T > 1.))
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " out of range [0, 1].");

    if ((n_chains > 1) && (rand_init_params == false))
    {
        tools::warning(
            "Unable to fix initial parameters when running multiple chains.");
        rand_init_params = true;
    }
    else if (n_chains < 1)
        throw std::out_of_range("N. of chains must be greater than 0.");

    if (thin == 0)
    {
        tools::warning("Thinning equal to 0 interpreted as 1.");
        thin = 1;
    }

    if (c_sample_size == 0)
    {
        if (n_steps < 1)
            throw std::out_of_range(
                "N. of steps must be greater than 0.");
        if (n_burnin >= n_steps)
            throw std::out_of_range(
                "Burn-in must be smaller than n. of steps.");
        if (thin > n_steps)
            throw std::out_of_range(
                "Thinning must be smaller than (or equal to) n. of steps. ");

        c_sample_size = (n_steps - n_burnin) / thin;

        if (c_sample_size == 0)
            throw std::out_of_range(
                "Choice of steps, burn-in and thin return a 0 sample size.");
    }
    else
        n_steps = n_burnin + c_sample_size * thin;

    // assign number of sweeps
    c_n_sweeps = n_adapt / thin;
    c_n_sweeps += n_burnin / thin;
    c_n_sweeps += c_sample_size;
}

/* The text to print in the Markov chian progress bar. */
const std::string MCMC::Impl::ImplSetup::fmt_progtxt()
{
    if (is_power_post())
        return "power post " + std::to_string(rung) + "/" + std::to_string(n_rungs);
    else if (n_chains > 1)
        return "sample " + std::to_string(chain) + "/" + std::to_string(n_chains);
    else
        return "sample";
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Assert dataset is not empty. */
void MCMC::Impl::assert_proper_dataset_()
{
    if (dataset.n_houses() == 0)
        throw std::length_error(
            "Unable to proceed. No database was found.");
}

/* Assert setup definitions are ok to proceed with mcmc sampling. */
void MCMC::Impl::init_monte_carlo_()
{
    setup.assert_mcmc_setup();

    // check if there is a dataset
    assert_proper_dataset_();

    // check if there is actually something to sample
    if (!(param.is_sample() || augm.is_sample()))
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented variable definitions.");

    // wheter augmented variable are generic or status specific
    if (augm.is_any())
        if (augm.is_sample())
        {
            if (augm.is_generic())
                move_augm_indiv_ = &MCMC::Impl::do_move_augm_indiv_generic_;
            else
                move_augm_indiv_ = &MCMC::Impl::do_move_augm_indiv_by_status_;
        }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC::Impl : build data  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign data after pre-processing. */
void MCMC::Impl::ImplDataSet::build(const dict_str_vec &input)
{
    try
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
                tools::warning(
                    "Unknown household structure: missing data column \"house_id\".");

            // safe load of households
            for (size_t i = 0; i < n_individuals_; ++i)
            {
                // read individual data
                std::string id = std::to_string(i + 1);                // default id
                std::string house_id = std::to_string(n_houses() + 1); // default house_id
                int status = MISSING_STATUS;
                data_key_map<double> data = {};

                for (const std::string &key : headers_)
                    if (key == "id")
                    {
                        id = input.at("id")[i];
                        if (tools::is_na(id))
                            throw std::invalid_argument(
                                "Missing id for individual in row n." + std::to_string(i + 2) + ".");
                    }
                    else if (key == "house_id")
                    {
                        house_id = input.at("house_id")[i];
                        if (tools::is_na(house_id))
                            throw std::invalid_argument(
                                "Missing house_id for individual in row n." + std::to_string(i + 2) + ".");
                    }
                    else if (key == "status")
                    {
                        // note: NA status is a missing status
                        if (!tools::is_na(input.at("status")[i]))
                            status = tools::stoi(input.at("status")[i], "individual status");
                    }
                    else if (tools::is_na(input.at(key)[i]))
                    {
                        data.emplace(key, HSP_NA_DOUBLE);
                    }
                    else
                    {
                        data.emplace(key, tools::stod(input.at(key)[i], "data \"" + key + "\""));
                    }

                // if household index doesn't exist, add new household
                if (hhid_to_index_.find(house_id) == hhid_to_index_.end())
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
                house.assign_member(i, id, status, data);
            }
            headers_ = headers_ - str_vec({"id", "house_id", "status"});

            // update report for data
            fill_report_data();
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing dataset:\n" +
                     std::string(e.what()));
    }
}

/* Use the loaded data to fill-in the respective report. */
void MCMC::Impl::ImplDataSet::fill_report_data()
{
    report_indiv_ = {{"total", n_individuals_},
                     {"missing", 0},
                     {"healthy", 0},
                     {"infected", 0}};

    report_house_ = {{"total", n_houses()},
                     {"no cases", 0},
                     {"at least 1 case", 0},
                     {"at least 2 cases", 0}};

    report_size_ = {{"max cases", {0, 0}},
                    {"min size", {1000, 0}},
                    {"max size", {0, 0}}};

    //===== Set Contact Rates and Perform Database Analysis =====//
    for (const Household &house : house_vector_)
    {
        const ushort house_size = house.size();
        const ushort n_infec = house.n_infec_status();

        report_indiv_["missing"] += house.n_missing_status();
        report_indiv_["healthy"] += house.n_healthy_status();
        report_indiv_["infected"] += n_infec;

        report_house_["no cases"] += (n_infec == 0);
        report_house_["at least 1 case"] += (n_infec > 0);
        report_house_["at least 2 cases"] += (n_infec > 1);

        // max number of house cases
        if (report_size_["max cases"].first < n_infec)
        {
            report_size_["max cases"].first = n_infec;
            report_size_["max cases"].second = 1;
        }
        else if (report_size_["max cases"].first == n_infec)
            ++report_size_["max cases"].second;

        // minimum house number of members
        if (report_size_["min size"].first > house_size)
        {
            report_size_["min size"].first = house_size;
            report_size_["min size"].second = 1;
        }
        else if (report_size_["min size"].first == house_size)
            ++report_size_["min size"].second;

        // maximum house number of members
        if (report_size_["max size"].first < house_size)
        {
            report_size_["max size"].first = house_size;
            report_size_["max size"].second = 1;
        }
        else if (report_size_["max size"].first == house_size)
            ++report_size_["max size"].second;
    }
}

/* Return a value from individaul data report. */
ulong MCMC::Impl::ImplDataSet::get_report_indiv(const std::string &key) const
{
    if (report_indiv_.count(key))
        return report_indiv_.at(key);
    else
        throw std::invalid_argument(
            "Unknown individual report entry \"" + key + "\".");
    return 0.;
}

/* Return a value from household data report. */
ulong MCMC::Impl::ImplDataSet::get_report_house(const std::string &key) const
{
    if (report_house_.count(key))
        return report_house_.at(key);
    else
        throw std::invalid_argument(
            "Unknown household report entry \"" + key + "\".");
    return 0;
}

/* Return a value from household sizes report. */
std::pair<ushort, ulong> MCMC::Impl::ImplDataSet::get_report_hh_size(const std::string &key) const
{
    if (report_size_.count(key))
        return report_size_.at(key);
    else
        throw std::invalid_argument(
            "Unknown household sizes report entry \"" + key + "\".");
    return std::pair<ushort, ulong>({0, 0});
}

/* Extra integer vectors for all households. */
void MCMC::Impl::ImplDataSet::assign_cache_integer_vector(const cache_key &e)
{
    for (Household &house : restricted())
        house.assign_cache_integer_vector(e);
}
/* Extra double vectors for all households. */
void MCMC::Impl::ImplDataSet::assign_cache_double_vector(const cache_key &e)
{
    for (Household &house : restricted())
        house.assign_cache_double_vector(e);
}
/* Extra double square matrix for all households. */
void MCMC::Impl::ImplDataSet::assign_cache_double_matrix(const cache_key &e)
{
    for (Household &house : restricted())
        house.assign_cache_double_matrix(e);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : build multiplicities  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Build multiplicities. */
void MCMC::Impl::ImplDataSet::build_multiplicity(const dict_str_vec &table)
{
    try
    {
        if (pImpl_.augm.is_any() || (table.size() == 0))
        {
            // note: multiplicities are not compatible with data augmentation
            if (!table.empty())
                tools::warning(
                    "Data augmentation is incompatible with multiplicity simplification.");

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
            if (n_houses() != table.at("house_id").size())
                throw std::out_of_range(
                    "Multiplicity definitions have different dimensions than dataset.");

            restrict_indexes_.clear();
            for (size_t r = 0; r < n_houses(); ++r)
            {
                size_t H = hhid_to_index_.at(table.at("house_id")[r]);
                size_t H0 = hhid_to_index_.at(table.at("class_id")[r]);
                ulong m = tools::stoul(table.at("multiplicity")[r]);

                house_vector_[H].set_multiplicity(m, H0);

                if (house_vector_[H].is_multiplicity_head_class())
                    restrict_indexes_.emplace_back(H);
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing multiplicities:\n" +
                     std::string(e.what()));
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC::Impl : build data list  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign data after pre-processing. */
void MCMC::Impl::ImplDataList::build(const dict_str_vec &input)
{
    try
    {
        if (!input.empty())
        {
            // data sequence headers
            headers_.clear();
            for (const auto &x : input)
                headers_.push_back(x.first);
            const size_t n_rows = input.at(headers_[0]).size();

            // check for individual IDs
            if (!(headers_ % "id"))
                throw std::invalid_argument(
                    "Missing individual identifier column \"id\" in datalist.");
            // remove "id"
            headers_ = headers_ - str_vec({"id"});

            // individual map for id -> index
            // note: this could be moved elsewhere if required
            std::unordered_map<std::string, size_t> id_to_index = {};
            size_t c = 0;
            for (const Household &house : pImpl_.dataset)
                for (const Individual &member : house)
                    if (id_to_index.find(member.id()) != id_to_index.end())
                        throw std::invalid_argument(
                            "Multiple individuals id equal to \"" + member.id() + "\"");
                    else
                        id_to_index.emplace(member.id(), c++);

            // datalist allocation
            if (indiv_datalist_.empty())
                indiv_datalist_.resize(pImpl_.dataset.n_individuals());

            for (data_key_map<dbl_vec> &this_map : indiv_datalist_)
                for (const std::string &key : headers_)
                {
                    if (this_map.count(key) > 0)
                        this_map[key].clear();
                    this_map.emplace(key, dbl_vec({}));
                }

            // load individual data list
            for (size_t r = 0; r < n_rows; ++r)
            {
                double x;
                std::string id = input.at("id")[r];
                if (tools::is_na(id))
                    throw std::invalid_argument(
                        "Invalid id for datalist in row n." + std::to_string(r + 2) + ".");
                else if (id_to_index.find(id) == id_to_index.end())
                {
                    if (pImpl_.dataset.is_empty())
                        throw std::invalid_argument(
                            "Datalist definition before dataset is not supported.");
                    else
                        throw std::invalid_argument(
                            "Id \"" + id + "\" from datalist has no match in loaded dataset.");
                }

                data_key_map<dbl_vec> &this_dict = indiv_datalist_.at(id_to_index.at(id));
                for (const std::string &key : headers_)
                {
                    x = HSP_NA_DOUBLE;
                    if (!tools::is_na(input.at(key)[r]))
                        x = tools::stod(input.at(key)[r], "datalist \"" + key + "\"");

                    // append (since vector is initialized to empty)
                    this_dict[key].push_back(x);
                }
            }

            // assign pointers to relative datalist to all individuals
            c = 0;
            for (Household &house : pImpl_.dataset)
                for (Individual &member : house)
                    member.assign_pointer_to_datalist(&indiv_datalist_.at(c++));
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing datalist:\n" + std::string(e.what()));
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : build group rates  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert to a dictionary of parameter pointers and pass the pointers to the individual */
void MCMC::Impl::ImplGroupRates::assign_rate_as_ptrs_(Individual &member,
                                                      const rate_key_mapvec<par_key> &these_rates)
{
    rate_key_mapvec<const double *> rate_ptrs = {};

    for (const auto &pair : these_rates)
    {
        const rate_key &r = pair.first;

        if (str_vec({CONTACT_RATE, MISSING_RATE}) % r)
            throw std::invalid_argument(
                "The group rate name \"" + r + "\" is reserved. Choose a different one.");

        else if (!tools::is_na(r))
        {
            bool is_all_reference = true;                 // all reference state
            std::vector<const double *> vec_par_ptr = {}; // parameter pointers

            for (const par_key &p : pair.second)
            {
                if ((tools::is_na(p)) || (p == MISSING_RATE))
                    std::invalid_argument(
                        "Missing rate \"" + r + "\" for individual \"" + member.id() + "\".");

                else if (p == REFERENCE)
                    vec_par_ptr.emplace_back(&REFERENCE_VALUE);

                else
                {
                    is_all_reference = false;
                    vec_par_ptr.emplace_back(pImpl_.param.get_ptr(p));
                }
            }

            // simplify for the combined reference values
            if (is_all_reference)
                rate_ptrs[r] = {&REFERENCE_VALUE};
            else
                rate_ptrs[r] = vec_par_ptr;
        }
    }

    member.assign_all_group_rates(rate_ptrs);
}

/* Update model report and set group rates associate to new individual. */
void MCMC::Impl::ImplGroupRates::build(const dict_str_vec &deftable)
{
    if (!deftable.empty())
    {
        // check that tempalte headers are present
        str_vec data_headers = tools::headers_from_deftable(deftable);
        tools::check_for_missing_entries(data_headers, file_headers_);

        // remove template headers to keep only dataset
        data_headers = data_headers - file_headers_;
        // check that required dataset columns are present
        tools::check_for_unused_entries(data_headers, pImpl_.dataset.headers());

        const size_t n_rows = deftable.at("group_rate").size();

        // store the names of the used group rates
        headers_ = {};
        for (const rate_key &g : deftable.at("group_rate"))
            if (!(headers_ % g))
                headers_.push_back(g);

        // loop over all household members
        for (Household &house : pImpl_.dataset)
            for (Individual &member : house)
            {
                rate_key_mapvec<par_key> these_rates = {};

                // obtain group rates for given member
                for (size_t r = 0; r < n_rows; ++r)
                {
                    // check if all conditions are satisfied
                    bool tag = true;
                    for (const data_key &d : data_headers)
                        if (!tools::is_na(deftable.at(d)[r]))
                            tag &= (member.data(d) == tools::stod(deftable.at(d)[r]));

                    if (tag)
                    {
                        const rate_key &g = deftable.at("group_rate")[r];
                        const par_key &p = deftable.at("param_name")[r];
                        if (these_rates.count(g) == 0)
                            these_rates.emplace(g, str_vec({p}));
                        else
                            these_rates[g].push_back(p);
                    }
                }

                // assign rates to the individual
                assign_rate_as_ptrs_(member, these_rates);
            }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Assign parameter definitions to dictionary. */
void MCMC::Impl::ImplParameter::build(const dict_str_vec &deftable)
{
    try
    {
        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
            tools::check_for_missing_entries(headers, file_headers);
            tools::check_for_unused_entries(headers, file_headers);

            // drop current list of parameters
            dictionary.clear();
            list_const.clear();
            list_sample.clear();
            list_sample_standard.clear();
            list_sample_cache.clear();

            //===== Load parameters from the row with options =====//
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &name = deftable.at("name")[k]; // parameter name
                double value = HSP_NA_DOUBLE;                     // init value
                par_key p = name;

                // Add entry to shared parameter dictionary -----//
                if (dictionary.count(p))
                {
                    // parameter is already defined
                    throw std::invalid_argument(
                        "Parameter \"" + p + "\" has multiple definitions.");
                }
                else if (tools::stobool(deftable.at("const")[k]))
                {
                    // assign constant parameter to the set of variables
                    if (tools::is_na(deftable.at("value")[k]))
                        throw std::invalid_argument(
                            "Undefined value for constant parameter \"" + name + "\".");
                    else
                        value = tools::stod(deftable.at("value")[k], "parameter");

                    // assign constant parameter to the set of variables
                    dictionary.emplace(p, Variable(name, value));
                    // update the lists of constant parameters
                    list_const.push_back(p);
                }
                else
                {
                    // load value
                    if (!tools::is_na(deftable.at("value")[k]))
                        value = tools::stod(deftable.at("value")[k], "parameter");

                    // load boolean flag for discrete variables
                    bool is_discrete = false;
                    if (!tools::is_na(deftable.at("discr")[k]))
                        is_discrete = tools::stobool(deftable.at("discr")[k]);

                    // load prior name
                    const std::string &prior_name = deftable.at("prior")[k];

                    // load rule name
                    std::string rule_name = "draw";
                    if (!tools::is_na(deftable.at("rule")[k]))
                        rule_name = deftable.at("rule")[k];

                    // load rule scale
                    double rule_scale = 1.;
                    if (!tools::is_na(deftable.at("scale")[k]))
                        rule_scale = tools::stod(deftable.at("scale")[k]);

                    //----- Create a new Parameter Variable -----//
                    dictionary.emplace(p, Variable(name, value, is_discrete));
                    Variable &new_param = dictionary.at(p);

                    // assign a prior distribution to the parameter
                    new_param.assign_distr(prior_name);

                    // assign scalar arguments to prior distribution
                    str_vec arg_key = tools::parse_list_args(deftable.at("arg_key")[k]);
                    str_vec arg_val = tools::parse_list_args(deftable.at("arg_val")[k]);
                    for (ushort c = 0; c < arg_key.size(); ++c)
                        // UPDATE: is this check really necessary?
                        if (!tools::is_na(arg_key[c]))
                        {
                            const std::string &x = arg_val[c];
                            new_param.assign_distr_attr(
                                arg_key[c], HSP_ATTR_ARGUMENT,
                                (tools::is_na(x) ? HSP_NA_DOUBLE : tools::stod(x)));
                        }

                    // load parameter range
                    double lower = HSP_NA_DOUBLE, upper = HSP_NA_DOUBLE;
                    if (!tools::is_na(deftable.at("min")[k]))
                        lower = tools::stod(deftable.at("min")[k], "parameter minimum");
                    if (!tools::is_na(deftable.at("max")[k]))
                        upper = tools::stod(deftable.at("max")[k], "parameter maximum");

                    // set the range of values to the parameter
                    if (!tools::is_na(lower))
                        new_param.assign_distr_attr("min", HSP_ATTR_RANGE, lower);
                    if (!tools::is_na(upper))
                        new_param.assign_distr_attr("max", HSP_ATTR_RANGE, upper);
                    new_param.init_distr_reshaping();

                    // assign the rule for the proposer to the parameter
                    new_param.assign_proposer(rule_name, rule_scale);

                    // update the lists of parameters to be sampled
                    list_sample.push_back(p);
                }
            }

            // initialize all parameters to standard
            list_sample_standard = list_sample;
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing parameters:\n" +
                     std::string(e.what()));
    }
}

/* Adjust lists of parameters to be inferred making parameter as dependent on cache. */
void MCMC::Impl::ImplParameter::make_cache_dependent(const par_key &p)
{
    // remove a parameter from standard if it was defined as such
    if (!(list_sample_standard % p))
        list_sample_standard = list_sample_standard - str_vec({p});

    // add the parameter as cache
    if (!(list_sample_cache % p))
        list_sample_cache.push_back(p);
}

/* Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::ImplParameter::init()
{
    for (const par_key &p : list_sample)
    {
        Variable &this_par = dictionary.at(p);
        // user initialization (for multiple chains)
        this_par.user_reset();

        if (pImpl_.setup.rand_init_params)
        {
            if (tools::is_na(this_par.show_distr_name()))
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                this_par.rand(pImpl_.setup.prng());
        }
        else if (tools::is_na(this_par.get()))
            throw std::invalid_argument(
                "Value of parameter \"" + p + "\" is undefined with no random initialization.");
    }
}
