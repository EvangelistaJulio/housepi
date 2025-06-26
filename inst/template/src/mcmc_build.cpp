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

void MCMC::Impl::Setup::build()
{
    if (deflist.size() > 0)
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

/* Assign setup. */
void MCMC::Impl::build_setup()
{
    try
    {
        setup.build();
    }
    catch (const std::exception &e)
    {
        tools::error("While building setup definitions:\n" +
                     std::string(e.what()));
    }
}

/* Set the number of threads for parrell computations (default is 1). */
void MCMC::Impl::Setup::set_n_threads()
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
void MCMC::Impl::Setup::set_prng_seed()
{
    // seed for the unique thread rand num generator
    _prng.seed(seed);
    _do_get_prng = &Impl::Setup::_prng_unique_thread;

#ifdef _OPENMP
    // seeds for the case of parallelization
    if (n_threads > 1)
    {
        // create different generators for different parallel threads
        // note: they are used only to parallelize data augmentantion
        _do_get_prng = &Impl::Setup::_prng_parallel_threads;

        std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

        _prng_parallel.reserve(n_threads);
        for (unsigned n = 0; n < n_threads; ++n)
            _prng_parallel[n].seed(seed_generator(_prng));
    }
#endif
}

void MCMC::Impl::Setup::assert_mcmc_setup()
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
const std::string MCMC::Impl::Setup::fmt_progtxt()
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
void MCMC::Impl::assert_proper_dataset()
{
    if (house_vector.size() == 0)
        throw std::length_error(
            "Unable to proceed. No database was found.");
}

/* Assert setup definitions are ok to proceed with mcmc sampling. */
void MCMC::Impl::assert_proper_mcmc_setup()
{
    setup.assert_mcmc_setup();

    // check if there is a dataset
    if (n_individuals == 0)
        throw std::length_error(
            "Dataset is empty or was not provided.");

    // check if there is actually something to sample
    if (((param_standard_infer_list + param_cache_infer_list).empty()) && (!is_augm_sample()))
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented variable definitions.");
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
void MCMC::Impl::build_dataset(const dict_str_vec &input)
{
    try
    {
        if (input.size() > 0)
        {
            // drop all previous data and prepare for new
            str_vec data_headers = {};
            for (const auto &x : input)
                data_headers.push_back(x.first);

            n_individuals = input.at(data_headers[0]).size();

            house_vector.clear();
            restrict_hh_indexes.clear();
            _hhid_to_index.clear();

            // assess the presence of house IDs
            if (!contains(data_headers, "house_id"))
                tools::warning(
                    "Unknown household structure: missing data column \"house_id\".");

            // safe load of households
            for (size_t i = 0; i < n_individuals; ++i)
            {
                // read individual data
                std::string id = std::to_string(i + 1);                // default id
                std::string house_id = std::to_string(n_houses() + 1); // default house_id
                int status = MISSING_STATUS;
                data_key_map<double> data = {};

                for (const std::string &key : data_headers)
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
                        tools::warning("Found NA values in the database.");
                    }
                    else
                        data.emplace(key, tools::stod(input.at(key)[i], "data \"" + key + "\""));

                // if household index doesn't exist, add new household
                if (_hhid_to_index.find(house_id) == _hhid_to_index.end())
                {
                    // map data index to household number
                    _hhid_to_index.emplace(house_id, n_houses());
                    // default choice for restricted indeces
                    restrict_hh_indexes.push_back(n_houses());
                    // create a new empty household
                    house_vector.emplace_back(house_id, n_houses());
                }

                // extract the household reference
                Household &house = house_vector.at(_hhid_to_index.at(house_id));

                // assign the individual the respective household
                house.assign_member(i, id, status, data);
            }

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
void MCMC::Impl::fill_report_data()
{
    _report_indiv_ = {{"total", n_individuals},
                      {"missing", 0},
                      {"healthy", 0},
                      {"infected", 0}};

    _report_house_ = {{"total", n_houses()},
                      {"no cases", 0},
                      {"at least 1 case", 0},
                      {"at least 2 cases", 0}};

    _report_hh_size_ = {{"max cases", {0, 0}},
                        {"min size", {1000, 0}},
                        {"max size", {0, 0}}};

    //===== Set Contact Rates and Perform Database Analysis =====//
    for (const Household &house : house_vector)
    {
        const ushort house_size = house.size();
        const ushort n_infec = house.n_infec_status();

        _report_indiv_["missing"] += house.n_missing_status();
        _report_indiv_["healthy"] += house.n_healthy_status();
        _report_indiv_["infected"] += n_infec;

        _report_house_["no cases"] += (n_infec == 0);
        _report_house_["at least 1 case"] += (n_infec > 0);
        _report_house_["at least 2 cases"] += (n_infec > 1);

        // max number of house cases
        if (_report_hh_size_["max cases"].first < n_infec)
        {
            _report_hh_size_["max cases"].first = n_infec;
            _report_hh_size_["max cases"].second = 1;
        }
        else if (_report_hh_size_["max cases"].first == n_infec)
            ++_report_hh_size_["max cases"].second;

        // minimum house number of members
        if (_report_hh_size_["min size"].first > house_size)
        {
            _report_hh_size_["min size"].first = house_size;
            _report_hh_size_["min size"].second = 1;
        }
        else if (_report_hh_size_["min size"].first == house_size)
            ++_report_hh_size_["min size"].second;

        // maximum house number of members
        if (_report_hh_size_["max size"].first < house_size)
        {
            _report_hh_size_["max size"].first = house_size;
            _report_hh_size_["max size"].second = 1;
        }
        else if (_report_hh_size_["max size"].first == house_size)
            ++_report_hh_size_["max size"].second;
    }
}

/* Return a value from individaul data report. */
ulong MCMC::Impl::get_report_indiv(const std::string &key) const
{
    if (_report_indiv_.count(key))
        return _report_indiv_.at(key);
    else
        throw std::invalid_argument(
            "Unknown individual report entry \"" + key + "\".");
    return 0.;
}

/* Return a value from household data report. */
ulong MCMC::Impl::get_report_house(const std::string &key) const
{
    if (_report_house_.count(key))
        return _report_house_.at(key);
    else
        throw std::invalid_argument(
            "Unknown household report entry \"" + key + "\".");
    return 0;
}

/* Return a value from household sizes report. */
std::pair<ushort, ulong> MCMC::Impl::get_report_hh_size(const std::string &key) const
{
    if (_report_hh_size_.count(key))
        return _report_hh_size_.at(key);
    else
        throw std::invalid_argument(
            "Unknown household sizes report entry \"" + key + "\".");
    return std::pair<ushort, ulong>({0, 0});
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
void MCMC::Impl::build_multiplicity(const dict_str_vec &table)
{
    try
    {
        if (is_augm_any() || (table.size() == 0))
        {
            // note: multiplicities are not compatible with data augmentation
            if (table.size() > 0)
                tools::warning(
                    "Data augmentation is incompatible with multiplicity simplification.");

            // default choice for restricted household indexes
            restrict_hh_indexes.clear();
            restrict_hh_indexes.resize(n_houses());
            for (size_t H = 0; H < n_houses(); ++H)
            {
                house_vector[H].set_multiplicity(1, H);
                restrict_hh_indexes[H] = H;
            }
        }
        else
        {
            if (n_houses() != table.at("house_id").size())
                throw std::out_of_range(
                    "Multiplicity definitions have different dimensions than dataset.");

            restrict_hh_indexes.clear();
            for (size_t r = 0; r < n_houses(); ++r)
            {
                size_t H = _hhid_to_index.at(table.at("house_id")[r]);
                size_t H0 = _hhid_to_index.at(table.at("class_id")[r]);
                ulong m = tools::stoul(table.at("multiplicity")[r]);

                house_vector[H].set_multiplicity(m, H0);

                if (house_vector[H].is_multiplicity_head_class())
                    restrict_hh_indexes.emplace_back(H);
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
void MCMC::Impl::build_datalist(const dict_str_vec &input)
{
    try
    {
        if (input.size() > 0)
        {
            // data sequence headers
            str_vec headers = {};
            for (const auto &x : input)
                headers.push_back(x.first);
            const size_t n_rows = input.at(headers[0]).size();

            // check for individual IDs
            if (!contains(headers, "id"))
                throw std::invalid_argument(
                    "Missing individual identifier column \"id\" in datalist.");
            // remove "id"
            headers.erase(std::remove(headers.begin(), headers.end(), "id"), headers.end());

            // individual map for id -> index
            // note: this could be moved elsewhere if required
            std::unordered_map<std::string, size_t> id_to_index = {};
            size_t c = 0;
            for (Household &house : house_vector)
                for (Individual &member : house)
                    if (id_to_index.find(member.id()) != id_to_index.end())
                        throw std::invalid_argument(
                            "Multiple individuals id equal to \"" + member.id() + "\"");
                    else
                        id_to_index.emplace(member.id(), c++);

            // datalist allocation
            if (indiv_datalist.empty())
                indiv_datalist.resize(n_individuals);

            for (data_key_map<dbl_vec> &this_map : indiv_datalist)
                for (const std::string &key : headers)
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
                    if (house_vector.empty())
                        throw std::invalid_argument(
                            "Datalist definition before dataset is not supported.");
                    else
                        throw std::invalid_argument(
                            "Id \"" + id + "\" from datalist has no match in loaded dataset.");
                }

                data_key_map<dbl_vec> &this_dict = indiv_datalist.at(id_to_index.at(id));
                for (const std::string &key : headers)
                {
                    x = HSP_NA_DOUBLE;
                    if (tools::is_na(input.at(key)[r]))
                        tools::warning("Found NA values in datalist.");
                    else
                        x = tools::stod(input.at(key)[r], "datalist \"" + key + "\"");

                    // append (since vector is initialized to empty)
                    this_dict[key].push_back(x);
                }
            }

            // assign pointers to relative datalist to all individuals
            c = 0;
            for (Household &house : house_vector)
                for (Individual &member : house)
                    member.assign_pointer_to_datalist(&indiv_datalist.at(c++));
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

/* Assign parameter definitions to dictionary. */
void MCMC::Impl::build_param()
{
    try
    {
        // reference to parameter definition table
        const dict_str_vec &deftable = _param_deftable;

        if (deftable.size() > 0)
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
            tools::check_for_missing_entries(headers, param_file_headers);
            tools::check_for_unused_entries(headers, param_file_headers);

            // drop current list of parameters
            param_map.clear();
            param_const_list.clear();
            param_standard_infer_list.clear();
            param_cache_infer_list.clear();

            //===== Load parameters from the row with options =====//
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &name = deftable.at("name")[k]; // parameter name
                double value = HSP_NA_DOUBLE;                     // init value
                par_key p = name;

                // Add entry to shared parameter dictionary -----//
                if (param_map.count(p))
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
                    param_map.emplace(p, Variable(name, value));
                    // update the lists of constant parameters
                    param_const_list.push_back(p);
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
                    param_map.emplace(p, Variable(name, value, is_discrete));
                    Variable &new_param = param_map.at(p);

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

                    // update the lists of parameters to be inferred
                    param_standard_infer_list.push_back(p);
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While building parameters:\n" +
                     std::string(e.what()));
    }
}

/* Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::init_parameters()
{
    for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
    {
        Variable &this_par = param_map.at(p);
        // user initialization (for multiple chains)
        this_par.user_reset();

        if (setup.rand_init_params)
        {
            if (tools::is_na(this_par.show_distr_name()))
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                this_par.rand(setup.prng());
        }
        else if (tools::is_na(this_par.get()))
            throw std::invalid_argument(
                "Value of parameter \"" + p + "\" is undefined with no random initialization.");
    }
}
