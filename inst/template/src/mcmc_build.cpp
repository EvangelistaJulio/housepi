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

#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (setter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign setup. */
void MCMC::Impl::build_setup()
{
    try
    {
        // reference to definition list
        const dict_str_str &deflist = _setup_deflist;

        if (!deflist.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deflist(_setup_deflist);
            tools::check_for_unused_entries(headers, setup_file_headers);

            // Random Initialization of Parameters
            if (deflist.count("progress") && (!tools::is_na(deflist.at("progress"))))
                setup_progress = tools::stoul(deflist.at("progress"));

            // N. of Threads
            if (deflist.count("n_threads") && (!tools::is_na(deflist.at("n_threads"))))
                setup_n_threads = tools::stoi(deflist.at("n_threads"), "n. of threads");

            // PRNG Seed
            if (deflist.count("seed") && (!tools::is_na(deflist.at("seed"))))
                setup_seed = tools::stoi32(deflist.at("seed"), "random seed");

            // N. of Chains
            if (deflist.count("n_chains") && (!tools::is_na(deflist.at("n_chains"))))
                setup_n_chains = tools::stoul(deflist.at("n_chains"), "n. of chains");

            // Single Chain Sample Size
            if (deflist.count("sample_size") && (!tools::is_na(deflist.at("sample_size"))))
                setup_c_sample_size = tools::stoul(deflist.at("sample_size"), "sample size");

            // N. of Steps
            if (deflist.count("n_steps") && (!tools::is_na(deflist.at("n_steps"))))
                setup_n_steps = tools::stoul(deflist.at("n_steps"), "n. of steps");

            // N. of Adaptive Steps
            if (deflist.count("n_adapt") && (!tools::is_na(deflist.at("n_adapt"))))
                setup_n_adapt = tools::stoul(deflist.at("n_adapt"), "n. of adaptive steps");

            // N. of Thinning Steps
            if (deflist.count("thin") && (!tools::is_na(deflist.at("thin"))))
                setup_thin = tools::stoul(deflist.at("thin"), "thin");

            // N. of Steps to "Burn-In"
            if (deflist.count("n_burnin") && (!tools::is_na(deflist.at("n_burnin"))))
                setup_n_burnin = tools::stoul(deflist.at("n_burnin"), "burn-in");

            // Random Initialization of Parameters
            if (deflist.count("rand_init") && (!tools::is_na(deflist.at("rand_init"))))
                setup_rand_init_params = tools::stobool(deflist.at("rand_init"));

            // Inverse Temperature
            if (deflist.count("inv_temper") && (!tools::is_na(deflist.at("inv_temper"))))
                setup_inv_T = tools::stod(deflist.at("inv_temper"), "inverse temperature");

            // Individual ID for Augmented Trace
            if (deflist.count("id_trace") && (!tools::is_na(deflist.at("id_trace"))))
                setup_id_trace = deflist.at("id_trace");
        }

        // define number of threads
        set_n_threads();

        // pseudo-random number generator seeding
        set_prng_seed();

        // reset internal states
        chain = 1;         // current chain
        rung = 0;          // current thermodynamic rung
        setup_n_rungs = 0; // number of rungs
    }
    catch (const std::exception &e)
    {
        tools::error("While building setup definitions:" +
                     std::string(e.what()));
    }
}

/* Set the number of threads for parrell computations (default is 1). */
void MCMC::Impl::set_n_threads()
{
    if (setup_n_threads == 0)
        throw std::invalid_argument("N. of threads must be greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_num_procs();

    if (setup_n_threads > max_threads)
    {
        tools::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available processors.");
        setup_n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(setup_n_threads);

#else
    if (setup_n_threads != 1)
        tools::warning("Parallelization with OpenMP not supported.");

    setup_n_threads = 1;
#endif
}

/* Assign seed for random generators. */
void MCMC::Impl::set_prng_seed()
{
    // seed for the unique thread rand num generator
    _prng.seed(setup_seed);
    _do_get_prng = &Impl::_prng_unique_thread;

#ifdef _OPENMP
    // seeds for the case of parallelization
    if (setup_n_threads > 1)
    {
        // create different generators for different parallel threads
        // note: they are used only to parallelize data augmentantion
        _do_get_prng = &Impl::_prng_parallel_threads;

        std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

        _prng_parallel.reserve(setup_n_threads);
        for (unsigned n = 0; n < setup_n_threads; ++n)
            _prng_parallel[n].seed(seed_generator(_prng));
    }
#endif
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
    if (house_vector.empty())
        throw std::length_error(
            "Unable to proceed. No database was found.");
}

/* Assert setup definitions are ok to proceed with mcmc sampling. */
void MCMC::Impl::assert_proper_mcmc_setup()
{
    if ((setup_inv_T < 0.) || (setup_inv_T > 1.))
        throw std::range_error(
            "Inverse temperature " + std::to_string(setup_inv_T) + " out of range [0, 1].");

    if ((setup_n_chains > 1) && (setup_rand_init_params == false))
    {
        tools::warning(
            "Unable to fix initial parameters when running multiple chains.");
        setup_rand_init_params = true;
    }
    else if (setup_n_chains < 1)
        throw std::out_of_range("N. of chains must be greater than 0.");

    if (setup_thin == 0)
    {
        tools::warning("Thinning equal to 0 interpreted as 1.");
        setup_thin = 1;
    }

    if (setup_c_sample_size == 0)
    {
        if (setup_n_steps < 1)
            throw std::out_of_range("N. of steps must be greater than 0.");

        if (setup_n_burnin >= setup_n_steps)
            throw std::out_of_range("Burn-in must be smaller than n. of steps.");

        if (setup_thin >= setup_n_steps)
            throw std::out_of_range("Thinning must be smaller than n. of steps. ");

        setup_c_sample_size = (setup_n_steps - setup_n_burnin) / setup_thin;

        if (setup_c_sample_size == 0)
            throw std::out_of_range("Choice of steps, burn-in and thin return a 0 sample size.");
    }
    else
        setup_n_steps = setup_n_burnin + setup_c_sample_size * setup_thin;

    // assign number of sweeps
    setup_c_n_sweeps = setup_n_adapt / setup_thin;
    setup_c_n_sweeps += setup_n_burnin / setup_thin;
    setup_c_n_sweeps += setup_c_sample_size;

    // check if there is a dataset
    if (n_individuals == 0)
        throw std::length_error("Dataset is empty or was not provided.");

    // check if there is actually something to sample
    if (((param_standard_infer_list + param_cache_infer_list).empty()) && (!is_augm_sample()))
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented variable definitions.");
}

/* The text to print in the Markov chian progress bar. */
const std::string MCMC::Impl::fmt_progtxt()
{
    if (is_power_post())
        return "power post " + std::to_string(rung) + "/" + std::to_string(setup_n_rungs);
    else if (setup_n_chains > 1)
        return "sample " + std::to_string(chain) + "/" + std::to_string(setup_n_chains);
    else
        return "sample";
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
        if (!input.empty())
        {
            // drop all previous data and prepare for new
            data_headers.clear();
            for (const auto &x : input)
                data_headers.push_back(x.first);

            n_individuals = input.at(data_headers[0]).size();
            n_households = 0;

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
                std::string id = std::to_string(i + 1);                  // default id
                std::string house_id = std::to_string(n_households + 1); // default house_id
                int status = MISSING_STATUS;
                dataIDmap<double> data = {};

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
                        data.emplace(key, HOUSEPI_NA_DOUBLE);
                        tools::warning("Found NA values in the database.");
                    }
                    else
                        data.emplace(key, tools::stod(input.at(key)[i], "data \"" + key + "\""));

                // if household index doesn't exist, add new household
                if (_hhid_to_index.find(house_id) == _hhid_to_index.end())
                {
                    // map data index to household number
                    _hhid_to_index.emplace(house_id, n_households);
                    // default choice for restricted indeces
                    restrict_hh_indexes.push_back(n_households);
                    // create a new empty household
                    house_vector.emplace_back(house_id, n_households++);
                }

                // extract the household reference
                Household &house = house_vector.at(_hhid_to_index.at(house_id));

                // assign the individual the respective household
                house.assign_member(i, id, status, data);
            }

            // update report for data
            fill_report_data();

            // clean up data entries for output files
            // eventually remove default headers to reorder them at the end
            const str_vec def_data_headers = {"house_id", "id", "status"};
            for (const std::string &col : def_data_headers)
                data_headers.erase(
                    std::remove(data_headers.begin(), data_headers.end(), col),
                    data_headers.end());
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing dataset:" +
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

    _report_house_ = {{"total", n_households},
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
        if (is_augm_any() || table.empty())
        {
            // note: multiplicities are not compatible with data augmentation
            if (!table.empty())
                tools::warning(
                    "Data augmentation is incompatible with multiplicity simplification.");

            // default choice for restricted household indexes
            restrict_hh_indexes.clear();
            restrict_hh_indexes.resize(n_households);
            for (size_t H = 0; H < n_households; ++H)
            {
                house_vector[H].set_multiplicity(1, H);
                restrict_hh_indexes[H] = H;
            }
        }
        else
        {
            if (n_households != table.at("house_id").size())
                throw std::out_of_range(
                    "Multiplicity definitions have different dimensions than dataset.");

            restrict_hh_indexes.clear();
            for (size_t r = 0; r < n_households; ++r)
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
        tools::error("While initializing multiplicities:" +
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
void MCMC::Impl::build_datalist(const dict_str_vec &input, const std::string &list_name)
{
    // alias for messages
    std::string alias = "data list";
    if (!tools::is_na(list_name))
        alias += " \"" + list_name + "\"";

    try
    {
        if (input.empty())
        {
            throw std::invalid_argument(
                "Nothing to build: " + alias + " is empty.");
        }
        else
        {
            // data sequence headers
            str_vec headers = {};
            for (const auto &x : input)
                headers.push_back(x.first);
            const size_t n_rows = input.at(headers[0]).size();

            // check for individual IDs
            if (!contains(headers, "id"))
                throw std::invalid_argument(
                    "Missing individual identifier column \"id\" in " + alias + ".");
            // remove "id"
            headers.erase(std::remove(headers.begin(), headers.end(), "id"), headers.end());

            // individual map for id -> index
            // note: this could be moved elsewhere if required
            std::unordered_map<std::string, size_t> id_to_index = {};
            for (Household &house : house_vector)
                for (Individual &member : house)
                    if (id_to_index.find(member.id()) != id_to_index.end())
                        throw std::invalid_argument(
                            "Multiple individuals id equal to \"" + member.id() + "\"");
                    else
                        id_to_index.emplace(member.id(), member.index());

            // datalist allocation
            if (indiv_datalist.empty())
                indiv_datalist.resize(n_individuals);

            for (dataIDmap<dbl_vec> &this_map : indiv_datalist)
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
                        "Invalid id for " + alias + " in row n." + std::to_string(r + 2) + ".");
                else if (id_to_index.find(id) == id_to_index.end())
                    throw std::invalid_argument(
                        "Id " + id + " from " + alias + " unmatched in loaded dataset.");

                dataIDmap<dbl_vec> &this_dict = indiv_datalist.at(id_to_index.at(id));
                for (const std::string &key : data_headers)
                {
                    x = HOUSEPI_NA_DOUBLE;
                    if (tools::is_na(input.at(key)[r]))
                        tools::warning("Found NA values in " + alias + ".");
                    else
                        x = tools::stod(input.at(key)[r], alias + " \"" + key + "\"");

                    // append (since vector is initialized to empty)
                    this_dict[key].push_back(x);
                }
            }
        }

        // assign pointers to relative datalist to all individuals
        for (Household &house : house_vector)
            for (Individual &member : house)
                member.assign_pointer_to_datalist(&indiv_datalist.at(member.index()));
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing " + alias + ":" +
                     std::string(e.what()));
    }
}
