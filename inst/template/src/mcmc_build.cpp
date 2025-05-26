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

/* Default seed of random number generator.  */
std::int32_t _default_setup_seed()
{
    auto time_raw = std::chrono::steady_clock::now().time_since_epoch().count();
    // ensure time_now is 64bit
    std::uint64_t time_now = static_cast<std::uint64_t>(time_raw);
    // xor between upper and lower part of time stamp
    std::uint64_t mixed = time_now ^ (time_now >> 32);
    // note: constrain to 32bit for compatibility with R (mask sign bit)
    return static_cast<std::int32_t>(mixed & 0x7FFFFFFF);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC::Impl : constructor  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::Impl::Impl(std::function<double(Household &)> callback_house_loglike,
                 std::function<void(Household &, indiv)> callback_update_cache)
    : _infec_profile(*this),
      _uniform_distr(0., 1.),
      _callback_house_loglike(std::move(callback_house_loglike)),
      _callback_update_cache(std::move(callback_update_cache))
{
    setup_seed = _default_setup_seed();
}
MCMC::Impl::~Impl() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (getter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Get the parameter value. */
double MCMC::Impl::get_param(const parID &p) const
{
    return param_map.at(p).get();
}

/* Get pointer to the parameter. */
const double *MCMC::Impl::get_param_pointer(const parID &p) const
{
    return param_map.at(p).get_pointer();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (setter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the pseudo random number generator. */
rnd_gen &MCMC::Impl::get_prng() { return (this->*_do_get_prng)(); }

/* Draw a random variable uniformly distributed. */
double MCMC::Impl::uniform() { return _uniform_distr(get_prng()); }

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
        tools::error({"While building setup definitions:",
                      std::string(e.what())});
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

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC::Impl : methods (parameter moves)  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set accept rates to zero. */
void MCMC::Impl::init_move_counts()
{
    for (const parID &p : param_standard_infer_list + param_cache_infer_list)
        param_map.at(p).reset_accept_rate();
}

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
bool MCMC::Impl::propose_valid_param_move(const parID &p)
{
    Variable &this_par = param_map.at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(_prng);

    // recover buffer if not valid move
    if (!this_par.valid_move())
    {
        this_par.recover_value();
        return false;
    }

    return true;
}

/* Restore the parameter value in the buffer. */
void MCMC::Impl::recover_parameter(const parID &p, bool cache_depen)
{
    // recover buffer value of parameter
    param_map.at(p).recover_value();

    // recover buffer value of full loglike
    full_loglike = full_buffer;

    // recover buffered loglikelihood in households
    for (Household &house : house_vector)
    {
        house.recover_loglike_value();
        if (cache_depen)
            house.recover_all_cache();
    }

    // FIXME: this should not be recoverd if constant for example...
    if (cache_depen)
        recover_infec_profile();
}

/* Update the acceptance rate of the paramter. */
void MCMC::Impl::update_param_accept_rate(const parID &p, bool is_accepted)
{
    param_map.at(p).update_accept_rate(is_accepted);
}

/* Adapt the parameter Metropolis move. */
void MCMC::Impl::adapt_param_rule(const parID &p, double ratio)
{
    param_map.at(p).adapt_rule(ratio);
}

/* Metropolis-Hastings to accept or not current parameter move. */
bool MCMC::Impl::accept_valid_param_move(const parID &p)
{
    const double D_ll = full_loglike - full_buffer;
    const double move_logratio = param_map.at(p).get_move_logratio();
    return uniform() < std::exp(setup_inv_T * D_ll + move_logratio);
}

/* Getter for full loglikelihood. */
double MCMC::Impl::get_full_loglike() const
{
    return full_loglike;
}

/* Set new full loglikelihood value and buffer old value. */
void MCMC::Impl::set_full_loglike(double new_ll)
{
    full_buffer = full_loglike;
    full_loglike = new_ll;
}

/* Increment full loglikelohood by a delta, buffer old value. */
void MCMC::Impl::increment_full_loglike(double delta)
{
    set_full_loglike(full_loglike + delta);
}

/* Compute statistics from sample of full loglikelihood. */
const std::pair<double, double> MCMC::Impl::smpl_full_loglike_stats() const
{
    return tools::compute_stats(smpl_full_loglike);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC::Impl : methods (diagnostics) ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init posterior sampling. */
void MCMC::Impl::init_diagnostics()
{
    ulong sample_size = setup_c_sample_size * setup_n_chains;

    // diagnostics over full loglikelihood
    smpl_count = 0;
    smpl_full_loglike.clear();
    smpl_full_loglike = dbl_vec(sample_size);

    // initialize  other diagnostic vectors
    // UPDATE: I might want to store also power posterior sample
    if (!is_power_post())
    {
        // diagnostics over parameters
        for (const parID &p : param_standard_infer_list + param_cache_infer_list)
            smpl_param_value[p] = dbl_vec(sample_size, HOUSEPI_NA_DOUBLE);

        // diagnostics over households
        // UPDATE: init a=only for selected households
        for (size_t H : restrict_hh_indexes)
            house_vector[H].init_diagnostics(sample_size);

        // diagnostics over augmented variables
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_list_infer + augm_list_func)
                smpl_augm_value[a] = dbl_vec(sample_size, HOUSEPI_NA_DOUBLE);
    }
}

/* Update Bayesian diagnostics */
void MCMC::Impl::update_diagnostics()
{
    // diagnostics for full loglike
    smpl_full_loglike[smpl_count] = get_full_loglike();

    if (!is_power_post())
    {
        // diagnostics for all parameters
        for (const parID &p : param_standard_infer_list + param_cache_infer_list)
            smpl_param_value.at(p)[smpl_count] = get_param(p);

        // diagnostics for selected households
        for (size_t H : restrict_hh_indexes)
            house_vector[H].update_diagnostics(smpl_count);

        // diagnostics for individual followup
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_list_infer + augm_list_func)
                smpl_augm_value.at(a)[smpl_count] = house_vector[smpl_H_star].augm(smpl_i_star, a);
    }

    ++smpl_count;
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
            "Impossible to proceed. No database was found.");
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

/* Tell if it's the case of power posteriors. */
bool MCMC::Impl::is_power_post() const { return setup_n_rungs > 0; }

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
void MCMC::Impl::build_dataset(const dict_str_vec &datatable)
{
    try
    {
        if (!datatable.empty())
        {
            // drop all previous data and prepare for new
            data_headers.clear();
            for (const auto &x : datatable)
                data_headers.push_back(x.first);

            n_individuals = datatable.at(data_headers[0]).size();
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
                        if (tools::is_na(datatable.at("id")[i]))
                            throw std::invalid_argument(
                                "Missing id for individual in row n." + std::to_string(i + 2) + ".");
                        else
                            id = datatable.at("id")[i];
                    }
                    else if (key == "house_id")
                    {
                        if (tools::is_na(datatable.at("house_id")[i]))
                            throw std::invalid_argument(
                                "Missing house_id for individual in row n." + std::to_string(i + 2) + ".");
                        else
                            house_id = datatable.at("house_id")[i];
                    }
                    else if (key == "status")
                    {
                        // note: NA status is a missing status
                        if (!tools::is_na(datatable.at("status")[i]))
                            status = tools::stoi(datatable.at("status")[i], "individual status");
                    }
                    else if (tools::is_na(datatable.at(key)[i]))
                    {
                        data.emplace(key, HOUSEPI_NA_DOUBLE);
                        tools::warning("Found NA value(s) in the database.");
                    }
                    else
                        data.emplace(key, tools::stod(datatable.at(key)[i], "data \"" + key + "\""));

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
        tools::error({"While initializing dataset:",
                      std::string(e.what())});
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
        tools::error({"While initializing multiplicities:",
                      std::string(e.what())});
    }
}
