/*
 * mcmc_core.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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
/////////////////////////////////  MCMC::Impl : methods (getter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Get the parameter value. */
double MCMC::Impl::get_param(const parID &p) const { return param_map.at(p).get_value(); }

/* Get pointer to the parameter. */
const double *MCMC::Impl::get_param_pointer(const parID &p) const { return param_map.at(p).get_pointer(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (setter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign seed for random generators. */
void MCMC::Impl::set_prng_seed(int seed)
{
    // seed for the unique thread rand num generator
    prng.seed(seed);

    // seeds for the case of parallelization
    if (setup_n_threads > 1)
    {
        std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);
        for (unsigned n = 0; n < setup_n_threads; ++n)
            prng_parallel[n].seed(seed_generator(prng));
    }
}

/* Return the pseudo random number generator. */
rnd_gen &MCMC::Impl::get_prng() { return (this->*_func_prng_)(); }

/* Draw a random variable uniformly distributed. */
double MCMC::Impl::uniform() { return _uniform_distr_(get_prng()); }

/* Set the number of threads for parrell computations (default is 1). */
void MCMC::Impl::set_n_threads()
{
    if (setup_n_threads == 0)
        throw std::invalid_argument("N. of threads must be greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_max_threads();
    if (setup_n_threads > max_threads)
    {
        tools::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available threads.");
        setup_n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(setup_n_threads);

    // create different generators for different parallel threads
    // note: they are used only to parallelize data augmentantion
    if (setup_n_threads > 1)
    {
        prng_parallel.reserve(setup_n_threads);
        _func_prng_ = &Impl::_prng_parallel_threads_;
    }
    else
        _func_prng_ = &Impl::_prng_unique_thread_;

#else
    if (setup_n_threads != 1)
        tools::warning("Parallelization with OpenMP not supported.");

    setup_n_threads = 1;
    _func_prng_ = &Impl::_prng_unique_thread_;
#endif
}

/* Set inverse temperature and its relative functional methods. */
void MCMC::Impl::set_inver_temper(double inv_T)
{
    if ((inv_T < 0.) || (inv_T > 1.))
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " out of range [0, 1].");
    else
        setup_inv_T = inv_T;
}

/* Assign setup. */
void MCMC::Impl::build_setup()
{
    try
    {
        //===== DEFAULT SETUP VALUES =====//
        auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();

        setup_seed = static_cast<unsigned long long>(time_now); // seed of random number generator
        setup_rand_init_params = true;                          // random initialization of parameters

        setup_n_chains = 1;    // n. chains
        setup_n_steps = 10000; // n. iterations
        setup_n_adapt = 0;     // n. (finite) adaptive steps
        setup_thin = 1;        // n. thinning iterations
        setup_n_burnin = 0;    // n. burnin iterations

        setup_c_sample_size = 0; // single chain sample size
        setup_c_n_sweeps = 0;    // single chain n. of sweeps

        setup_inv_T = 1.0;   // inverse temperature
        setup_n_threads = 1; // number of OMP threads

        setup_id_trace = NULL_STRING; // augmented trace to save

        // reference to definition list
        const dict_str_str &deflist = _setup_deflist_;

        if (!deflist.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deflist(_setup_deflist_);
            tools::check_for_unused_entries(headers, setup_file_headers);

            // N. of Threads
            if (deflist.count("n_threads") && (!tools::is_na(deflist.at("n_threads"))))
                setup_n_threads = tools::stoi(deflist.at("n_threads"), "n. of threads");

            // PRNG Seed
            if (deflist.count("seed") && (!tools::is_na(deflist.at("seed"))))
                setup_seed = tools::stoull(deflist.at("seed"), "random seed");

            // N. of Chains
            if (deflist.count("n_chains") && (!tools::is_na(deflist.at("n_chains"))))
                setup_n_chains = tools::stoul(deflist.at("n_chains"), "n. of chains");

            // Single Chain Sample Size
            if (deflist.count("size") && (!tools::is_na(deflist.at("size"))))
                setup_c_sample_size = tools::stoul(deflist.at("size"), "sample size");

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
                set_inver_temper(tools::stod(deflist.at("inv_temper"), "inverse temperature"));

            // Individual ID for Augmented Trace
            if (deflist.count("id_trace") && (!tools::is_na(deflist.at("id_trace"))))
                setup_id_trace = deflist.at("id_trace");
        }

        // define number of threads
        set_n_threads();

        // pseudo-random number generator seeding
        set_prng_seed(setup_seed);

        // reset internal states
        _chain_ = 1;   // current chain
        _rung_ = 0;    // current thermodynamic rung
        _n_rungs_ = 0; // number of rungs
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While building setup definitions: " + std::string(e.what()));
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Assign parameter definitions to dictionary. */
void MCMC::Impl::build_params()
{
    try
    {
        // reference to parameter definition table
        const dict_str_vec &deftable = _param_deftable_;

        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
            tools::check_for_missing_entries(headers, param_file_headers);
            tools::check_for_unused_entries(headers, param_file_headers);

            // drop current list of parameters
            param_map.clear();
            param_const_list.clear();
            param_infer_list.clear();

            //===== Load parameters from the row with options =====//
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &name = deftable.at("name")[k]; // parameter name
                double value = NULL_DOUBLE;                       // init value

                // UPDATE : temporary parameter identifier
                // NOTE : id parID changes this should change too
                parID p = name;

                //----- Add entry to shared parameter dictionary -----//
                if (param_map.count(p))
                {
                    // parameter is already defined
                    throw std::invalid_argument(
                        "Parameter \"" + p + "\" has multiple definitions.");
                }
                else if (tools::stobool(deftable.at("constant")[k]))
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
                    if (!tools::is_na(deftable.at("discrete")[k]))
                        is_discrete = tools::stobool(deftable.at("discrete")[k]);

                    // load prior name
                    const std::string &prior_name = deftable.at("prior")[k];

                    // load range
                    double lower = NULL_DOUBLE;
                    if (!tools::is_na(deftable.at("min")[k]))
                        lower = tools::stod(deftable.at("min")[k], "parameter minimum");

                    double upper = NULL_DOUBLE;
                    if (!tools::is_na(deftable.at("max")[k]))
                        upper = tools::stod(deftable.at("max")[k], "parameter maximum");

                    // load rule name
                    std::string rule_name = "gibbs";
                    if (!tools::is_na(deftable.at("rule")[k]))
                        rule_name = deftable.at("rule")[k];

                    // load rule scale
                    double rule_scale = 1.;
                    if (!tools::is_na(deftable.at("rule_scale")[k]))
                        rule_scale = tools::stod(deftable.at("rule_scale")[k]);

                    //----- Create a new Parameter Variable -----//
                    param_map.emplace(p, Variable(name, value, is_discrete));
                    Variable &new_param = param_map.at(p);

                    // assign a prior distribution to the parameter
                    new_param.assign_distribution(prior_name);

                    // assign scalar arguments to prior distribution
                    str_vec arg_name = tools::parse_list_args(deftable.at("arg_name")[k]);
                    str_vec arg_value = tools::parse_list_args(deftable.at("arg_value")[k]);
                    for (ushort c = 0; c < arg_name.size(); ++c)
                        if (!tools::is_na(arg_name[c]))
                            new_param.assign_distr_scalar_arg(arg_name[c],
                                                              tools::natod(arg_value[c]));

                    // set the range of values to the parameter
                    if (!tools::is_na(lower))
                        new_param.assign_distr_scalar_range_min(lower);
                    if (!tools::is_na(upper))
                        new_param.assign_distr_scalar_range_max(upper);
                    new_param.init_distr_reshaping();

                    // assign the rule for the proposer to the parameter
                    new_param.assign_proposer(rule_name, rule_scale);

                    // update the lists of parameters to be inferred
                    param_infer_list.push_back(p);
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While building parameter definitions: " + std::string(e.what()));
    }
}

/* Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::init_parameters()
{
    for (const parID &p : param_infer_list)
    {
        // restore user initialization (meaningful for multiple chains)
        param_map.at(p).user_reset();

        if (setup_rand_init_params)
        {
            if (tools::is_na(param_map.at(p).show_distr_name()))
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                param_map.at(p).rand(prng);
        }
        else if (tools::is_na(param_map.at(p).get_value()))
            throw std::invalid_argument(
                "Value of parameter \"" + p + "\" is undefined with no random initialization.");
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC::Impl : methods (type-defined variables)  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed std::string to a double */
double MCMC::Impl::type_to_double(const std::string &value,
                                  const std::string &type,
                                  const Household &house, indiv i)
{
    double output = NULL_DOUBLE;

    try
    {
        if (!(tools::is_na(value) && tools::is_na(type)))
        {
            if (!tools::is_na(value) && tools::is_na(type))
                output = tools::stod(value, "argument without type");

            else if (tools::is_na(value) && !tools::is_na(type))
                tools::warning(
                    "Ignoring value type without a paired value.");

            else if (tools::stoupper(type) == TYPE_SCALAR)
                output = tools::stod(value, "argument with scalar type");

            else if (tools::stoupper(type) == TYPE_DATA)
                output = house.data(i, value);

            else if (tools::stoupper(type) == TYPE_PARAM)
                // UPDATE: here it won't be straightforward to drop maps
                if (param_map.at(value).is_constant())
                    output = get_param(value);
                else
                    throw std::invalid_argument(
                        "Parameter to be inferred \"" + value + "\" assigned as constant argument.");

            else
                throw std::invalid_argument(
                    "Unknown variable type specifier \"" + type + "\".");
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return output;
}

/* Convert a std::string with type-dictionary to std::string-to-double */
void MCMC::Impl::assign_augm_type_args(Household &house, indiv i, augmID a,
                                       const str_vec &keys,
                                       const str_vec &vals,
                                       const str_vec &typs)
{
    // check size validity of keys, values and type
    if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
        throw std::invalid_argument(
            "Mismatching definitions in arguments.");

    // loop over all arguments
    for (size_t x = 0; x < keys.size(); ++x)
    {
        std::string key = keys[x];
        std::string value = vals[x];
        std::string type = typs[x];

        if (tools::is_na(value) && tools::is_na(type))
            ; // pass

        else if (!tools::is_na(value) && tools::is_na(type))
            throw std::invalid_argument(
                "Impossible to interpret default value of  \"" + value + "\". Please provide a type.");

        else if (tools::is_na(value) && !tools::is_na(type))
            tools::warning(
                "Found a type without a paired value. Falling back to default.");

        else if (tools::stoupper(type) == TYPE_DATA)
            house.assign_augm_model_scalar_arg(i, a, key, house.data(i, value));

        else if (tools::stoupper(type) == TYPE_PARAM)
        {
            if (param_map.at(value).is_constant())
                house.assign_augm_model_scalar_arg(i, a, key, get_param(value));
            else
                house.assign_augm_model_pointer_arg(i, a, key, get_param_pointer(value));
        }
        else if (tools::stoupper(type) == TYPE_SCALAR)
            house.assign_augm_model_scalar_arg(i, a, key, tools::stod(value, "argument with scalar type"));

        else
            throw std::invalid_argument(
                "Unsopported type \"" + type + "\" found while loading \"" + value + "\" for pointer value.");
    }
}

/* Convert a std::string with type-dictionary to std::string-to-double */
void MCMC::Impl::assign_augm_type_range(Household &house, indiv i, augmID a,
                                        const std::string &min_value,
                                        const std::string &min_type,
                                        const std::string &max_value,
                                        const std::string &max_type)
{
    //-------------------//
    //----- Minimum -----//
    //-------------------//

    if (tools::is_na(min_value) && tools::is_na(min_type))
        ; // pass

    else if (!tools::is_na(min_value) && tools::is_na(min_type))
        throw std::invalid_argument(
            "Impossible to interpret default value of  \"" + min_value + "\". Please provide a type.");

    else if (tools::is_na(min_value) && !tools::is_na(min_type))
        tools::warning(
            "Found a type without a paired value. Falling back to default.");

    else if (tools::stoupper(min_type) == TYPE_DATA)
        house.assign_augm_model_scalar_range_min(i, a, house.data(i, min_value));

    else if (tools::stoupper(min_type) == TYPE_PARAM)
    {
        if (param_map.at(min_value).is_constant())
            house.assign_augm_model_scalar_range_min(i, a, get_param(min_value));
        else
            house.assign_augm_model_pointer_range_min(i, a, get_param_pointer(min_value));
    }

    else if (tools::stoupper(min_type) == TYPE_SCALAR)
        house.assign_augm_model_scalar_range_min(i, a,
                                                 tools::stod(min_value, "scalar type argument"));

    else
        throw std::invalid_argument(
            "Unsopported type \"" + min_type + "\" found while loading \"" + min_value + "\" for pointer value.");

    //-------------------//
    //----- Maximum -----//
    //-------------------//

    if (tools::is_na(max_value) && tools::is_na(max_type))
        ; // pass

    else if (!tools::is_na(max_value) && tools::is_na(max_type))
        throw std::invalid_argument(
            "Impossible to interpret default value of  \"" + max_value + "\". Please provide a type.");

    else if (tools::is_na(max_value) && !tools::is_na(max_type))
        tools::warning(
            "Found a type without a paired value. Falling back to default.");

    else if (tools::stoupper(max_type) == TYPE_DATA)
        house.assign_augm_model_scalar_range_max(i, a, house.data(i, max_value));

    else if (tools::stoupper(max_type) == TYPE_PARAM)
    {
        if (param_map.at(max_value).is_constant())
            house.assign_augm_model_scalar_range_max(i, a, get_param(max_value));
        else
            house.assign_augm_model_pointer_range_max(i, a, get_param_pointer(max_value));
    }

    else if (tools::stoupper(max_type) == TYPE_SCALAR)
        house.assign_augm_model_scalar_range_max(i, a,
                                                 tools::stod(max_value, "argument with scalar type"));

    else
        throw std::invalid_argument(
            "Unsopported type \"" + max_type + "\" found while loading \"" + max_value + "\" for pointer value.");

    //---------------------//
    //----- Reshaping -----//
    //---------------------//

    house.init_augm_model_reshaping(i, a);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC::Impl : build augmented  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Produce auxiliary variable to load definitions : a -> status lists */
const dict_str_vec MCMC::Impl::_build_augm_name_status_check_()
{
    const dict_str_vec &deftable = _augm_deftable_;
    const size_t table_length = deftable.at("name").size();

    dict_str_vec output = {};

    //===== Check Augmented Variables lists of usage =====//
    for (size_t k = 0; k < table_length; ++k)
    {
        const std::string &a = deftable.at("name")[k];
        const std::string &status = deftable.at("status")[k];

        // check that augmID and status are correct
        if (tools::is_na(a))
            throw std::invalid_argument(
                "Augmented variable name can not be undefined.");

        else if (output.count(a))
        {
            // check no multiple augmID with same status
            const str_vec &mylist = output.at(a);
            if (std::find(mylist.begin(), mylist.end(), status) == mylist.end())
                output.at(a).push_back(status);
            else
                throw std::invalid_argument(
                    "Incompatible definitions for augmented variable \"" + a + "\".");
        }

        else // new augmented variable status
            output.emplace(a, str_vec{status});
    }

    return output;
}

/* Check univocal definitions and allocate new lists. */
void MCMC::Impl::_build_augm_allocate_univocal_defs_(const dict_str_vec &check)
{
    // drop previous definitions
    augm_list_full.clear();
    augm_status_map_infer.clear();
    augm_status_map_model.clear();

    for (const auto &pair : check)
    {
        const std::string &a = pair.first;
        const str_vec &mylist = pair.second;

        for (const std::string &status : mylist)
        {
            if (!tools::is_na(status))
            {
                int s = tools::stoi(status, "augmented status");
                if (s >= INFECTED_STATUS)
                {
                    augm_status_map_infer.emplace(s, std::unordered_set<std::string>{});
                    augm_status_map_model.emplace(s, std::unordered_set<std::string>{});
                }
                else
                    throw std::invalid_argument(
                        "Impossible to associate augmented \"" + a + "\" for non infectious status.");
            }
            else if (mylist.size() > 1)
                throw std::invalid_argument(
                    "Conflicting status-specific and general definitions for augmented \"" + a + "\".");

            augm_list_full.emplace(a);
        }
    }

    // set generic or status-specifc augmented variables
    if (augm_status_map_infer.empty())
        // only generic definitions
        _func_augm_model_by_definition_ = &MCMC::Impl::_augm_model_generic_;
    else
        // status specific definitions
        _func_augm_model_by_definition_ = &MCMC::Impl::_augm_model_by_status_;
}

/* Allocate memory for augmented variables to each individual.
Search individual id for augmented traces follow-up. */
void MCMC::Impl::_build_augm_clear_houses_and_set_id_trace_()
{
    // cleaer all augmented definitions
    for (Household &house : house_vector)
        for (indiv i = 0; i < house.size(); ++i)
            clear_all_augm(house, i);

    if (is_augm_sample())
    {
        if (!tools::is_na(setup_id_trace))
        {
            dgn_H_star = n_households;
            for (Household &house : house_vector)
                for (indiv i = 0; i < house.size(); ++i)
                    if (!tools::is_na(setup_id_trace))
                        if (house.id(i) == setup_id_trace)
                        {
                            dgn_H_star = house.index();
                            dgn_i_star = i;
                        }

            if (dgn_H_star == n_households)
            {
                tools::warning(
                    "Ignoring id trace: no individual \"" + setup_id_trace + "\" was found.");
                setup_id_trace = NULL_STRING;
            }
            else if (house_vector[dgn_H_star].status(dgn_i_star) < INFECTED_STATUS)
                tools::warning(
                    "Individual id \"" + setup_id_trace + "\" correspond to a healthy individual.");
        }
    }
    else if (!tools::is_na(setup_id_trace))
    {
        tools::warning(
            "Ignoring id trace: no augmented variables to be inferred was provided.");
        setup_id_trace = NULL_STRING;
    }
}

/* Assign augmented variables to all individuals. */
void MCMC::Impl::build_augm()
{
    try
    {
        const dict_str_vec &deftable = _augm_deftable_;
        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(_augm_deftable_);
            tools::check_for_missing_entries(headers, augm_file_headers);
            tools::check_for_unused_entries(headers, augm_file_headers);

            // auxiliary variable to load definitions : a -> status lists
            const dict_str_vec name_status_check = _build_augm_name_status_check_();
            _build_augm_allocate_univocal_defs_(name_status_check);

            //===========================================//
            //===== Assign new definitions to lists =====//
            //===========================================//

            // drop previous definitions
            augm_infer_list.clear();
            augm_list_model.clear();

            // note : I didn't manage to collapse this with the previous because of generic status
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &a = deftable.at("name")[k];
                const std::string &status = deftable.at("status")[k];

                // read option for constant variable
                bool is_constant = false;
                if (!tools::is_na(deftable.at("constant")[k]))
                    is_constant = tools::stobool(deftable.at("constant")[k]);

                //----- Augmented Variables to be Inferred -----//
                if (!is_constant)
                {
                    // check random initialization option
                    if (!tools::is_na(deftable.at("init_rand")[k]))
                        if (tools::stobool(deftable.at("init_rand")[k]))
                            if (tools::is_na(deftable.at("model")[k]))
                                throw std::invalid_argument(
                                    "Impossible to initialize \"" + a + "\" without a model.");

                    // read option for monte carlo proposer
                    if (!tools::is_na(deftable.at("rule")[k]))
                    {
                        augm_infer_list.emplace(a);

                        if (!tools::is_na(status)) // status-specific definition
                            augm_status_map_infer.at(tools::stoi(status, "augmented status")).emplace(a);

                        else if (!is_augm_generic()) // generic definition
                            for (const auto &pair : augm_status_map_infer)
                                augm_status_map_infer.at(pair.first).emplace(a);
                    }
                }

                //----- Augmented Variables with a Model -----//
                if (!tools::is_na(deftable.at("model")[k]))
                {
                    augm_list_model.emplace(a);

                    if (!tools::is_na(status)) // status-specific definition
                        augm_status_map_model.at(tools::stoi(status)).emplace(a);

                    else if (!is_augm_generic()) // generic definition
                        for (const auto &pair : augm_status_map_model)
                            augm_status_map_model.at(pair.first).emplace(a);
                }
            }

            // finalize assigment in the households
            _build_augm_clear_houses_and_set_id_trace_();
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While building augmented variable definitions: " + std::string(e.what()));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC::Impl : state checkers   /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if there are ant augmented variable definitions. */
bool MCMC::Impl::is_augm_any() const { return !augm_list_full.empty(); }

/* Check if there are non-constant augmented variable definitions. */
bool MCMC::Impl::is_augm_sample() const { return !augm_infer_list.empty(); }

/* Check if there are augmented variable definitions with a model. */
bool MCMC::Impl::is_augm_model() const { return !augm_list_model.empty(); }

/* Check if augmented variable definition are generic or not (i. e. individual status specific). */
bool MCMC::Impl::is_augm_generic() const { return _func_augm_model_by_definition_ == &MCMC::Impl::_augm_model_generic_; }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC::Impl : assign augmented variables  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::Impl::assign_all_augm(Household &house, indiv i)
{
    try
    {
        const dict_str_vec &deftable = _augm_deftable_;

        for (size_t k = 0; k < deftable.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = deftable.at("status")[k];

            // assign augmented variable of the corresponding status or undefined status
            if ((tools::is_na(augm_status)) || (house.status(i) == tools::stoi(augm_status, "augmented status")))
            {
                //----- Index -----//
                // |====> UPDATE: to be implemented <===|
                augmID a = deftable.at("name")[k];

                //----- Name -----//
                const std::string &name = deftable.at("name")[k];

                //----- is Constant -----//
                bool is_constant = tools::stobool(deftable.at("constant")[k]);

                //----- Value -----//
                double value = type_to_double(deftable.at("value")[k],
                                              deftable.at("value_type")[k],
                                              house, i);

                //----- Model -----//
                const std::string &model_name = deftable.at("model")[k];

                //----- is Discrete -----//
                bool is_discrete = tools::stobool(deftable.at("discrete")[k]);

                if (is_constant)
                {
                    //----- Assign Constant Augmented Variable -----//
                    house.assign_augm(i, a, name, value,
                                      is_discrete, model_name);
                }
                else
                {
                    //----- Proposer Rule -----//
                    const std::string &rule_name = deftable.at("rule")[k];
                    double rule_scale = type_to_double(deftable.at("rule_scale")[k],
                                                       deftable.at("rule_scale_type")[k],
                                                       house, i);
                    double rule_from = type_to_double(deftable.at("rule_from")[k],
                                                      deftable.at("rule_from_type")[k],
                                                      house, i);

                    //----- Assign Augmented Variable -----//
                    house.assign_augm(i, a, name, value,
                                      is_discrete, model_name,
                                      rule_name, rule_scale, rule_from);
                }

                //----- Arguments -----//
                assign_augm_type_args(house, i, a,
                                      tools::parse_list_args(deftable.at("arg_name")[k]),
                                      tools::parse_list_args(deftable.at("arg_value")[k]),
                                      tools::parse_list_args(deftable.at("arg_type")[k]));

                //----- Range -----//
                assign_augm_type_range(house, i, a,
                                       deftable.at("min")[k],
                                       deftable.at("min_type")[k],
                                       deftable.at("max")[k],
                                       deftable.at("max_type")[k]);

                // random initialization at assignment
                if (!is_constant)
                    if (tools::stobool(deftable.at("init_rand")[k]))
                        house.augm_model_rand(i, name, get_prng());
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While initializing augmented variables: " + std::string(e.what()));
    }
}

/* Clear all augmented definition of individual. */
void MCMC::Impl::clear_all_augm(Household &house, indiv i)
{
    for (const augmID &a : augm_list_full)
        house.clear_augm(i, a);
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
    for (const parID &p : param_infer_list)
        param_map.at(p).reset_accept_rate();
}

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
bool MCMC::Impl::propose_valid_param_move(const parID &p)
{
    bool output = true;

    // set new value of the parameter and store the old to buffer
    param_map.at(p).propose(prng);

    // restore buffer if not valid move
    if (!param_map.at(p).valid_move())
    {
        param_map.at(p).restore_value();
        output = false;
    }

    return output;
}

/* Restore the parameter value in the buffer. */
void MCMC::Impl::restore_simple_param_value(const parID &p)
{
    // restore buffer value of parameter
    param_map.at(p).restore_value();

    // restore buffer value of global loglike
    global_loglike = buffer;

    // restore buffered loglikelihood in households
    for (Household &house : house_vector)
        house.restore_loglike();
}

/* Restore buffered parameter and its extra dependecies. */
void MCMC::Impl::restore_extra_param_value(const parID &p)
{
    // restore buffer value of parameter
    param_map.at(p).restore_value();

    // restore buffer value of global loglike
    global_loglike = buffer;

    // restore buffered loglikelihood in households
    for (Household &house : house_vector)
    {
        house.restore_loglike();
        house.restore_all_extra();
    }
}

/* Update the acceptance rate of the paramter. */
void MCMC::Impl::update_param_accept_rate(const parID &p, bool is_accepted) { param_map.at(p).update_accept_rate(is_accepted); }

/* Adapt the parameter Metropolis move. */
void MCMC::Impl::adapt_param_rule(const parID &p, double ratio) { param_map.at(p).adapt_rule(ratio); }

/* Metropolis-Hastings to accept or not current parameter move. */
bool MCMC::Impl::accept_valid_param_move(const parID &p)
{
    return uniform() < exp(setup_inv_T * (global_loglike - buffer) + param_map.at(p).get_move_logratio());
}

/* Extract average and variance from vector of global loglike values. */
void MCMC::Impl::dgn_glob_ll_stats(double &avg, double &var) const { tools::compute_stats(dgn_glob_loglike, avg, var); }

double MCMC::Impl::get_global_loglike() const { return global_loglike; };

void MCMC::Impl::increment_global_loglike(double delta) { global_loglike += delta; };

void MCMC::Impl::set_global_loglike(double new_ll)
{
    buffer = global_loglike;
    global_loglike = new_ll;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void MCMC::Impl::init_diagnostics()
{
    ulong sample_size = setup_c_sample_size * setup_n_chains;

    // diagnostics over global loglikelihood
    dgn_count = 0;
    dgn_glob_loglike = dbl_vec(sample_size);

    // initialize  other diagnostic vectors
    if (!is_power_post())
    {
        // diagnostics over parameters
        for (const parID &p : param_infer_list + param_extra_infer_list)
            dgn_param_value[p] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());

        // diagnostics over households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].init_diagnostics(sample_size);

        // diagnostics over augmented variables
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_infer_list)
                dgn_augm_value[a] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());
    }
}

/* Update Bayesian diagnostics */
void MCMC::Impl::update_diagnostics()
{
    // diagnostics for global loglike
    dgn_glob_loglike[dgn_count] = get_global_loglike();

    if (!is_power_post())
    {
        // diagnostics for all parameters
        for (const parID &p : param_infer_list + param_extra_infer_list)
            dgn_param_value.at(p)[dgn_count] = get_param(p);

        // diagnostics for all households
        for (size_t H = 0; H < n_households; ++H)
            house_vector[H].update_diagnostics(dgn_count);

        // diagnostics for individual followup
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_infer_list)
                dgn_augm_value.at(a)[dgn_count] = house_vector[dgn_H_star].augm(dgn_i_star, a);
    }

    ++dgn_count;
}

/* Compute the loglikelihood term associated to augmented variables. */
double MCMC::Impl::compute_augm_model_loglike(const Household &house) const
{
    double hh_loglike = 0.;

    // Contribution of model associated to augmented variables
    for (const indiv &i : house.infectious_list())
        (this->*_func_augm_model_by_definition_)(house, i, hh_loglike);

    return hh_loglike;
}

void MCMC::Impl::_augm_model_by_status_(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_status_map_model.at(house.status(i)))
        delta += house.augm_model_logpdf(i, a);
}

void MCMC::Impl::_augm_model_generic_(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_list_model)
        delta += house.augm_model_logpdf(i, a);
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
    if ((param_infer_list.empty()) && (!is_augm_sample()))
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented variable definitions.");
}

/* Tell if it's the case of power posteriors. */
bool MCMC::Impl::is_power_post() const { return _n_rungs_ > 0; }

/* The text to print in the Markov chian progress bar. */
const std::string MCMC::Impl::fmt_progtxt()
{
    if (is_power_post())
        return "power post " + std::to_string(_rung_) + "/" + std::to_string(_n_rungs_);
    else if (setup_n_chains > 1)
        return "sample " + std::to_string(_chain_) + "/" + std::to_string(setup_n_chains);
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
///////////////////////////////////////  MCMC : build data  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign data after pre-processing. */
void MCMC::Impl::build_data(const dict_str_vec &datatable)
{
    try
    {
        if (!datatable.empty())
        {
            //=====  Drop all Previous Data and Prepare for New  =====//
            data_headers.clear();
            for (const auto &x : datatable)
                data_headers.push_back(x.first);
            n_individuals = datatable.at(data_headers[0]).size();
            n_households = 0;
            house_vector.clear();

            //===== Assess the presence of house IDs =====//
            if (std::find(data_headers.begin(), data_headers.end(), "house_id") == data_headers.end())
                tools::warning(
                    "Unknown household structure: missing data column \"house_id\".");

            //===== Safe Load of Households =====//
            std::unordered_map<std::string, size_t> hh_ID_to_index = {};
            for (size_t i = 0; i < n_individuals; ++i)
            {
                //----- Read individual data -----//
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
                            status = tools::stoul(datatable.at("status")[i], "individual status");
                    }
                    else if (tools::is_na(datatable.at(key)[i]))
                    {
                        data.emplace(key, std::numeric_limits<double>::quiet_NaN());
                        tools::warning("Found NA value(s) in the database.");
                    }
                    else
                        data.emplace(key, tools::stod(datatable.at(key)[i], "data \"" + key + "\""));


                // if household index doesn't exist, add new household
                if (hh_ID_to_index.find(house_id) == hh_ID_to_index.end())
                {
                    // map data index to household number
                    hh_ID_to_index.emplace(house_id, n_households);
                    // create a new empty household
                    house_vector.emplace_back(house_id, n_households++);
                }

                // extract the household reference
                Household &house = house_vector[hh_ID_to_index.at(house_id)];

                // Assign the individual the respective household
                house.assign_individual(id, status, data);
            }

            //===== Update Report for Data =====//
            update_report_data();

            //===== Clean up data entries for output files =====//
            // eventually remove standard headers to reorder them at the end
            for (const std::string &col : {"house_id", "id", "status"})
                data_headers.erase(
                    std::remove(data_headers.begin(), data_headers.end(), col),
                    data_headers.end());
        }
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While initializing dataset: " + std::string(e.what()));
    }
}

/* Use the loaded data to fill-in the respective report. */
void MCMC::Impl::update_report_data()
{
    _report_data_ = {{"n.cases", 0},
                     {"n.miss.stat", 0},
                     {"no.HH.case", 0},
                     {"lst.1.case", 0},
                     {"lst.2.case", 0},
                     {"max.n.HH.cases", 0},
                     {"n.HH.with.max.cases", 0},
                     {"min.HH.size", 1000},
                     {"n.HH.with.min.size", 0},
                     {"max.HH.size", 0},
                     {"n.HH.with.max.size", 0}};

    //===== Set Contact Rates and Perform Database Analysis =====//
    for (const Household &house : house_vector)
    {
        const ushort house_size = house.size();
        const ushort n_infec = house.n_infectious();

        _report_data_["n.cases"] += n_infec;
        _report_data_["n.miss.stat"] += house.n_missing_status();

        _report_data_["no.HH.case"] += (n_infec == 0);
        _report_data_["lst.1.case"] += (n_infec > 0);
        _report_data_["lst.2.case"] += (n_infec > 1);

        // max number of house cases
        if (_report_data_["max.n.HH.cases"] < n_infec)
        {
            _report_data_["n.HH.with.max.cases"] = 1;
            _report_data_["max.n.HH.cases"] = n_infec;
        }
        else if (_report_data_["max.n.HH.cases"] == n_infec)
            ++_report_data_["n.HH.with.max.cases"];

        // minimum house number of members
        if (_report_data_["min.HH.size"] > house_size)
        {
            _report_data_["n.HH.with.min.size"] = 1;
            _report_data_["min.HH.size"] = house_size;
        }
        else if (_report_data_["min.HH.size"] == house_size)
            ++_report_data_["n.HH.with.min.size"];

        // maximum house number of members
        if (_report_data_["max.HH.size"] < house_size)
        {
            _report_data_["n.HH.with.max.size"] = 1;
            _report_data_["max.HH.size"] = house_size;
        }
        else if (_report_data_["max.HH.size"] == house_size)
            ++_report_data_["n.HH.with.max.size"];
    }
}

/* Return a value from report data. */
size_t MCMC::Impl::report_data(const std::string &key) const
{
    size_t output = 0;
    if (_report_data_.count(key))
        output = _report_data_.at(key);

    else
        throw std::invalid_argument(
            "Unknown summary entry \"" + key + "\".");

    return output;
}
