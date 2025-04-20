/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_core.cpp
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

MCMC::Impl::Impl() : _infec_profile(*this), _uniform_distr(0., 1.)
{
    setup_seed = _default_setup_seed();
}
MCMC::Impl::~Impl() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (getter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Get the parameter value. */
double MCMC::Impl::get_param(const parID &p) const { return param_map.at(p).get(); }

/* Get pointer to the parameter. */
const double *MCMC::Impl::get_param_pointer(const parID &p) const { return param_map.at(p).get_pointer(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : methods (setter)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the pseudo random number generator. */
rnd_gen &MCMC::Impl::get_prng() { return (this->*_func_get_prng)(); }

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
    _func_get_prng = &Impl::_prng_unique_thread;

#ifdef _OPENMP
    // seeds for the case of parallelization
    if (setup_n_threads > 1)
    {
        // create different generators for different parallel threads
        // note: they are used only to parallelize data augmentantion
        _func_get_prng = &Impl::_prng_parallel_threads;

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

/* Assign parameter definitions to dictionary. */
void MCMC::Impl::build_param()
{
    try
    {
        // reference to parameter definition table
        const dict_str_vec &deftable = _param_deftable;

        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
            tools::check_for_missing_entries(headers, param_file_headers);
            tools::check_for_unused_entries(headers, param_file_headers);

            // drop current list of parameters
            param_map.clear();
            param_const_list.clear();
            param_standard_infer_list.clear();
            param_extra_infer_list.clear();

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
                    param_standard_infer_list.push_back(p);
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error({"While building parameter definitions:",
                      std::string(e.what())});
    }
}

/* Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::init_parameters()
{
    for (const parID &p : param_standard_infer_list + param_extra_infer_list)
    {
        // restore user initialization (for multiple chains)
        param_map.at(p).user_reset();

        if (setup_rand_init_params)
        {
            if (tools::is_na(param_map.at(p).show_distr_name()))
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                param_map.at(p).rand(_prng);
        }
        else if (tools::is_na(param_map.at(p).get()))
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
                                        const std::string &min_at,
                                        const std::string &max_value,
                                        const std::string &max_at)
{
    //-------------------//
    //----- Minimum -----//
    //-------------------//

    if (tools::is_na(min_value) && tools::is_na(min_at))
        ; // pass

    else if (!tools::is_na(min_value) && tools::is_na(min_at))
        throw std::invalid_argument(
            "Impossible to interpret default value of  \"" + min_value + "\". Please provide a type.");

    else if (tools::is_na(min_value) && !tools::is_na(min_at))
        tools::warning(
            "Found a type without a paired value. Falling back to default.");

    else if (tools::stoupper(min_at) == TYPE_DATA)
        house.assign_augm_model_scalar_range_min(i, a, house.data(i, min_value));

    else if (tools::stoupper(min_at) == TYPE_PARAM)
    {
        if (param_map.at(min_value).is_constant())
            house.assign_augm_model_scalar_range_min(i, a, get_param(min_value));
        else
            house.assign_augm_model_pointer_range_min(i, a, get_param_pointer(min_value));
    }

    else if (tools::stoupper(min_at) == TYPE_SCALAR)
        house.assign_augm_model_scalar_range_min(i, a,
                                                 tools::stod(min_value, "scalar type argument"));

    else
        throw std::invalid_argument(
            "Unsopported type \"" + min_at + "\" found while loading \"" + min_value + "\" for pointer value.");

    //-------------------//
    //----- Maximum -----//
    //-------------------//

    if (tools::is_na(max_value) && tools::is_na(max_at))
        ; // pass

    else if (!tools::is_na(max_value) && tools::is_na(max_at))
        throw std::invalid_argument(
            "Impossible to interpret default value of  \"" + max_value + "\". Please provide a type.");

    else if (tools::is_na(max_value) && !tools::is_na(max_at))
        tools::warning(
            "Found a type without a paired value. Falling back to default.");

    else if (tools::stoupper(max_at) == TYPE_DATA)
        house.assign_augm_model_scalar_range_max(i, a, house.data(i, max_value));

    else if (tools::stoupper(max_at) == TYPE_PARAM)
    {
        if (param_map.at(max_value).is_constant())
            house.assign_augm_model_scalar_range_max(i, a, get_param(max_value));
        else
            house.assign_augm_model_pointer_range_max(i, a, get_param_pointer(max_value));
    }

    else if (tools::stoupper(max_at) == TYPE_SCALAR)
        house.assign_augm_model_scalar_range_max(i, a,
                                                 tools::stod(max_value, "argument with scalar type"));

    else
        throw std::invalid_argument(
            "Unsopported type \"" + max_at + "\" found while loading \"" + max_value + "\" for pointer value.");

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
const dict_str_vec MCMC::Impl::_build_augm_name_status_check()
{
    const dict_str_vec &deftable = _augm_deftable;
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
void MCMC::Impl::_build_augm_allocate_univocal_defs(const dict_str_vec &check)
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
        _func_augm_model_by_definition = &MCMC::Impl::__augm_model_generic;
    else
        // status specific definitions
        _func_augm_model_by_definition = &MCMC::Impl::__augm_model_by_status;
}

/* Allocate memory for augmented variables to each individual.
Search individual id for augmented traces follow-up. */
void MCMC::Impl::_build_augm_clear_houses_and_set_id_trace()
{
    // cleaer all augmented definitions
    for (Household &house : house_vector)
        for (indiv i = 0; i < house.size(); ++i)
            clear_all_augm(house, i);

    if (is_augm_sample())
    {
        if (!tools::is_na(setup_id_trace))
        {
            smpl_H_star = n_households;
            for (Household &house : house_vector)
                for (indiv i = 0; i < house.size(); ++i)
                    if (!tools::is_na(setup_id_trace))
                        if (house.id(i) == setup_id_trace)
                        {
                            smpl_H_star = house.index();
                            smpl_i_star = i;
                        }

            if (smpl_H_star == n_households)
            {
                tools::warning(
                    "Ignoring id trace: no individual \"" + setup_id_trace + "\" was found.");
                setup_id_trace = NULL_STRING;
            }
            else if (house_vector[smpl_H_star].status(smpl_i_star) < INFECTED_STATUS)
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
        const dict_str_vec &deftable = _augm_deftable;
        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(_augm_deftable);
            tools::check_for_missing_entries(headers, augm_file_headers);
            tools::check_for_unused_entries(headers, augm_file_headers);

            // auxiliary variable to load definitions : a -> status lists
            const dict_str_vec name_status_check = _build_augm_name_status_check();
            _build_augm_allocate_univocal_defs(name_status_check);

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
            _build_augm_clear_houses_and_set_id_trace();
        }
    }
    catch (const std::exception &e)
    {
        tools::error({"While building augmented variable definitions:",
                      std::string(e.what())});
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
bool MCMC::Impl::is_augm_generic() const { return _func_augm_model_by_definition == &MCMC::Impl::__augm_model_generic; }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC::Impl : assign augmented variables  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to single individual. */
void MCMC::Impl::assign_all_augm(Household &house, indiv i)
{
    try
    {
        const dict_str_vec &deftable = _augm_deftable;

        for (size_t k = 0; k < deftable.at("name").size(); ++k)
        {
            // assign augmented variable of the corresponding status
            const std::string &augm_status = deftable.at("status")[k];

            // assign augmented variable of the corresponding status or undefined status
            if ((tools::is_na(augm_status)) || (house.status(i) == tools::stoi(augm_status, "augmented status")))
            {
                // index
                // |====> UPDATE: to be implemented <===|
                augmID a = deftable.at("name")[k];
                const std::string &name = deftable.at("name")[k];

                // initial value
                double value = type_to_double(deftable.at("value")[k],
                                              deftable.at("value_at")[k],
                                              house, i);

                // constant flag
                bool is_constant = tools::stobool(deftable.at("constant")[k]);
                // discretization flag
                bool is_discrete = tools::stobool(deftable.at("discrete")[k]);
                // model name
                const std::string &model_name = deftable.at("model")[k];

                // FIXME: what if constant = T, discrete = T and value != int?

                // assign augmented variable
                house.assign_augm(i, a, name, value,
                                  is_discrete, model_name);

                // model arguments
                assign_augm_type_args(house, i, a,
                                      tools::parse_list_args(deftable.at("arg_name")[k]),
                                      tools::parse_list_args(deftable.at("arg_value")[k]),
                                      tools::parse_list_args(deftable.at("arg_at")[k]));

                // validity range
                assign_augm_type_range(house, i, a,
                                       deftable.at("min")[k],
                                       deftable.at("min_at")[k],
                                       deftable.at("max")[k],
                                       deftable.at("max_at")[k]);

                // non constant definitions
                if (!is_constant)
                {
                    // assign proposer rule
                    const std::string &rule_name = deftable.at("rule")[k];
                    double rule_scale = type_to_double(deftable.at("scale")[k],
                                                       deftable.at("scale_at")[k],
                                                       house, i);
                    double rule_from = type_to_double(deftable.at("from")[k],
                                                      deftable.at("from_at")[k],
                                                      house, i);
                    house.assign_augm_proposer(i, a,
                                               rule_name, rule_scale, rule_from);

                    // random initialization
                    // UPDATE: I would like a more general switch for this
                    if (tools::stobool(deftable.at("init_rand")[k]))
                        house.augm_model_rand(i, name, get_prng());
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error({"While initializing augmented variables:",
                      std::string(e.what())});
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
    for (const parID &p : param_standard_infer_list + param_extra_infer_list)
        param_map.at(p).reset_accept_rate();
}

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
bool MCMC::Impl::propose_valid_param_move(const parID &p)
{
    // set new value of the parameter and store the old to buffer
    param_map.at(p).propose(_prng);

    // restore buffer if not valid move
    if (!param_map.at(p).valid_move())
    {
        param_map.at(p).restore_value();
        return false;
    }

    return true;
}

/* Restore the parameter value in the buffer. */
void MCMC::Impl::restore_parameter(const parID &p, bool extra_depen)
{
    // restore buffer value of parameter
    param_map.at(p).restore_value();

    // restore buffer value of full loglike
    full_loglike = full_buffer;

    // restore buffered loglikelihood in households
    for (Household &house : house_vector)
    {
        house.restore_loglike();
        if (extra_depen)
            house.restore_all_extra();
    }

    // FIXME: this should not be restored if constant for example...
    if (extra_depen)
        restore_infec_profile();
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
    return uniform() < std::exp(setup_inv_T * D_ll + param_map.at(p).get_move_logratio());
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
    smpl_full_loglike = dbl_vec(sample_size);

    // initialize  other diagnostic vectors
    // UPDATE: I might want to store also power posterior sample
    if (!is_power_post())
    {
        // diagnostics over parameters
        for (const parID &p : param_standard_infer_list + param_extra_infer_list)
            smpl_param_value[p] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());

        // diagnostics over households
        // UPDATE: init a=only for selected households
        for (size_t H : restrict_hh_indexes)
            house_vector[H].init_diagnostics(sample_size);

        // diagnostics over augmented variables
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_infer_list)
                smpl_augm_value[a] = dbl_vec(sample_size, std::numeric_limits<double>::quiet_NaN());
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
        for (const parID &p : param_standard_infer_list + param_extra_infer_list)
            smpl_param_value.at(p)[smpl_count] = get_param(p);

        // diagnostics for selected households
        for (size_t H : restrict_hh_indexes)
            house_vector[H].update_diagnostics(smpl_count);

        // diagnostics for individual followup
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_infer_list)
                smpl_augm_value.at(a)[smpl_count] = house_vector[smpl_H_star].augm(smpl_i_star, a);
    }

    ++smpl_count;
}

/* Compute the loglikelihood term associated to augmented variables. */
double MCMC::Impl::compute_augm_model_loglike(const Household &house) const
{
    double hh_loglike = 0.;

    // Contribution of model associated to augmented variables
    for (const indiv &i : house.infected_status_list())
        (this->*_func_augm_model_by_definition)(house, i, hh_loglike);

    return hh_loglike;
}

void MCMC::Impl::__augm_model_by_status(const Household &house, indiv i, double &delta) const
{
    for (const augmID &a : augm_status_map_model.at(house.status(i)))
        delta += house.augm_model_logpdf(i, a);
}

void MCMC::Impl::__augm_model_generic(const Household &house, indiv i, double &delta) const
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
    if (((param_standard_infer_list + param_extra_infer_list).empty()) && (!is_augm_sample()))
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
            if (std::find(data_headers.begin(), data_headers.end(), "house_id") == data_headers.end())
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
                house.assign_individual(id, status, data);
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
