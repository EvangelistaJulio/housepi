/*
 * mcmc_core.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc_build.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : methods  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init the values for the Bayesian diagnostics. */
void Parameter::init_diagnostics(ulong sample_size)
{
    dgn_count = 0;
    dgn_value.clear();
    dgn_value = std::vector<double>(sample_size, std::numeric_limits<double>::quiet_NaN());
}

/* Update the values for the Bayesian diagnostics. */
void Parameter::update_diagnostics() { dgn_value[dgn_count++] = get_value(); }

/* Message to be used at the moment of display diagnostics. (TEMPORARY) */
std::string Parameter::show_diagnostics() const
{
    using namespace tools;
    std::ostringstream oss;

    double avg = 0, var = 0;
    tools::compute_stats(dgn_value, avg, var);

    oss << avg << C_DELIM << std::sqrt(var) << C_DELIM << truth(get_infer())
        << C_DELIM << get_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
        << C_DELIM << get_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
        << C_DELIM << (tools::is_na(show_range_min()) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(show_range_min()))
        << C_DELIM << (tools::is_na(show_range_max()) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(show_range_max()));
    return oss.str();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Build : methods (getter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Default functional method to return parameter. */
double Monte_Carlo_Build::_get_param_default_(const parID &p) const { return param_map.at(p).get_value(); }

/* Safe method to get a value with previous assertion. */
double Monte_Carlo_Build::_get_param_assert_(const parID &p) const
{
    return (assert_param(p)) ? _get_param_default_(p) : INF_DOUBLE;
}

/* Default functional method to return pointer to parameter. */
const double *Monte_Carlo_Build::_get_param_pointer_default_(const parID &p) const { return param_map.at(p).get_pointer(); }

/* Safe method to get a pointer to a parameter with previous assertion. */
const double *Monte_Carlo_Build::_get_param_pointer_assert_(const parID &p) const
{
    return (assert_param(p)) ? _get_param_pointer_default_(p) : nullptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Build : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a new parameter to the dictionary. */
void Monte_Carlo_Build::assign_param(const parID &p, const Parameter &new_param)
{
    // create the new parameter if it's not already defined
    if (param_map.find(p) != param_map.end())
        throw std::invalid_argument(
            "Parameter '" + p + "' has multiple definitions.");
    else
    {
        // emplace the new parameter
        param_map.emplace(p, new_param);

        // update the lists of parameters to be inferred, and these to be static
        if (new_param.get_infer())
            param_to_be_inferred.push_back(p);
        else
            param_to_be_static.push_back(p);
    }
}

/* Add a new scalar for augmented typedict. TODO : is there a more elegant option? */
// FIXME : this might break because of dynamic allocation
const double *Monte_Carlo_Build::assign_scalar_pointer(double val)
{
    if (_scalar_type_params_.find(val) == _scalar_type_params_.end())
    {
        Parameter new_param(UNDEFINED, val);
        _scalar_type_params_.emplace(val, std::move(new_param));
    }

    return _scalar_type_params_.at(val).get_pointer();
}

/* Set the number of threads for parrell computations (default is 1). */
void Monte_Carlo_Build::set_n_threads(unsigned n)
{
    if (n == 0)
        throw std::invalid_argument("Number of threads must be equal or greater than 1.");

#ifdef _OPENMP
    if (n > _max_n_available_threads)
    {
        std::ostringstream oss;
        oss << "Number of threads (" << n << ") changed to number of available threads (" << _max_n_available_threads << ").";
        tools::warning(oss.str());
        n = _max_n_available_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n);
    setup_n_threads = n;

    // create different generators for different parallel threads
    // UPDATE: they are used only to parallelize data augmentantion
    if (n > 1)
    {
        prng_parallel.reserve(n);
        _func_prng_ = &Monte_Carlo_Build::_prng_parallel_threads_;
    }
    else
    {
        _func_prng_ = &Monte_Carlo_Build::_prng_unique_thread_;
    }
#else
    if (n != 1)
        tools::warning("Parallelization with OpenMP not supported.");

    setup_n_threads = 1;
    _func_prng_ = &Monte_Carlo_Build::_prng_unique_thread_;
#endif
}

/* Assign seed for random generators. */
void Monte_Carlo_Build::set_prng_seed(int seed)
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

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Build : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
bool Monte_Carlo_Build::assert_param(const parID &p) const
{
    if (param_map.find(p) == param_map.end())
        throw std::invalid_argument(
            "Parameter '" + p + "' not found. Check your model and parameter definitions.");
    else
        return true;

    return false;
}

/* Wrapper of functional method to return the parameter value. */
double Monte_Carlo_Build::get_param(const parID &p) const { return (this->*_func_get_param_)(p); }

/* Wrapper of functional method to return a pointer to the parameter. */
const double *Monte_Carlo_Build::get_param_pointer(const parID &p) const { return (this->*_func_get_param_pointer_)(p); }

/* Getter for the infer attribute of parameter (mutable/constant). */
bool Monte_Carlo_Build::get_infer(const parID &p) const { return param_map.at(p).get_infer(); }

/* Randomly initialize values of the parameters to be inferred */
void Monte_Carlo_Build::init_parameters()
{
    for (const parID &p : param_to_be_inferred)
    {
        // restore user initialization (meaningful for multiple chains)
        param_map.at(p).user_reset();

        if (setup_rnd_init_params)
        {
            if (param_map.at(p).get_distr_name() == UNDEFINED)
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                param_map.at(p).rand(prng);
        }
        else
        {
            if (param_map.at(p).get_value() == INF_DOUBLE)
                throw std::invalid_argument(
                    "Undefined value for parameter to infer '" + p + "' while no random initialization.");
        }
    }
}

/* Set accept rates to zero. */
void Monte_Carlo_Build::init_move_counts()
{
    for (const parID &p : param_to_be_inferred)
        param_map.at(p).reset_accept_rate();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Monte_Carlo_Build : methods (parameter moves)  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
// FIXME : Gibbs sampling is missing here!
// 1) always valid 2) move_logratio == 0
bool Monte_Carlo_Build::propose_valid_param_move(const parID &p)
{
    bool output = true;

    // set new value of the parameter and store the old to buffer
    param_map.at(p).propose(prng);

    // restore buffer if not valid move
    if (!param_map.at(p).valid_move())
    {
        restore_param_value(p);
        output = false;
    }

    return output;
}

/* Restore the parameter value in the buffer. */
void Monte_Carlo_Build::restore_param_value(const parID &p) { param_map.at(p).restore_value(); }

/* Update the acceptance rate of the paramter. */
void Monte_Carlo_Build::update_param_accept_rate(const parID &p, bool is_accepted) { param_map.at(p).update_accept_rate(is_accepted); }

/* Metropolis logratio associate to the current move. */
double Monte_Carlo_Build::param_move_logratio(const parID &p) const { return param_map.at(p).get_move_logratio(); }

/* Adapt the parameter Metropolis move. */
void Monte_Carlo_Build::adapt_param_rule(const parID &p, double ratio) { param_map.at(p).adapt_rule(ratio); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Read Setup file and produce a definition list. */
void Monte_Carlo_Build::read_setup(const std::string &pathtofile)
{
    tools::read_csv_list(pathtofile, _setup_deflist_);
    build_setup();
    _report_pathtofile_["setup"] = pathtofile;
}

/* Assign setup. */
void Monte_Carlo_Build::build_setup()
{
    try
    {
        // default number of threads
        unsigned n_threads = 1;

        // default prng seed
        //----- default prng seed -----//
        auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();
        setup_seed = static_cast<unsigned long long>(time_now);

        if (!_setup_deflist_.empty())
        {
            //===== Check that all provided setups are correct =====//
            str_vec file_entries = {"n_threads", "seed", "n_chains",
                                    "n_steps", "n_adapt", "thin",
                                    "n_burnin", "rnd_init_params",
                                    "n_rungs", "inv_temper"};

            str_vec wrong_entries = {};

            for (const auto &pair : _setup_deflist_)
                if (find(file_entries.begin(), file_entries.end(), pair.first) == file_entries.end())
                    wrong_entries.push_back(pair.first);

            if (wrong_entries.size() > 0)
            {
                std::ostringstream oss;
                oss << "Unknown setup option(s): '" << wrong_entries[0] << "'";
                for (size_t k = 1; k < wrong_entries.size(); ++k)
                    oss << ", '" << wrong_entries[k] << "'";
                oss << ". Available setup options are: '" << file_entries[0] << "'";
                for (size_t k = 1; k < file_entries.size(); ++k)
                    oss << ", '" << file_entries[k] << "'";
                tools::warning(oss.str() + ".");
            }

            //===== ASSIGN SETUP CONFIGURATIONS =====//

            //===== NUMBER OF THREADS =====//
            if (_setup_deflist_.find("n_threads") != _setup_deflist_.end())
                n_threads = tools::stoi(_setup_deflist_.at("n_threads"), "number of threads");

            //===== PRNG SEED =====//
            if ((_setup_deflist_.find("seed") != _setup_deflist_.end()) && (!tools::is_na(_setup_deflist_.at("seed"))))
                setup_seed = tools::stoull(_setup_deflist_.at("seed"), "random seed");

            //===== NUMBER OF CHAINS =====//
            if (_setup_deflist_.find("n_chains") != _setup_deflist_.end())
            {
                setup_n_chains = tools::stoul(_setup_deflist_.at("n_chains"), "number of chains");
                if (setup_n_steps < 1)
                    throw std::out_of_range(
                        "Number of chains must be greater than 0.");
            }

            //===== NUMBER OF ITERATIONS =====//
            if (_setup_deflist_.find("n_steps") != _setup_deflist_.end())
            {
                setup_n_steps = tools::stoul(_setup_deflist_.at("n_steps"), "number of iterations");
                if (setup_n_steps < 2)
                    throw std::out_of_range(
                        "Number of iterations must be greater than 1.");
            }

            //===== NUMBER OF (FINITE) ADAPTIVE ITERATIONS =====//
            if (_setup_deflist_.find("n_adapt") != _setup_deflist_.end())
                setup_n_adapt = tools::stoul(_setup_deflist_.at("n_adapt"), "number of adaptive steps");

            //===== THINNING =====//
            if (_setup_deflist_.find("thin") != _setup_deflist_.end())
            {
                setup_thin = tools::stoul(_setup_deflist_.at("thin"), "thin");
                if (setup_thin >= setup_n_steps)
                    throw std::out_of_range(
                        "Thinning must be smaller than number of iterations. ");
                else if (setup_thin == 0)
                    throw std::out_of_range(
                        "Thinning must be greater than 0.");
            }

            //===== BURN-IN =====//
            if (_setup_deflist_.find("n_burnin") != _setup_deflist_.end())
            {
                setup_n_burnin = tools::stoul(_setup_deflist_.at("n_burnin"), "burn-in");
                if (setup_n_burnin >= setup_n_steps)
                    throw std::out_of_range(
                        "Burn-in must be smaller than number of iterations.");
            }

            //===== RANDOM INITIALIZATION OF PARAMETERS =====//
            if (_setup_deflist_.find("rnd_init_params") != _setup_deflist_.end())
            {
                setup_rnd_init_params = tools::stobool(_setup_deflist_.at("rnd_init_params"));
                if ((setup_n_chains > 1) && (setup_rnd_init_params == false))
                {
                    tools::warning(
                        "Unable to fix initial parameters when running multiple chains.");
                    setup_rnd_init_params = true;
                }
            }

            //===== INVERSE TEMPERATURE =====//
            if (_setup_deflist_.find("inv_temper") != _setup_deflist_.end())
            {
                double inv_T = tools::stod(_setup_deflist_.at("inv_temper"), "inverse temperature");
                if ((inv_T < 0.) || (inv_T > 1.))
                    throw std::out_of_range(
                        "Inverse temperature " + std::to_string(inv_T) + " out of range [0, 1].");
                else
                    setup_inv_T = inv_T;
            }
        }

        //===== Define Number of Threads (Monte_Carlo_Build) =====//
        set_n_threads(n_threads);

        //===== Pseudo-Random Number Generator Seeding =====//
        set_prng_seed(setup_seed);

        setup_c_sample_size = (setup_n_steps - setup_n_burnin) / setup_thin;
        setup_c_n_sweeps = setup_n_adapt / setup_thin;
        setup_c_n_sweeps += setup_n_burnin / setup_thin;
        setup_c_n_sweeps += setup_c_sample_size;
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Read parameter definitions and run assignment. */
void Monte_Carlo_Build::read_params(const std::string &pathtofile)
{
    tools::read_csv_table(pathtofile, _param_deftable_, param_file_headers);
    build_params();
    _report_pathtofile_["params"] = pathtofile;
}

/* Assign parameter definitions to dictionary. */
void Monte_Carlo_Build::build_params()
{
    try
    {
        if (_param_deftable_.empty())
        {
            // error if nothing was previously defined
            if (param_to_be_static.empty() && param_to_be_inferred.empty())
                throw std::invalid_argument(
                    "Undefined parameter definitions.");
        }
        else
        {
            //===== Read Parameter Definitions =====//
            size_t table_length = _param_deftable_.at("infer").size();

            //===== Empty current list of parameters =====//
            param_map.clear();
            param_to_be_static.clear();
            param_to_be_inferred.clear();

            //===== Load parameters from the row with options =====//
            for (size_t k = 0; k < table_length; ++k)
            {
                const std::string &name = _param_deftable_.at("name")[k]; // parameter name
                double value = INF_DOUBLE;                                // init value

                // FIXME : temporary
                // parameter identifier
                parID p = name;

                // add entry to shared parameter dictionary
                if (tools::stobool(_param_deftable_.at("infer")[k]))
                {
                    // load value
                    if (!tools::is_na(_param_deftable_.at("value")[k]))
                        value = tools::stod(_param_deftable_.at("value")[k], "parameter");

                    // load boolean flag for discrete variables
                    bool is_discrete = false;
                    if (!tools::is_na(_param_deftable_.at("discrete")[k]))
                        is_discrete = tools::stobool(_param_deftable_.at("discrete")[k]);

                    // load prior name
                    const std::string &prior_name = _param_deftable_.at("prior")[k];
                    const dict_args &prior_args = tools::merge_list_args(_param_deftable_.at("args")[k],
                                                                         _param_deftable_.at("argsval")[k]);

                    // load range
                    dict_args range = {};
                    if (!tools::is_na(_param_deftable_.at("min")[k]))
                        range.emplace("min", tools::stod(_param_deftable_.at("min")[k], "parameter minimum"));

                    if (!tools::is_na(_param_deftable_.at("max")[k]))
                        range.emplace("max", tools::stod(_param_deftable_.at("max")[k], "parameter maximum"));

                    // load rule name
                    const std::string &rule_name = _param_deftable_.at("rule")[k];

                    // load rule scale
                    double rule_scale = INF_DOUBLE;
                    if (!tools::is_na(_param_deftable_.at("rule_scale")[k]))
                        rule_scale = tools::stod(_param_deftable_.at("rule_scale")[k]);

                    // assign parameter to the set of variables
                    Parameter new_param(name, value, is_discrete, prior_name,
                                        prior_args, range, rule_name, rule_scale);
                    assign_param(p, std::move(new_param));
                }
                else
                {
                    // assign constant parameter to the set of variables
                    if (tools::is_na(_param_deftable_.at("value")[k]))
                        throw std::invalid_argument(
                            "Undefined value for constant parameter '" + name + "'.");
                    else
                        value = tools::stod(_param_deftable_.at("value")[k], "parameter");
                    Parameter new_param(name, value);
                    assign_param(p, std::move(new_param));
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }
}

/* Extract average and variance from vector of global loglike values. */
void Monte_Carlo_Build::dgn_glob_ll_stats(double &avg, double &var) const { tools::compute_stats(dgn_glob_loglike, avg, var); }

double Monte_Carlo_Build::get_global_loglike() const { return global_loglike; };
double Monte_Carlo_Build::get_delta_global_loglike() const { return global_loglike - buffer; }
void Monte_Carlo_Build::restore_global_loglike() { global_loglike = buffer; }
void Monte_Carlo_Build::increment_global_loglike(double delta) { global_loglike += delta; };

void Monte_Carlo_Build::set_global_loglike(double new_ll)
{
    buffer = global_loglike;
    global_loglike = new_ll;
}

/* Return the pseudo random number generator. */
rnd_gen &Monte_Carlo_Build::get_prng() { return (this->*_func_prng_)(); }

/* Draw a random variable uniformly distributed. */
double Monte_Carlo_Build::uniform() { return _uniform_distr_(get_prng()); }

/* String with parameter definitions at the current state. */
std::string Monte_Carlo_Build::show_param(const parID &p) const { return p + C_DELIM + param_map.at(p).show(); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : show setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the complete choices for the setup. */
void Monte_Carlo_Build::show_setup() const
{
    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nn_threads" << C_DELIM << setup_n_threads
        << "\nn_chains" << C_DELIM << setup_n_chains
        << "\nn_adapt" << C_DELIM << setup_n_adapt
        << "\nn_burnin" << C_DELIM << setup_n_burnin
        << "\nn_steps" << C_DELIM << setup_n_steps
        << "\nthin" << C_DELIM << setup_thin
        << "\nrnd_init_params" << C_DELIM << tools::truth(setup_rnd_init_params)
        << "\ninv_temper" << C_DELIM << setup_inv_T
        << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void Monte_Carlo_Build::show_all_params() const
{
    std::ostringstream oss;
    oss << param_file_headers[0];

    for (size_t k = 1; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];

    oss << "\n";

    for (const parID &p : param_to_be_inferred)
        oss << show_param(p) << "\n";

    for (const parID &p : param_to_be_static)
        oss << show_param(p) << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void Monte_Carlo_Build::show_results() const
{
    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";

    for (size_t k = 2; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];

    oss << "\n";

    for (const parID &p : param_to_be_inferred)
        oss << p << C_DELIM << param_map.at(p).show_diagnostics() << "\n";

    tools::print_stdout(oss.str());
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void Monte_Carlo_Build::print_infer_state(const std::string &output_filename) const
{
    if (!tools::is_na(output_filename))
    {
        std::ofstream inferstate_ofstream;
        inferstate_ofstream.open(output_filename.c_str());

        if (inferstate_ofstream.is_open())
        {
            //----- HEADERS -----//
            // header loglikelihood
            inferstate_ofstream << "loglikelihood";

            // header parameters
            for (const parID &p : param_to_be_inferred)
                inferstate_ofstream << C_DELIM << p;

            // header chain
            if (setup_n_chains > 1)
                inferstate_ofstream << C_DELIM << "CHAIN__";

            inferstate_ofstream << std::endl;

            //----- VALUES -----//

            // FIXME : // HERE : to be checked
            for (size_t c = 0; c < setup_n_chains; ++c)
            {
                for (size_t k = 0; k < setup_c_sample_size; k++)
                {
                    size_t i = k + c * setup_c_sample_size;
                    // print loglikelihood
                    inferstate_ofstream << std::setprecision(2) << std::fixed << dgn_glob_loglike[i];

                    // print parameter values
                    for (const parID &p : param_to_be_inferred)
                        inferstate_ofstream << C_DELIM << std::setprecision(5)
                                            << std::fixed << param_map.at(p).dgn_value[i];

                    // print which chain
                    if (setup_n_chains > 1)
                        inferstate_ofstream << C_DELIM << (c + 1);

                    inferstate_ofstream << std::endl;
                }
            }

            inferstate_ofstream.close();
        }
        else
        {
            throw std::runtime_error(
                "Unable to open file: " + output_filename);
        }
    }
}
