/*
 * mcmc_core.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

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

/* Set the number of threads for parrell computations (default is 1). */
void MCMC::Impl::set_n_threads(unsigned n)
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
    // note: they are used only to parallelize data augmentantion
    if (n > 1)
    {
        prng_parallel.reserve(n);
        _func_prng_ = &Impl::_prng_parallel_threads_;
    }
    else
    {
        _func_prng_ = &Impl::_prng_unique_thread_;
    }
#else
    if (n != 1)
        tools::warning("Parallelization with OpenMP not supported.");

    setup_n_threads = 1;
    _func_prng_ = &Impl::_prng_unique_thread_;
#endif
}

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

/* Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::init_parameters()
{
    for (const parID &p : param_to_be_inferred)
    {
        // restore user initialization (meaningful for multiple chains)
        param_map.at(p).user_reset();

        if (setup_rand_init_params)
        {
            if (tools::is_na(param_map.at(p).get_distr_name()))
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                param_map.at(p).rand(prng);
        }
        else
        {
            if (tools::is_na(param_map.at(p).get_value()))
                throw std::invalid_argument(
                    "Undefined value for parameter to infer \"" + p + "\" while no random initialization.");
        }
    }
}

/* Set accept rates to zero. */
void MCMC::Impl::init_move_counts()
{
    for (const parID &p : param_to_be_inferred)
        param_map.at(p).reset_accept_rate();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  MCMC::Impl : methods (parameter moves)  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
bool MCMC::Impl::propose_valid_param_move(const parID &p)
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
void MCMC::Impl::restore_param_value(const parID &p) { param_map.at(p).restore_value(); }

/* Update the acceptance rate of the paramter. */
void MCMC::Impl::update_param_accept_rate(const parID &p, bool is_accepted) { param_map.at(p).update_accept_rate(is_accepted); }

/* Metropolis logratio associate to the current move. */
double MCMC::Impl::param_move_logratio(const parID &p) const { return param_map.at(p).get_move_logratio(); }

/* Adapt the parameter Metropolis move. */
void MCMC::Impl::adapt_param_rule(const parID &p, double ratio) { param_map.at(p).adapt_rule(ratio); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Set inverse temperature and its relative functional methods. */
void MCMC::Impl::set_inver_temper(double inv_T)
{
    if ((inv_T < 0.) || (inv_T > 1.))
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " out of range [0, 1].");
    else
        setup_inv_T = inv_T;
}

/* Read Setup file and produce a definition list. */
void MCMC::Impl::read_setup(const std::string &pathtofile)
{
    tools::read_csv_list(pathtofile, _setup_deflist_);
    build_setup();
    _report_pathtofile_["setup"] = pathtofile;
}

/* Assign setup. */
void MCMC::Impl::build_setup()
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
                                    "n_burnin", "rand_init",
                                    "n_rungs", "inv_temper", "id_trace"};

            str_vec wrong_entries = {};

            for (const auto &pair : _setup_deflist_)
                if (find(file_entries.begin(), file_entries.end(), pair.first) == file_entries.end())
                    wrong_entries.push_back(pair.first);

            if (wrong_entries.size() > 0)
            {
                std::ostringstream oss;
                oss << "Unknown setup option(s): \"" << wrong_entries[0] << "\"";
                for (size_t k = 1; k < wrong_entries.size(); ++k)
                    oss << ", \"" << wrong_entries[k] << "\"";
                oss << ". Available setup options are: \"" << file_entries[0] << "\"";
                for (size_t k = 1; k < file_entries.size(); ++k)
                    oss << ", \"" << file_entries[k] << "\"";
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
                if (setup_n_steps < 1)
                    throw std::out_of_range(
                        "Number of iterations must be greater than 0.");
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
            if (_setup_deflist_.find("rand_init") != _setup_deflist_.end())
            {
                setup_rand_init_params = tools::stobool(_setup_deflist_.at("rand_init"));
                if ((setup_n_chains > 1) && (setup_rand_init_params == false))
                {
                    tools::warning(
                        "Unable to fix initial parameters when running multiple chains.");
                    setup_rand_init_params = true;
                }
            }

            //===== INVERSE TEMPERATURE =====//
            if (_setup_deflist_.find("inv_temper") != _setup_deflist_.end())
            {
                double tmp = tools::stod(_setup_deflist_.at("inv_temper"), "inverse temperature");
                set_inver_temper(tmp);
            }

            //===== INDIVIDUAL ID AUGMENTED TRACE =====//
            if (_setup_deflist_.find("id_trace") != _setup_deflist_.end())
            {
                setup_id_trace = _setup_deflist_.at("id_trace");
            }
        }

        //===== Define Number of Threads (MCMC::Impl) =====//
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
void MCMC::Impl::read_params(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_table(pathtofile, _param_deftable_, param_file_headers);
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While reading parameter definitions: " + std::string(e.what()));
    }

    // use loaded definitions to build
    build_params();
    // store the file name in the report
    _report_pathtofile_["params"] = pathtofile;
}

/* Assign parameter definitions to dictionary. */
void MCMC::Impl::build_params()
{
    try
    {
        if (_param_deftable_.empty())
        {
            // error if nothing was previously defined
            if (param_to_be_constant.empty() && param_to_be_inferred.empty())
                throw std::invalid_argument(
                    "Undefined parameter definitions.");
        }
        else
        {
            //===== Read Parameter Definitions =====//
            size_t table_length = _param_deftable_.at("name").size();

            //===== Empty current list of parameters =====//
            param_map.clear();
            param_to_be_constant.clear();
            param_to_be_inferred.clear();

            //===== Load parameters from the row with options =====//
            for (size_t k = 0; k < table_length; ++k)
            {
                const std::string &name = _param_deftable_.at("name")[k]; // parameter name
                double value = NULL_DOUBLE;                          // init value

                // UPDATE : temporary parameter identifier
                // NOTE : id parID changes this should change too
                parID p = name;                 

                //----- Add entry to shared parameter dictionary -----//
                if (param_map.find(p) != param_map.end())
                {
                    // create the new parameter if it's not already defined
                    throw std::invalid_argument(
                        "Parameter \"" + p + "\" has multiple definitions.");
                }
                else if (tools::stobool(_param_deftable_.at("constant")[k]))
                {
                    // assign constant parameter to the set of variables
                    if (tools::is_na(_param_deftable_.at("value")[k]))
                        throw std::invalid_argument(
                            "Undefined value for constant parameter \"" + name + "\".");
                    else
                        value = tools::stod(_param_deftable_.at("value")[k], "parameter");

                    // assign constant parameter to the set of variables
                    param_map.emplace(p, Base_Variable(name, value));
                    // update the lists of constant parameters
                    param_to_be_constant.push_back(p);
                }
                else
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

                    // load range
                    double lower = NULL_DOUBLE;
                    if (!tools::is_na(_param_deftable_.at("min")[k]))
                        lower = tools::stod(_param_deftable_.at("min")[k], "parameter minimum");

                    double upper = NULL_DOUBLE;
                    if (!tools::is_na(_param_deftable_.at("max")[k]))
                        upper = tools::stod(_param_deftable_.at("max")[k], "parameter maximum");

                    // load rule name
                    std::string rule_name = "gibbs";
                    if (!tools::is_na(_param_deftable_.at("rule")[k]))
                        rule_name = _param_deftable_.at("rule")[k];

                    // load rule scale
                    double rule_scale = 1.;
                    if (!tools::is_na(_param_deftable_.at("rule_scale")[k]))
                        rule_scale = tools::stod(_param_deftable_.at("rule_scale")[k]);

                    //----- Create a new Parameter Variable -----//
                    param_map.emplace(p, Base_Variable(name, value, is_discrete));
                    Base_Variable &new_param = param_map.at(p);

                    // assign a prior distribution to the parameter
                    new_param.assign_distribution(prior_name);

                    // assign scalar arguments to prior distribution
                    str_vec arg_name = tools::parse_list_args(_param_deftable_.at("arg_name")[k]);
                    str_vec arg_value = tools::parse_list_args(_param_deftable_.at("arg_value")[k]);
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
                    param_to_be_inferred.push_back(p);
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

/* Extract average and variance from vector of global loglike values. */
void MCMC::Impl::dgn_glob_ll_stats(double &avg, double &var) const { tools::compute_stats(dgn_glob_loglike, avg, var); }

double MCMC::Impl::get_global_loglike() const { return global_loglike; };
double MCMC::Impl::get_delta_global_loglike() const { return global_loglike - buffer; }
void MCMC::Impl::restore_global_loglike() { global_loglike = buffer; }
void MCMC::Impl::increment_global_loglike(double delta) { global_loglike += delta; };

void MCMC::Impl::set_global_loglike(double new_ll)
{
    buffer = global_loglike;
    global_loglike = new_ll;
}

/* Return the pseudo random number generator. */
rnd_gen &MCMC::Impl::get_prng() { return (this->*_func_prng_)(); }

/* Draw a random variable uniformly distributed. */
double MCMC::Impl::uniform() { return _uniform_distr_(get_prng()); }

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
void MCMC::Impl::show_setup() const
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
        << "\nrand_init" << C_DELIM << tools::truth(setup_rand_init_params)
        << "\nid_trace" << C_DELIM << setup_id_trace
        << "\ninv_temper" << C_DELIM << setup_inv_T
        << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void MCMC::Impl::show_all_params() const
{
    std::ostringstream oss;
    oss << param_file_headers[0];

    for (size_t k = 1; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];

    oss << "\n";

    for (const parID &p : param_to_be_inferred)
        oss << param_map.at(p).show() << "\n";

    for (const parID &p : param_to_be_constant)
        oss << param_map.at(p).show() << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void MCMC::Impl::show_results() const
{
    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";

    for (size_t k = 2; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];

    oss << "\n";

    for (const parID &p : param_to_be_inferred)
    {
        const Base_Variable &this_param = param_map.at(p);
        double avg = 0, var = 0;
        tools::compute_stats(dgn_param_value.at(p), avg, var);

        oss << p
            << C_DELIM << avg
            << C_DELIM << std::sqrt(var)
            << C_DELIM << tools::truth(this_param.is_constant())
            << C_DELIM << tools::truth(this_param.is_discrete())
            << C_DELIM << this_param.get_distr_name()
            << C_DELIM << this_param.show_distr_arg_names()
            << C_DELIM << this_param.show_distr_arg_values()
            << C_DELIM << this_param.get_prop_name()
            << C_DELIM << this_param.show_prop_scale()
            << C_DELIM << tools::natod(this_param.show_range_min())
            << C_DELIM << tools::natod(this_param.show_range_max())
            << "\n";
    }

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
void MCMC::Impl::print_infer_state(const std::string &output_filename) const
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

            // augmented variables
            if (is_augm_infer() && (!tools::is_na(setup_id_trace)))
                for (const augmID &a : augm_list_to_be_inferred)
                    inferstate_ofstream << C_DELIM << a;

            // header chain
            if (setup_n_chains > 1)
                inferstate_ofstream << C_DELIM << "CHAIN__";

            inferstate_ofstream << std::endl;

            //----- VALUES -----//

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
                                            << std::fixed << dgn_param_value.at(p)[i];

                    // print selected augmented variables
                    if (is_augm_infer() && (!tools::is_na(setup_id_trace)))
                        for (const augmID &a : augm_list_to_be_inferred)
                            inferstate_ofstream << C_DELIM << std::setprecision(5)
                                                << std::fixed << dgn_augm_value.at(a)[i];

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
