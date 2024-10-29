/*
 * core.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "core.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Welford online algorithm for mean and variance. */
void compute_stats(const std::vector<double> &data, double &avg, double &var)
{
    size_t N = data.size();
    double x, delta, delta2;
    double sqr_diff = 0.0;

    avg = 0.0;
    for (size_t i = 0; i < N; ++i)
    {
        x = data[i];
        delta = x - avg;
        avg += delta / (i + 1);
        delta2 = x - avg;
        sqr_diff += delta * delta2;
    }
    var = sqr_diff / (N - 1);
}

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
void Parameter::update_diagnostics()
{
    dgn_value[dgn_count++] = get_value();
}

/* Message to be used at the moment of display diagnostics. (TEMPORARY) */
std::string Parameter::show_diagnostics() const
{
    using namespace tools;
    std::ostringstream oss;

    double avg = 0, var = 0;
    compute_stats(dgn_value, avg, var);

    oss << avg << C_DELIM << std::sqrt(var) << C_DELIM << truth(infer)
        << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
        << C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
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
//////////////////////////////  Monte_Carlo_Core : methods (getter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Safe method to get a value with previous assertion. */
double Monte_Carlo_Core::_get_param_assert_(const std::string &parID) const
{
    assert_param(parID);
    return _get_param_default_(parID);
}

/* Safe method to get a pointer to a parameter with previous assertion. */
const double *Monte_Carlo_Core::_get_param_pointer_assert_(const std::string &parID) const
{
    assert_param(parID);
    return _get_param_pointer_default_(parID);
}

/* Disable parameter assertion. */
void Monte_Carlo_Core::disable_assert_param()
{
    _get_param_method_ = &Monte_Carlo_Core::_get_param_default_;
    _get_param_pointer_method_ = &Monte_Carlo_Core::_get_param_pointer_default_;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a new parameter to the dictionary. */
void Monte_Carlo_Core::assign_param(const std::string &parID, const Parameter &new_param)
{
    // create the new parameter if it's not already defined
    if (params.find(parID) != params.end())
    {
        throw std::invalid_argument("Parameter '" + parID + "' has multiple definitions.");
    }
    else
    {
        // emplace the new parameter
        params.emplace(parID, new_param);

        // update the lists of parameters to be inferred, and these to be static
        if (new_param.get_infer())
        {
            param_to_be_inferred.push_back(parID);
        }
        else
        {
            param_to_be_static.push_back(parID);
        }
    }
}

/* Add a new scalar for augmented typedict. TODO : is there a more elegant option? */
const double *Monte_Carlo_Core::assign_scalar_pointer(double val)
{
    if (_scalar_type_params_.find(val) == _scalar_type_params_.end())
    {
        Parameter new_param(UNDEFINED, val);
        _scalar_type_params_.emplace(val, new_param);
    }
    /* Add a new scalar for augmented typedict FIXME : temporary. */
    return _scalar_type_params_.at(val).get_pointer();
}

/* Set the number of threads for parrell computations (default is 1). */
void Monte_Carlo_Core::set_num_threads(unsigned n)
{
    if (n == 0)
    {
        throw std::invalid_argument("Number of threads must be equal or greater than 1.");
    }

#ifdef _OPENMP
    if (n > _max_num_available_threads)
    {
        std::ostringstream oss;
        oss << "Number of threads (" << n << ") changed to number of available threads (" << _max_num_available_threads << ").";
        tools::warning(oss.str());
        n = _max_num_available_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n);
    setup_num_threads = n;

    // create different generators for different parallel threads
    // UPDATE: they are used only to parallelize data augmentantion
    if (n > 1)
    {
        prng_parallel.reserve(n);
        _prng_method_ = &Monte_Carlo_Core::_prng_parallel_threads_;
    }
    else
    {
        _prng_method_ = &Monte_Carlo_Core::_prng_unique_thread_;
    }
#else
    if (n != 1)
    {
        tools::warning("Parallelization with OpenMP not supported.");
    }
    setup_num_threads = 1;
    _prng_method_ = &Monte_Carlo_Core::_prng_unique_thread_;
#endif
}

/* Assign seed for random generators. */
void Monte_Carlo_Core::set_prng_seed(int seed)
{
    // seed for the unique thread rand num generator
    prng.seed(seed);

    // seeds for the case of parallelization
    if (setup_num_threads > 1)
    {
        std::uniform_int_distribution<unsigned long int> seed_generator(0, 0XFFFFFFFF);
        for (unsigned n = 0; n < setup_num_threads; ++n)
        {
            prng_parallel[n].seed(seed_generator(prng));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Core : methods (setter)  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
void Monte_Carlo_Core::assert_param(const std::string &parID) const
{
    if (params.find(parID) == params.end())
    {
        throw std::invalid_argument("Parameter '" + parID + "' not found. Check your model and parameter definitions.");
    }
}

/* randomly initialize values of the parameters to be inferred */
void Monte_Carlo_Core::init_rand_params()
{
    if (setup_rnd_init_params)
    {
        for (const std::string &parID : param_to_be_inferred)
        {
            params.at(parID).rand(prng);
        }
    }
    else
    {
        if (_is_undef_params_)
        {
            tools::error("Parameters with undefined values are not supported.");
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Monte_Carlo_Core : methods (parameter moves)  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Monte_Carlo_Core::init_param_move_counts()
{
    for (const std::string &parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }
}

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
// FIXME : Gibbs sampling is missing here!
// 1) always valid 2) move_logratio == 0
bool Monte_Carlo_Core::propose_valid_param_move(const std::string &parID)
{
    // set new value of the parameter and store the old to buffer
    params.at(parID).propose(prng);

    // proceed only if valid move, otherwise restore buffer
    if (params.at(parID).valid_move())
    {
        return true;
    }
    else
    {
        // reassign buffered value to param
        restore_param_value(parID);
        return false;
    }
}

/* Init diagnostics for all parameters */
void Monte_Carlo_Core::init_param_diagnostics(ulong sample_size)
{
    for (const std::string &parID : param_to_be_inferred)
    {
        params.at(parID).init_diagnostics(sample_size);
    }
}

/* Update diagnostics for all parameters */
void Monte_Carlo_Core::update_param_diagnostics()
{
    for (const std::string &parID : param_to_be_inferred)
    {
        params.at(parID).update_diagnostics();
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef USE_RCPP
/* Import setup from R list and run build. */
void Monte_Carlo_Core::import_setup(Rcpp::List setup)
{
    tools::rlist_to_map(setup, _setup_deflist_);
    build_setup();
}

/* Import params from R DataFrame and run build. */
void Monte_Carlo_Core::import_params(Rcpp::DataFrame params)
{
    tools::rdataframe_to_map(params, _param_deftable_);
    build_params();
}

#endif

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Read Setup file and produce a definition list. */
void Monte_Carlo_Core::read_setup(const std::string &pathtofile)
{
    tools::read_csv_list(pathtofile, _setup_deflist_);
    build_setup();
}

/* Assign setup. */
void Monte_Carlo_Core::build_setup()
{
    try
    {
        // default number of threads
        unsigned num_threads = 1;

        // default prng seed
        //----- default prng seed -----//
        auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();
        setup_seed = static_cast<unsigned long long>(time_now);

        if (!_setup_deflist_.empty())
        {
            //===== Check that all provided setups are correct =====//
            std::vector<std::string> file_entries = {"num_threads", "seed", "num_chains",
                                                     "num_iter", "num_adapt", "thinning",
                                                     "burn_in", "rnd_init_params",
                                                     "num_rungs", "inv_temper"};
            std::vector<std::string> wrong_entries = {};
            for (const auto &pair : _setup_deflist_)
            {
                if (find(file_entries.begin(), file_entries.end(), pair.first) == file_entries.end())
                {
                    wrong_entries.push_back(pair.first);
                }
            }
            if (wrong_entries.size() > 0)
            {
                std::ostringstream oss;
                oss << "Unknown setup option(s): '" << wrong_entries[0] << "'";
                for (size_t k = 1; k < wrong_entries.size(); ++k)
                {
                    oss << ", '" << wrong_entries[k] << "'";
                }
                oss << ". Available setup options are: '" << file_entries[0] << "'";
                for (size_t k = 1; k < file_entries.size(); ++k)
                {
                    oss << ", '" << file_entries[k] << "'";
                }
                tools::warning(oss.str() + ".");
            }

            //===== ASSIGN SETUP CONFIGURATIONS =====//

            //===== NUMBER OF THREADS =====//
            if (_setup_deflist_.find("num_threads") != _setup_deflist_.end())
            {
                num_threads = tools::stoi(_setup_deflist_.at("num_threads"), "number of threads");
            }

            //===== PRNG SEED =====//
            if ((_setup_deflist_.find("seed") != _setup_deflist_.end()) && (!tools::is_na(_setup_deflist_.at("seed"))))
            {
                setup_seed = tools::stoull(_setup_deflist_.at("seed"), "random seed");
            }

            //===== NUMBER OF CHAINS =====//
            if (_setup_deflist_.find("num_chains") != _setup_deflist_.end())
            {
                setup_num_chains = tools::stoul(_setup_deflist_.at("num_chains"), "number of chains");
                if (setup_num_iter < 1)
                {
                    throw std::out_of_range("Number of chains must be greater than 0.");
                }
            }

            //===== NUMBER OF ITERATIONS =====//
            if (_setup_deflist_.find("num_iter") != _setup_deflist_.end())
            {
                setup_num_iter = tools::stoul(_setup_deflist_.at("num_iter"), "number of iterations");
                if (setup_num_iter < 2)
                {
                    throw std::out_of_range("Number of iterations must be greater than 1.");
                }
            }

            //===== NUMBER OF (FINITE) ADAPTIVE ITERATIONS =====//
            if (_setup_deflist_.find("num_adapt") != _setup_deflist_.end())
            {
                setup_num_adapt = tools::stoul(_setup_deflist_.at("num_adapt"), "number of adaptive steps");
            }

            //===== THINNING =====//
            if (_setup_deflist_.find("thinning") != _setup_deflist_.end())
            {
                setup_thinning = tools::stoul(_setup_deflist_.at("thinning"), "thinning");
                if (setup_thinning >= setup_num_iter)
                {
                    throw std::out_of_range("Thinning must be smaller than number of iterations. ");
                }
                else if (setup_thinning == 0)
                {
                    throw std::out_of_range("Thinning must be greater than 0.");
                }
            }

            //===== BURN-IN =====//
            if (_setup_deflist_.find("burn_in") != _setup_deflist_.end())
            {
                setup_burn_in = tools::stoul(_setup_deflist_.at("burn_in"), "burn-in");
                if (setup_burn_in >= setup_num_iter)
                {
                    throw std::out_of_range("Burn-in must be smaller than number of iterations.");
                }
            }

            //===== RANDOM INITIALIZATION OF PARAMETERS =====//
            if (_setup_deflist_.find("rnd_init_params") != _setup_deflist_.end())
            {
                setup_rnd_init_params = tools::stobool(_setup_deflist_.at("rnd_init_params"));
                if ((setup_num_chains > 1) && (setup_rnd_init_params == false))
                {
                    tools::warning("Unable to fix initial parameters when running multiple chains.");
                    setup_rnd_init_params = true;
                }
            }

            //===== INVERSE TEMPERATURE =====//
            if (_setup_deflist_.find("inv_temper") != _setup_deflist_.end())
            {
                double inv_T = tools::stod(_setup_deflist_.at("inv_temper"), "inverse temperature");
                if ((inv_T < 0.) || (inv_T > 1.))
                {
                    std::ostringstream oss;
                    oss << "Inverse temperature " << inv_T << " out of range [0, 1].";
                    throw std::out_of_range(oss.str());
                }
                else
                {
                    setup_inv_T = inv_T;
                }
            }
        }

        //===== Define Number of Threads (Monte_Carlo_Core) =====//
        set_num_threads(num_threads);

        //===== Pseudo-Random Number Generator Seeding =====//
        set_prng_seed(setup_seed);

        setup_c_sample_size = (setup_num_iter - setup_burn_in) / setup_thinning;
        setup_c_num_sweeps = setup_num_adapt / setup_thinning;
        setup_c_num_sweeps += setup_burn_in / setup_thinning;
        setup_c_num_sweeps += setup_c_sample_size;
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
void Monte_Carlo_Core::read_params(const std::string &pathtofile)
{
    tools::read_csv_table(pathtofile, _param_deftable_, param_file_headers);
    build_params();
}

/* Assign parameter definitions to dictionary. */
void Monte_Carlo_Core::build_params()
{
    try
    {
        if (_param_deftable_.empty())
        {
            // error if nothing was previously defined
            if (param_to_be_static.empty() && param_to_be_inferred.empty())
            {
                throw std::invalid_argument("Undefined parameter definitions.");
            }
        }
        else
        {
            //===== Read Parameter Definitions =====//
            size_t table_length = _param_deftable_.at("infer").size();

            //===== Empty current list of parameters =====//
            params.clear();
            param_to_be_static.clear();
            param_to_be_inferred.clear();

            //===== Load parameters from the row with options =====//
            for (size_t p = 0; p < table_length; ++p)
            {
                std::string parID = _param_deftable_.at("name")[p];
                double value = INF_DOUBLE;

                // add entry to shared parameter dictionary
                if (tools::stobool(_param_deftable_.at("infer")[p]))
                {
                    // load value
                    if (tools::is_na(_param_deftable_.at("value")[p]))
                    {
                        if (setup_rnd_init_params)
                        {
                            _is_undef_params_ = true;
                        }
                        else
                        {
                            std::ostringstream oss;
                            oss << "Undefined value without random initialization for parameter to infer '" << parID << "'.";
                            throw std::invalid_argument(oss.str());
                        }
                    }
                    else
                    {
                        value = tools::stod(_param_deftable_.at("value")[p], "parameter");
                    }

                    // load range
                    dictargs range = {};
                    if (!tools::is_na(_param_deftable_.at("min")[p]))
                    {
                        range.emplace("min", tools::stod(_param_deftable_.at("min")[p], "parameter minimum"));
                    }
                    if (!tools::is_na(_param_deftable_.at("max")[p]))
                    {
                        range.emplace("max", tools::stod(_param_deftable_.at("max")[p], "parameter maximum"));
                    }

                    // assign parameter to the set of variables
                    Parameter new_param(
                        parID, value,
                        _param_deftable_.at("prior")[p],
                        tools::merge_list_args(_param_deftable_.at("args")[p],
                                               _param_deftable_.at("argsval")[p]),
                        range,
                        _param_deftable_.at("rule")[p],
                        tools::merge_list_args(_param_deftable_.at("rule_args")[p],
                                               _param_deftable_.at("rule_argsval")[p]));
                    assign_param(parID, new_param);
                }
                else
                {
                    // load value
                    if (tools::is_na(_param_deftable_.at("value")[p]))
                    {
                        throw std::invalid_argument("Undefined value for immutable parameter '" + parID + "'.");
                    }
                    else
                    {
                        value = tools::stod(_param_deftable_.at("value")[p], "parameter");
                    }

                    // assign immutable parameter to the set of variables
                    Parameter new_param(parID, value);
                    assign_param(parID, new_param);
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
void Monte_Carlo_Core::dgn_glob_ll_stats(double &avg, double &var) const
{
    compute_stats(dgn_glob_loglike, avg, var);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef USE_RCPP

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::String stoRString(const std::string &x)
{
    return (tools::is_na(x)) ? NA_STRING : Rcpp::String(x);
}

/* Export current setup configurations to R List. */
Rcpp::List Monte_Carlo_Core::export_setup() const
{
    using namespace Rcpp;

    return List::create(
        _["seed"] = setup_seed,
        _["num_threads"] = setup_num_threads,
        _["num_chains"] = setup_num_chains,
        _["num_adapt"] = setup_num_adapt,
        _["burn_in"] = setup_burn_in,
        _["num_iter"] = setup_num_iter,
        _["thinning"] = setup_thinning,
        _["rnd_init_params"] = setup_rnd_init_params,
        _["inv_temper"] = setup_inv_T);
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame Monte_Carlo_Core::export_params() const
{
    using namespace Rcpp;

    std::vector<std::string> V_names = param_to_be_inferred;
    V_names.insert(V_names.end(), param_to_be_static.begin(), param_to_be_static.end());

    NumericVector V_value = {};
    LogicalVector V_infer = {};
    CharacterVector V_prior = {}, V_args = {}, V_argsval = {};
    CharacterVector V_rule = {}, V_rule_args = {}, V_rule_argsval = {};
    NumericVector V_min = {}, V_max = {};

    for (const std::string &parID : V_names)
    {
        V_value.push_back(get_param(parID));
        V_infer.push_back(get_infer(parID));

        V_prior.push_back(stoRString(params.at(parID).show_distr_name()));
        V_args.push_back(stoRString(params.at(parID).show_distr_arg_names()));
        V_argsval.push_back(stoRString(params.at(parID).show_distr_arg_values()));

        V_rule.push_back(stoRString(params.at(parID).show_prop_name()));
        V_rule_args.push_back(stoRString(params.at(parID).show_prop_arg_names()));
        V_rule_argsval.push_back(stoRString(params.at(parID).show_prop_arg_values()));

        const std::string min_val = params.at(parID).show_range_min();
        const std::string max_val = params.at(parID).show_range_max();
        V_min.push_back(tools::is_na(min_val) ? NA_REAL : tools::stod(min_val));
        V_max.push_back(tools::is_na(max_val) ? NA_REAL : tools::stod(max_val));
    }

    return DataFrame::create(_["name"] = V_names, _["value"] = V_value, _["infer"] = V_infer,
                             _["prior"] = V_prior, _["args"] = V_args, _["argsval"] = V_argsval,
                             _["rule"] = V_rule, _["rule_args"] = V_rule_args, _["rule_argsval"] = V_rule_argsval,
                             _["min"] = V_min, _["max"] = V_max);
}

/* Export MCMC traces with posterior samples. */
Rcpp::DataFrame Monte_Carlo_Core::export_traces() const
{
    using namespace Rcpp;

    DataFrame output = DataFrame::create(
        _["loglikelihood"] = wrap(dgn_glob_loglike));

    for (const std::string &parID : param_to_be_inferred)
    {
        output.push_back(wrap(params.at(parID).get_dgn_value()), parID);
    }

    // retrieve chain identities
    if (setup_num_chains > 1)
    {
        IntegerVector df_chain(setup_c_sample_size * setup_num_chains);
        for (ulong c = 0; c < setup_num_chains; ++c)
        {
            std::fill(df_chain.begin() + c * setup_c_sample_size,
                      df_chain.begin() + (c + 1) * setup_c_sample_size, c + 1);
        }
        output.push_back(df_chain, "CHAIN__");
    }

    return output;
}

#else

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
void Monte_Carlo_Core::show_setup() const
{
    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nnum_threads" << C_DELIM << setup_num_threads
        << "\nnum_chains" << C_DELIM << setup_num_chains
        << "\nnum_adapt" << C_DELIM << setup_num_adapt
        << "\nburn_in" << C_DELIM << setup_burn_in
        << "\nnum_iter" << C_DELIM << setup_num_iter
        << "\nthinning" << C_DELIM << setup_thinning
        << "\nrnd_init_params" << C_DELIM << tools::truth(setup_rnd_init_params)
        << "\ninv_temper" << C_DELIM << setup_inv_T
        << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void Monte_Carlo_Core::show_all_params() const
{
    std::ostringstream oss;
    oss << param_file_headers[0];
    for (size_t k = 1; k < param_file_headers.size(); ++k)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    oss << "\n";

    for (const std::string &parID : param_to_be_inferred)
    {
        oss << show_param(parID) << "\n";
    }
    for (const std::string &parID : param_to_be_static)
    {
        oss << show_param(parID) << "\n";
    }

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void Monte_Carlo_Core::show_results() const
{
    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";
    for (size_t k = 2; k < param_file_headers.size(); ++k)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    oss << "\n";
    for (const std::string &parID : param_to_be_inferred)
    {
        oss << parID << C_DELIM << params.at(parID).show_diagnostics() << "\n";
    }

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void Monte_Carlo_Core::print_infer_state(const std::string &output_filename) const
{
    if (!tools::is_na(output_filename))
    {
        std::ofstream inferstate_ofstream;
        inferstate_ofstream.open(output_filename.c_str());

        if (inferstate_ofstream.is_open())
        {
            //----- HEADERS -----//
            // header : loglikelihood
            inferstate_ofstream << "loglikelihood";
            // header : parameters
            for (const std::string &parID : param_to_be_inferred)
            {
                inferstate_ofstream << C_DELIM << parID;
            }
            // header : chain
            if (setup_num_chains > 1)
            {
                inferstate_ofstream << C_DELIM << "CHAIN__";
            }
            inferstate_ofstream << std::endl;

            //----- VALUES -----//

            // FIXME : // FLAG : // HERE :
            for (size_t c = 0; c < setup_num_chains; ++c)
            {
                for (size_t k = 0; k < setup_c_sample_size; k++)
                {
                    size_t i = k + c * setup_c_sample_size;
                    // print : loglikelihood
                    inferstate_ofstream << std::setprecision(2) << std::fixed << dgn_glob_loglike[i];
                    // print : parameter values
                    for (const std::string &parID : param_to_be_inferred)
                    {
                        // FIXME
                        inferstate_ofstream << C_DELIM << std::setprecision(5) << std::fixed << params.at(parID).get_dgn_value(i);
                    }
                    // print : which chain
                    if (setup_num_chains > 1)
                    {
                        inferstate_ofstream << C_DELIM << (c + 1);
                    }
                    inferstate_ofstream << std::endl;
                }
            }

            inferstate_ofstream.close();
        }
        else
        {
            tools::error("Unable to open file: " + output_filename);
        }
    }
}

#endif