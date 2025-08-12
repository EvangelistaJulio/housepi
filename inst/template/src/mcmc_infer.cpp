/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build_param.cpp
 *
 */

#include <atomic>

#include "mcmc_impl.h"

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
void MCMC::Impl::init_move_counts_()
{
    for (const std::string &p : param.list_sample)
        param.at(p).reset_accept_rate();
}

/* propose a valid move (retrun true), otherwise update refusal (return false) */
bool MCMC::Impl::propose_valid_param_move_(const std::string &p)
{
    Variable &this_par = param.at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(setup.prng());

    // recover buffer if not valid move
    if (!this_par.is_valid_move())
    {
        this_par.recover_value();
        return false;
    }

    return true;
}

/* Restore the parameter value in the buffer. */
void MCMC::Impl::recover_parameter_(const std::string &p, bool is_cache_depen)
{
    // recover buffer value of parameter
    param.at(p).recover_value();

    // recover buffer value of full loglike
    full_loglike_ = full_buffer_;

    // recover buffered loglikelihood in households
    for (Household &house : dataset.restricted())
    {
        house.recover_buffer_values(is_cache_depen);
        if (augm.is_any_sample())
            for (Individual &member : house.infected())
                member.recover_augm_model_loglike();
    }

    infec_profile.recover(is_cache_depen);
}

/* Metropolis-Hastings to accept or not current parameter move. */
bool MCMC::Impl::accept_valid_param_move_(const std::string &p)
{
    const double D_ll = full_loglike_ - full_buffer_;
    const double move_logratio = param.at(p).get_move_logratio();
    return setup.unif() < std::exp(setup.inv_T * D_ll + move_logratio);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : posterior sample  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize posterior sample of parameters. */
void MCMC::Impl::Parameter::init_draws(unsigned long sample_size)
{
    smpl_values.clear();
    smpl_values.reserve(list_sample.size());
    for (const std::string &p : list_sample)
    {
        smpl_values[p].clear();
        smpl_values[p].reserve(sample_size);
        smpl_values[p].assign(sample_size, HSP_NA_DOUBLE);
    }
}

/* Store current sample of parameters. */
void MCMC::Impl::Parameter::update_draws(unsigned long smpl_count)
{
    for (const std::string &p : list_sample)
        smpl_values[p][smpl_count] = get(p);
}

/* Clear posterior sample of parameters. */
void MCMC::Impl::Parameter::clear_draws()
{
    for (auto &pair : smpl_values)
    {
        pair.second.clear();
        pair.second.shrink_to_fit();
    }
    smpl_values.clear();
    smpl_values.rehash(0);
}

/* Initialize posterior sample of parameters. */
void MCMC::Impl::DataSet::init_draws(unsigned long sample_size)
{
    // sample over households
    smpl_values_.clear();
    smpl_values_.resize(restrict_indexes_.size());
    for (size_t R = 0; R < restrict_indexes_.size(); ++R)
    {
        smpl_values_[R].clear();
        smpl_values_[R].reserve(sample_size);
        smpl_values_[R].assign(sample_size, HSP_NA_DOUBLE);
    }
}

/* Store current sample of loglikelihood for restricted households. */
void MCMC::Impl::DataSet::update_draws(unsigned long smpl_count)
{
    // store sample for restricted households
    for (size_t R = 0; R < restrict_indexes_.size(); ++R)
    {
        Household &house = house_vector_[restrict_indexes_[R]];
        double x = house.get_obs_loglike() / house.multiplicity();
        for (const Individual &member : house.infected())
            x += member.get_augm_model_loglike();
        smpl_values_[R][smpl_count] = x;
    }
}

/* Clear posterior sample of dataset. */
void MCMC::Impl::DataSet::clear_draws()
{
    for (size_t R = 0; R < restrict_indexes_.size(); ++R)
    {
        smpl_values_[R].clear();
        smpl_values_[R].shrink_to_fit();
    }
    smpl_values_.clear();
    smpl_values_.shrink_to_fit();
}

/* Init posterior sampling. */
void MCMC::Impl::init_posterior_sample_()
{
    unsigned long sample_size = setup.c_sample_size * setup.n_chains;

    // diagnostics over full loglikelihood
    smpl_count_ = 0;
    smpl_full_loglike_.clear();
    smpl_full_loglike_.reserve(sample_size);
    smpl_full_loglike_.assign(sample_size, HSP_NA_DOUBLE);

    // initialize  other diagnostic vectors
    // UPDATE: I might want to store also power posterior sample
    if (!setup.is_power_post())
    {
        dataset.init_draws(sample_size);
        param.init_draws(sample_size);
        augm.init_draws(sample_size);
    }
}

/* Update Bayesian posterior sample. */
void MCMC::Impl::update_posterior_sample_()
{
    // store sample of the full loglikelihood
    smpl_full_loglike_[smpl_count_] = full_loglike_;

    if (!setup.is_power_post())
    {
        param.update_draws(smpl_count_);
        dataset.update_draws(smpl_count_);
        augm.update_draws(smpl_count_);
    }

    ++smpl_count_;
}

/* Clear posterior sampling. */
void MCMC::Impl::clear_posterior_sample_()
{
    smpl_count_ = 0;
    smpl_full_loglike_.clear();
    smpl_full_loglike_.shrink_to_fit();

    dataset.clear_draws();
    param.clear_draws();
    augm.clear_draws();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute new loglikelihood, store old value in buffer. Same for households. */
void MCMC::Impl::update_full_loglike_(bool cache_depen)
{
    const bool is_augm_ll = augm.is_any_model();
    double new_ll = 0.0;
    std::atomic<bool> has_error(false);

    // buffer infectivity profile if lookup tables are active
    infec_profile.update(cache_depen);

#pragma omp parallel for reduction(+ : new_ll)
    for (size_t R = 0; R < dataset.n_restricted_houses(); ++R)
    {
        Household &house = dataset.restricted(R);

        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            new_ll += update_house_loglike_(house, is_augm_ll, cache_depen);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    error_message = "While computing loglikelihood for house \"" +
                                    house.house_id() + "\":\n" +
                                    std::string(e.what());
            }
        }
    }
    if (!error_message.empty())
        throw std::runtime_error(error_message);

    set_full_loglike_(new_ll);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update and return contribution of given household to full logliklehood.
 * NOTE: house multiplicity doesn't work with data augmentation
 */
double MCMC::Impl::update_house_loglike_(Household &house, bool is_augm_ll, bool is_cache_depen)
{
    // buffer all household cache dependencies
    if (is_cache_depen)
    {
        house.buffer_all_cache();
        for (Individual &member : house.infected())
            update_cache_dependecies_(house, member);
    }

    // update household observation loglikelihood
    double output = compute_and_set_household_obs_loglike_(house);

    // update individual augmented loglikelihood
    if (is_augm_ll)
        for (Individual &member : house.infected())
            output += member.compute_and_set_augm_model_loglike();

    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : parameter move  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
bool MCMC::Impl::move_single_parameter_(const std::string &p, bool cache_depen)
{
    // proceed only if the new value is consistent
    // otherwise automatically refuse move
    if (propose_valid_param_move_(p))
    {
        // if parameter is cache dependent, update cache
        update_full_loglike_(cache_depen);

        // perform Metropolis-Hastings to accept or not move
        if (accept_valid_param_move_(p))
            return true; // accepted move
        else
            recover_parameter_(p, cache_depen);
    }
    return false; // refused move
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Move all model parameters. */
void MCMC::Impl::move_all_()
{
    // move augmented variables if any
    move_all_augm_();

    // move all simple model parameters
    for (const std::string &p : param.list_sample_standard)
        move_single_parameter_(p, false);

    // move parameters with cache dependecies
    for (const std::string &p : param.list_sample_cache)
        move_single_parameter_(p, true);
}

/* Move all augmented variables and model parameters and adapt their rule arguments. */
void MCMC::Impl::move_all_and_adapt_params_(unsigned long step)
{
    // move augmented variables if
    move_all_augm_();

    // move all model parameters and adapt their rule arguments
    const double ratio = static_cast<double>(step) / setup.n_adapt;

    // UPDATE: loops can be improved with viewers

    // move all simple model parameters
    for (const std::string &p : param.list_sample_standard)
    {
        Variable &this_param = param.at(p);
        this_param.update_accept_rate(move_single_parameter_(p, false));
        this_param.adapt_rule(ratio);
    }

    // move parameters with cache dependecies
    for (const std::string &p : param.list_sample_cache)
    {
        Variable &this_param = param.at(p);
        this_param.update_accept_rate(move_single_parameter_(p, true));
        this_param.adapt_rule(ratio);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run markov chain  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::Impl::run_monte_carlo_()
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    unsigned long &chain = setup.chain; // current chain
    init_posterior_sample_();           // initialize diagnostics

    // run multiple MCMC chains&
    for (chain = 1; chain <= setup.n_chains; ++chain)
    {
        // update progressbar text
        setup.prog_set_text();      // update progress bar text
        param.init();               // initialize augmented variables
        augm.init();                // initialize augmented variables
        update_full_loglike_(true); // init cache + full loglikelihood
        init_move_counts_();        // init move counts for adaptive

        // (Finite) Adaptative Iterations
        for (unsigned long i = 0; i < setup.n_adapt; i += setup.thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.thin; ++t)
                move_all_and_adapt_params_(i + t);

        // Burn-in Iterations
        for (unsigned long i = 0; i < setup.n_burnin; i += setup.thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.thin; ++t)
                move_all_();

        // Posterior sampling
        for (unsigned long i = setup.n_burnin; i < setup.n_steps; i += setup.thin, setup.prog_update())
        {
            for (unsigned long t = 1; t <= setup.thin; ++t)
                move_all_();
            update_posterior_sample_(); // store sample
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC::Impl : evaluate loglikelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize model and compute the loglikelihood. */
double MCMC::Impl::evaluate()
{
    try
    {
        setup.n_chains = 0;             // use 0 chains as a state variable
        assert_proper_dataset_();       // check if dataset is empty
        setup.rand_init_params = false; // UPDATE: this should also be true..
        param.init();                   // check no undefined parameters
        augm.init();                    // check no undefined augmented variables
        init_model_for_inference_();    // initialize model for inference
        update_full_loglike_(true);     // initialize cache + full loglikelihood

        return full_loglike_;
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return HSP_NA_DOUBLE;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC::Impl : inference  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
void MCMC::Impl::infer()
{
    try
    {
        setup.prog_declare();        // declare current progress bar
        init_monte_carlo_();         // initialize mcmc sampling
        init_model_for_inference_(); // user initialization model

        setup.prog_open(setup.c_n_sweeps * setup.n_chains); // open progress bar
        run_monte_carlo_();                                 // run Monte Carlo
        setup.prog_close();                                 // close progress bar
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        tools::error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : power posteriors  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
std::vector<dbl_vec> MCMC::Impl::power_posteriors_(const std::string &scheme, unsigned long n_rungs)
{
    try
    {
        //==========================//
        //===== INITIALIZATION =====//
        //==========================//

        const str_vec options = {"adaptive", "power-fractions"};
        if (!(options % scheme))
            throw std::invalid_argument("Unknown power posterior scheme \"" + scheme + "\"." +
                                        tools::str::find_closest(scheme, options, 3));

        init_monte_carlo_(); // initialize mcmc sampling

        // fix internal number of rungs (for progress bar)
        if (n_rungs >= 2)
            setup.n_rungs = n_rungs;
        else
            throw std::invalid_argument("N. of rungs must be greater than 1.");
        unsigned long &rung = setup.rung; // current rung

        // declare current progress bar
        setup.prog_declare();

        // initialize model for inference
        init_model_for_inference_();

        // implementation setup aliases
        double &t = setup.inv_T; // current inverse temperature

        //----- Buffer user setup -----//
        double buffer_inv_T = t;
        bool buffer_rand_init_params = setup.rand_init_params;

        // thermodynamic ladder estimates
        std::pair<double, double> stats = {}; // results of stats computation
        dbl_vec itm = {};                     // inverse temperature
        dbl_vec avg = {};                     // mean loglikelihood
        dbl_vec var = {};                     // variance loglikelihood

        //=================================//
        //===== INVERSE TEMPERATURE 0 =====//
        //=================================//

        setup.prog_open((setup.c_n_sweeps * setup.n_chains) * (setup.n_rungs + 1));

        // UPDATE: should I use previous inference for itm = 1 if available?
        // UPDATE: always under-estimate time if starting from itm = 0

        setup.rand_init_params = true;

        rung = 0;
        t = 0.0;
        run_monte_carlo_();
        stats = smpl_full_loglike_stats_();

        itm.push_back(t);
        avg.push_back(stats.first);
        var.push_back(stats.second);

        setup.rand_init_params = false;

        if (scheme == "adaptive")
        {
            //===========================//
            //===== ADAPTIVE LADDER =====//
            //===========================//

            //----- Inverse Temperature 1 -----//
            rung = 1;
            t = 1.0;
            run_monte_carlo_();
            stats = smpl_full_loglike_stats_();

            itm.push_back(t);
            avg.push_back(stats.first);
            var.push_back(stats.second);

            //----- Inverse Temperatures (0, 1) -----//
            size_t k = 0; // index to interval lower bound
            bool is_right = false;
            while (rung < n_rungs)
            {
                ++rung;

                //----- adaptive choice of inverse temperature -----//
                // intersection between tangents
                double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
                // otherwise weighted average
                if ((new_t <= itm[k]) || (new_t >= itm[k + 1]))
                    new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

                //----- power posterior sampling -----//
                t = new_t;
                run_monte_carlo_();
                stats = smpl_full_loglike_stats_();

                //----- find larger rectangular discretization error -----//
                is_right = (itm[k + 1] - t) * (avg[k + 1] - stats.first) > (new_t - itm[k]) * (stats.first - avg[k]);

                //----- insert in-between for all lists -----//
                itm.insert(itm.begin() + k + 1, t);
                avg.insert(avg.begin() + k + 1, stats.first);
                var.insert(var.begin() + k + 1, stats.second);

                // if right interval is chosen, index points to inserted as lower bound
                k += (is_right);
            }
        }
        else if (scheme == "power-fractions")
        {
            //==================================//
            //===== POWER FRACTIONS LADDER =====//
            //==================================//

            //----- Inverse Temperatures (0, 1] -----//
            for (rung = 1; rung <= n_rungs; ++rung)
            {
                // inverse tempearture = (rung / #rungs)^5
                t = std::pow(static_cast<double>(rung) / n_rungs, 5);
                run_monte_carlo_();
                stats = smpl_full_loglike_stats_();

                itm.push_back(t);
                avg.push_back(stats.first);
                var.push_back(stats.second);
            }
        }

        //----- Close progress bar -----//
        setup.prog_close();

        //----- Restore user setup -----//
        t = buffer_inv_T;
        setup.rand_init_params = buffer_rand_init_params;

        return {itm, avg, var};
    }
    catch (const std::exception &e)
    {
        setup.prog_crash();
        throw e;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : impl : simulation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary method to initialize simulations. */
void MCMC::Impl::aux_init_simulate_()
{
    assert_proper_dataset_();       // check if dataset is empty
    setup.rand_init_params = false; // disable random parameter initialization
    setup.n_chains = 0;             // use 0 chains as a state variable
    setup.n_threads = 1;            // impose no parallelization
    setup.set_n_threads();          // change accordingly
    param.init();                   // check no undefined parameters
    init_model_for_simulation_();   // initialize model for simulation
}

#ifdef _HOUSEPI_RCPP
//-- Perform a unique repetition. --//
Rcpp::DataFrame MCMC::Impl::simulate_repetition()
{
    using namespace Rcpp;

    const size_t n_cols = dataset.n_individuals();

    int_vec df_status = int_vec(n_cols);
    std::vector<dbl_vec> df_augm(augm.list_export.size(), dbl_vec(n_cols));

    size_t row = 0;
    for (Household house : dataset) // copies of households
        try
        {
            simulate_single_household_(house); // inplace simulation
            for (const Individual &member : house)
            {
                df_status[row] = member.status();
                for (size_t i = 0; i < augm.list_export.size(); ++i)
                    df_augm[i][row] = member.augm(augm.list_export[i]);
                ++row;
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While simulating house \"" + house.house_id() + "\":\n" + std::string(e.what()));
        }

    List output = List::create(_["status"] = wrap(std::move(df_status)));
    for (size_t i = 0; i < augm.list_export.size(); ++i)
        output[augm.list_export[i]] = wrap(std::move(df_augm[i]));
    return DataFrame(output);
}

//-- Perform simulation and return simulated dataset. --//
Rcpp::List MCMC::Impl::simulate(size_t n_rep)
{
    using namespace Rcpp;
    std::atomic<bool> has_error(false);

    try
    {
        // initialization
        aux_init_simulate_();
        setup.prog_open(dataset.n_houses() * n_rep);

        // run simulations
        std::vector<DataFrame> results(n_rep);

#pragma omp parallel for
        for (size_t rep = 0; rep < n_rep; ++rep)
        {
            if (has_error.load(std::memory_order_relaxed))
                continue;
            try
            {
                results[rep] = simulate_repetition();
            }
            catch (const std::exception &e)
            {
#pragma omp critical
                {
                    if (!has_error.exchange(true))
                        error_message = e.what();
                }
            }
        }

        // convert to Rcpp List
        List output(n_rep);
        for (size_t rep = 0; rep < n_rep; ++rep)
            output[rep] = std::move(results[rep]);

        // close progress bar
        setup.prog_close();

        return output;
    }
    catch (std::exception &e)
    {
        tools::error(e.what());
    }

    return List::create();
}

Rcpp::List MCMC::Impl::evidence(const std::string &scheme, unsigned long n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors_(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(_["estimate"] = output,
                            _["scheme"] = scheme,
                            _["n_rungs"] = n_rungs,
                            _["inv_temper"] = wrap(inv_tem),
                            _["mean_loglike"] = wrap(avg_vec),
                            _["var_loglike"] = wrap(var_vec));
    }
    catch (const std::exception &e)
    {
        tools::error("While estimating evidence with power posteriors:\n" +
                     std::string(e.what()));
    }
    return List();
}

#else
/* Print the simulated dataset to the provided file. */
void MCMC::Impl::print_simulated_dataset(const std::string &pathtofile)
{
    try
    {
        //----- initialization -----//
        setup.prog_declare("simulate");
        dataset.print_open(pathtofile);
        aux_init_simulate_();

        setup.prog_open(dataset.n_houses()); // open progress bar

        //----- run simulation -----//
        // note: deep copy of the household
        for (Household house : dataset)
            try
            {
                simulate_single_household_(house); // inplace simulation
                dataset.print_update(house);       // print the simulated dataset
                setup.prog_update();               // update progress bar
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating house \"" + house.house_id() + "\":\n" + e.what());
            }

        setup.prog_close();    // close progress bar
        dataset.print_close(); // finalize output
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        tools::error(e.what());
    }
}

/* Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
void MCMC::Impl::print_evidence(const std::string &scheme, unsigned long n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors_(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        // print result to stdout
        std::cout << output << std::endl;
    }
    catch (const std::exception &e)
    {
        tools::error("While estimating evidence with power posteriors:\n" +
                     std::string(e.what()));
    }
}
#endif
