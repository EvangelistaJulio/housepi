/*
 * mcmc.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Large_p_WAIC_Warning = false; // unreliable WAIC

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (user defined)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void MCMC::augm_store_dependecies(int h, int i)
{
    UNUSED(h);
    UNUSED(i);
}

/* Compute the loglikelihood for a single household. */
double MCMC::compute_house_loglike(int h) const
{
    UNUSED(h);
    return 0.;
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : methods (general)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* House loglikelihood including augmented modl contribution. */
double MCMC::compute_augm_loglike(int h) const
{
    double hh_loglike = 0.;

    // Contribution of model associated to augmented variables
    for (const string &augmID : augm_with_model_set)
    {
        for (int i : house[h].infectious_set)
        {
            hh_loglike += house[h].indiv[i].augm_model_logpdf(augmID);
            // FIXME add augm_model_lognorm
        }
    }

    return hh_loglike;
}

/* Compute global loglikelihood */
double MCMC::compute_global_loglike()
{
    double output = 0.0;

#pragma omp parallel for reduction(+ : output)
    for (int h = 0; h < num_house; ++h)
    {
        house[h].set_loglike(compute_house_loglike(h));
        output += house[h].get_loglike() + compute_augm_loglike(h);
    }

    return output;
}

/* Set accept rates to zero. */
void MCMC::init_move_counts()
{
    // for parameters
    for (const string &parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }

    // for augmented variables
    for (const string &augmID : augm_to_be_inferred)
    {
        for (int h = 0; h < num_house; ++h)
        {
            house[h].reset_augm_accept_rate(augmID);
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
//////////////////////////////////////  MCMC : parameter move  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
void MCMC::move_single_param(const string &parID)
{
    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (propose_valid_param_move(parID))
    {
        // Compute the new gobla loglikelihood
        // set loglike in households and store old values in buffer
        double new_global_loglike = compute_global_loglike();

        //===== Metropolis-Hasting algorithm =====//

        double logQ = setup_inv_T * (new_global_loglike - get_global_loglike());
        logQ += get_move_logratio(); // logratio priors - logratio proposal

        if (uniform() < exp(logQ))
        {
            accept_param_move(parID);
            set_global_loglike(new_global_loglike);
        }
        else
        {
            // restore buffered parameter valie
            refuse_param_move(parID);

            // restore buffered loglikelihood in households
#pragma omp parallel for default(shared)
            for (int h = 0; h < num_house; ++h)
            {
                house[h].refuse();
            }
        }
    }
}

/* Move all model parameters. */
void MCMC::move_all_parameters()
{
    for (const string &parID : param_to_be_inferred)
    {
        move_single_param(parID);
    }
}

/* Move all model parameters and adapt their rule arguments. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / setup_num_adapt;

    for (const string &parID : param_to_be_inferred)
    {
        move_single_param(parID);
        adapt(parID, adapt_ratio);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Change the status and add augmented variable if infectious.
  NOTE: this approach slows down enormously. // FIXME
*/
void MCMC::change_status_and_augm(int h, int i, int new_status)
{
    if (house[h].change_status(i, new_status))
    {
        if (_is_any_augm_)
        {
            if (house[h].is_infectious(i))
            {
                assign_augm(h, i);
            }
            /*
            // note : this is probably not necessary.
            else
            {
                clear_augm(h, i);
            }
            */
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented variable for all individuals */
void MCMC::init_augmented_variables()
{
    for (int h = 0; h < num_house; ++h)
    {
        // initialization of augmented variables for infectious
        for (int i : house[h].infectious_set)
        {
            augm_store_dependecies(h, i);
        }
    }
}

/* Move augmented variables for infectious members of a single household.
The move is performed one member at a time, one augmented variable at at time.
// FIXME adjust with buffer as for param.
*/
double MCMC::move_single_house_augm(int h)
{
    double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    // FIXME : truncated model_logpdf are not taken into account
    // FIXME : Gibbs sampling should be always valid but it's not now
    for (const int &i : house[h].infectious_set)
    {
        // Loop over augmented variables of infectious
        for (const string &augmID : augm_to_be_inferred)
        {
            bool is_accepted = false;
            double old_hh_ll = house[h].get_loglike();

            // propose new value and store old in buffer
            house[h].indiv[i].augm_propose(augmID, get_prng());

            // proceed only if the new value is consistent
            if (house[h].indiv[i].augm_valid(augmID))
            {
                // update dependencies to perform computations
                augm_store_dependecies(h, i);
                double new_hh_ll = compute_house_loglike(h);

                //===== Metropolis-Hasting algorithm =====//

                // local change in loglikelihood
                double delta_ll = new_hh_ll - old_hh_ll;
                delta_ll += house[h].indiv[i].augm_model_logratio(augmID);

                // total acceptance ratio
                // FIXME : Gibbs sampling doesn't work here if inverse temperature...
                double logQ = setup_inv_T * delta_ll;
                logQ -= house[h].indiv[i].augm_propose_logratio(augmID);

                if (uniform() < exp(logQ))
                {
                    is_accepted = true;
                    house[h].set_loglike(new_hh_ll);
                    delta_global_loglike += delta_ll;
                }
                else
                {
                    // restore value and its dependencies
                    house[h].indiv[i].augm_restore_value(augmID);
                    augm_store_dependecies(h, i);
                }
            }
            else
            {
                house[h].indiv[i].augm_restore_value(augmID);
            }

            house[h].update_augm_accept_rate(augmID, is_accepted);
        }
    }

    return delta_global_loglike;
}

/* Move augmented variables from all households. */
void MCMC::move_all_augm_variables()
{
    if (_is_move_augm_)
    {
        double variation = 0.;

#pragma omp parallel for reduction(+ : variation)
        for (int h = 0; h < num_house; ++h)
        {
            variation += move_single_house_augm(h);
        }
        set_global_loglike(get_global_loglike() + variation);
    }
}

/* Return the value of the acceptance rate of the given augmented variable */
double MCMC::get_augm_accept_rate(const string &augmID)
{
    int accept_count = 0, propos_count = 0;

#pragma omp parallel for reduction(+ : accept_count, propos_count)
    for (int h = 0; h < num_house; ++h)
    {
        accept_count += house[h].get_augm_accept_counts(augmID);
        propos_count += house[h].get_augm_propos_counts(augmID);
    }

    return static_cast<double>(accept_count) / propos_count;
}
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : display likelihood  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method to only compute the loglikelihood and display it. */
double MCMC::evaluate()
{
    //===== Initialization =====//

    //----- Initialize model -----//
    init_model_for_inference();

    //-----random initialization of augmented variables with model-----//
    init_augmented_variables();

    return compute_global_loglike();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void MCMC::init_diagnostics()
{
    ulong sample_size = (setup_num_iter - setup_burn_in) / setup_thinning;
    sample_size *= setup_num_chains;

    dgn_count = 0;
    dgn_glob_loglike = vector<double>(sample_size, std::numeric_limits<double>::quiet_NaN());

    if (!_is_power_post_)
    {
        // diagnostics for all parameters
        init_param_diagnostics(sample_size);

// diagnostics over households
#pragma omp parallel for default(shared)
        for (int h = 0; h < num_house; ++h)
        {
            house[h].init_diagnostics(sample_size);
        }
    }
}

/* Update Bayesian diagnostics */
void MCMC::update_diagnostics()
{
    // diagnostics for global loglike
    dgn_glob_loglike[dgn_count++] = get_global_loglike();

    if (!_is_power_post_)
    {
        // diagnostics for all parameters
        update_param_diagnostics();

        // diagnostics over households
#pragma omp parallel for default(shared)
        for (int h = 0; h < num_house; ++h)
        {
            house[h].update_diagnostics();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run markov chain  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The text to print in the Markov chian progress bar. */
string MCMC::get_progtxt(const string &txt)
{
    std::ostringstream oss;
    if (_is_power_post_)
    {
        oss << "power post " << _rung_ << "/" << _num_rungs_;
    }
    else if (setup_num_chains > 1)
    {
        oss << txt << " " << _chain_ << "/" << setup_num_chains;
    }
    else
    {
        oss << txt;
    }
    return oss.str();
}

/* Disable all assert methods for user definitions. */
void MCMC::disable_assert()
{
    disable_assert_param();

    // disable augmented variable assertion for all individuals
    for (int h = 0; h < num_house; ++h)
    {
        for (int i = 0; i < house[h].house_size; ++i)
        {
            house[h].indiv[i].disable_assert();
        }
    }
}

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::run_monte_carlo(bool progbar, ulong sweep)
{
    using namespace basics;

    //==========================//
    //===== Initialization =====//
    //==========================//

    //----- Initialize model -----//
    init_model_for_inference();

    //----- Initialize diagnostics -----//
    init_diagnostics();

    //===== Run Multiple MCMC Chains =====//
    for (_chain_ = 1; _chain_ <= setup_num_chains; ++_chain_)
    {
        //-----Initialize parameter inference-----//
        init_rand_params(setup_rnd_init_params);

        //-----Initialize augmented variables with model-----//
        init_augmented_variables();

        // set the value of the global loglikelihood
        set_global_loglike(compute_global_loglike());

        // avoid asserting user definitions after first round
        // UPDATE : this can be improved, e. g. with a debug mode
        disable_assert();

        // initialize the acceptance counts
        init_move_counts();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        for (ulong i = 0; i < setup_num_adapt; i += setup_thinning)
        {
            for (ulong t = 1; t <= setup_thinning; ++t)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(i + t);
            }
            ++sweep;
            progress_bar_run(get_progtxt("adapt"), sweep, _tot_thin_sweeps_, progbar);
        }

        // reset the acceptance counts
        // UPDATE : move counts become useless since they are not even printed : improve.
        if (setup_num_adapt > 0)
        {
            init_move_counts();
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        for (ulong i = 0; i < setup_burn_in; i += setup_thinning)
        {
            for (ulong t = 1; t <= setup_thinning; ++t)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            ++sweep;
            progress_bar_run(get_progtxt("burn"), sweep, _tot_thin_sweeps_, progbar);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        for (ulong i = setup_burn_in; i < setup_num_iter; i += setup_thinning)
        {
            // thinning mcmc moves
            for (ulong t = 1; t <= setup_thinning; ++t)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            ++sweep;

            // compute diagnostics for the previous iterations
            update_diagnostics();
            progress_bar_run(get_progtxt("sample"), sweep, _tot_thin_sweeps_, progbar);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
// FIXME : Gelman Rubin could be integrated
void MCMC::infer(bool progbar)
{
    //----- Check if data is meaningful -----//
    if (report_data["n.miss.stat"] == num_indiv)
    {
        basics::raise_error("No infectious status has been provided in the data file.");
    }

    //----- total number of mcmc sweeps -----//
    _tot_thin_sweeps_ = _thin_sweeps_;

    //----- Start progress bar and chrono -----//
    basics::progress_bar_init("init", progbar);

    //----- Run Markov Chain -----//
    run_monte_carlo(progbar);

    //----- Close progress bar -----//
    basics::progress_bar_end("end", progbar);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : power posteriors  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Evaluate the logarithm of marginal likelihood with power posteriors. */
double MCMC::evaluate_marginal(int num_rungs, bool progbar, bool power_fractions)
{
    /* Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).*/

    if (num_rungs < 2)
    {
        basics::raise_error("Number of rungs must be greater than 1.");
    }
    else
    {
        // fix internal number of rungs (for progress bar)
        _num_rungs_ = num_rungs;
    }

    //----- Backup user setup to restore them later -----//
    bool bkp_rnd_init_params = setup_rnd_init_params;
    double bkp_inv_T = setup_inv_T;

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    double avg, var;
    vector<double> inv_t = {};      // inverse temperature
    vector<double> vec_avg_ll = {}; // mean loglikelihood
    vector<double> vec_var_ll = {}; // variance loglikelihood

    //----- Fix internal states -----//
    _is_power_post_ = true;

    //----- Total number of mcmc sweeps -----//
    _tot_thin_sweeps_ = _thin_sweeps_ * (_num_rungs_ + 1);

    //----- Start progress bar and chrono -----//
    basics::progress_bar_init("init", progbar);

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    setup_rnd_init_params = true;
    _rung_ = 0;
    setup_inv_T = 0.;
    run_monte_carlo(progbar, 0);
    basics::statistics(dgn_glob_loglike, avg, var);

    inv_t.push_back(setup_inv_T);
    vec_avg_ll.push_back(avg);
    vec_var_ll.push_back(var);

    // thermodynamic ladder: use modified power fractions or adaptive scheme
    if (power_fractions)
    {
        //==================================//
        //===== POWER FRACTIONS LADDER =====//
        //==================================//

        //----- Inverse Temperatures (0, 1] -----//
        setup_rnd_init_params = false;
        for (_rung_ = 1; _rung_ < _num_rungs_; ++_rung_)
        {
            setup_inv_T = std::pow(static_cast<long double>(_rung_) / (_num_rungs_ - 1), 5);
            run_monte_carlo(progbar, _rung_ * _thin_sweeps_);
            basics::statistics(dgn_glob_loglike, avg, var);

            inv_t.push_back(setup_inv_T);
            vec_avg_ll.push_back(avg);
            vec_var_ll.push_back(var);
        }
    }
    else
    {
        //===========================//
        //===== ADAPTIVE LADDER =====//
        //===========================//

        //----- Inverse Temperature 1 -----//
        setup_rnd_init_params = false;
        _rung_ = 1;
        setup_inv_T = 1.;
        run_monte_carlo(progbar, _thin_sweeps_);
        basics::statistics(dgn_glob_loglike, avg, var);

        inv_t.push_back(setup_inv_T);
        vec_avg_ll.push_back(avg);
        vec_var_ll.push_back(var);

        //----- Inverse Temperatures (0, 1) -----//
        size_t k = 0; // index to interval lower bound
        while (_rung_ < _num_rungs_)
        {
            ++_rung_;

            //----- adaptive choice of inverse temperature -----//
            // intersection between tangents
            setup_inv_T = (vec_avg_ll[k + 1] - inv_t[k + 1] * vec_var_ll[k + 1] - vec_avg_ll[k] + inv_t[k] * vec_var_ll[k]) / (vec_var_ll[k] - vec_var_ll[k + 1]);

            // otherwise weighted average
            if ((setup_inv_T <= inv_t[k]) || (setup_inv_T >= inv_t[k + 1]))
            {
                setup_inv_T = inv_t[k] + vec_var_ll[k + 1] * (inv_t[k + 1] - inv_t[k]) / (vec_var_ll[k] + vec_var_ll[k + 1]);
            }

            //----- power posterior sampling -----//
            run_monte_carlo(progbar, _rung_ * _thin_sweeps_);
            basics::statistics(dgn_glob_loglike, avg, var);

            //----- find larger rectangular discretization error -----//
            is_right = (inv_t[k + 1] - setup_inv_T) * (vec_avg_ll[k + 1] - avg) > (setup_inv_T - inv_t[k]) * (avg - vec_avg_ll[k]);

            //----- insert in-between for all lists -----//
            inv_t.insert(inv_t.begin() + k + 1, setup_inv_T);
            vec_avg_ll.insert(vec_avg_ll.begin() + k + 1, avg);
            vec_var_ll.insert(vec_var_ll.begin() + k + 1, var);

            // if right interval is chosen, index points to inserted as lower bound
            k += (is_right);
        }
    }

    //----- trapezoidal rule to compute marginal likelihood -----//
    double d_inv_t, output = 0.;
    for (int k = 0; k < num_rungs; ++k)
    {
        d_inv_t = inv_t[k + 1] - inv_t[k];
        output += (vec_avg_ll[k + 1] + vec_avg_ll[k]) * d_inv_t / 2. - (vec_var_ll[k + 1] - vec_var_ll[k]) * std::pow(d_inv_t, 2) / 12.;
    }

    //----- Close progress bar -----//
    basics::progress_bar_end("end", progbar);

    //----- Restore user setup -----//
    setup_inv_T = bkp_inv_T;
    setup_rnd_init_params = bkp_rnd_init_params;

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : methods (simulation)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the time of infection of a patient. */
void MCMC::simulate_single_house(int h)
{
    UNUSED(h);
}

/* Perform simulation as provided by the user */
void MCMC::simulate()
{
    //===== Initialization =====//

    //----- Initialize model -----//
    init_model_for_simulation();

    //-----random initialization of augmented variables with model-----//
    init_augmented_variables();

    bool progbar = true;
    basics::progress_bar_init("simul", progbar);

#ifdef _OPENMP
    int thin_num_house = num_house / get_num_threads();
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_house; ++h)
    {
        simulate_single_house(h);

        if (omp_get_thread_num() == 0)
        {
            basics::progress_bar_run("simul", h, thin_num_house, progbar);
        }
    }
#else
    for (int h = 0; h < num_house; ++h)
    {
        simulate_single_house(h);
        basics::progress_bar_run("simul", h, num_house, progbar);
    }
#endif

    basics::progress_bar_end("simul", progbar);

    update_report_data();
}