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

/* Return the time of infection of a patient. */
double MCMC::get_infec_time(int h, int i) const
{
    UNUSED(h);
    UNUSED(i);
    return INF_DOUBLE;
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void MCMC::augm_store_dependecies(int h, int i)
{
    UNUSED(h);
    UNUSED(i);
}

/* Compute the loglikelihood for a single household. */
double MCMC::compute_house_loglike(int h)
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

/* Compute global loglikelihood */
long double MCMC::compute_global_loglike()
{
    long double output = 0.0;

#pragma omp parallel for reduction(+ : output)
    for (int h = 0; h < num_househ; h++)
    {
        // UPDATE : household loglikelihood is not long double, just double
        // to modify this traces to parameter values,
        // i. e. all program doubles to be cast to long double
        househ[h].set_loglike(compute_house_loglike(h));
        output += househ[h].get_loglike();
    }

    return output;
}

/* Set accept rates to zero. */
void MCMC::init_move_counts()
{
    // for parameters
    for (const std::string &parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }

    // for augmented variables
    for (const std::string &augmID : augm_to_be_inferred)
    {
        for (int h = 0; h < num_househ; h++)
        {
            househ[h].reset_augm_accept_rate(augmID);
        }
    }
}

/* Update household first case and first infection time.
Note: this structure assumes a single first infector.
*/
void MCMC::update_first_case(int h)
{
    double t_infec;
    double new_t0 = INF_DOUBLE; // NOTE: init to impossible values
    househ[h].first_case = -1;  // NOTE: init to impossible values

    for (const int &i : househ[h].infectious_set)
    {
        t_infec = get_infec_time(h, i);
        new_t0 = std::min(t_infec, new_t0);
        if (t_infec == new_t0)
        {
            househ[h].first_case = i;
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
void MCMC::move_single_param(std::string parID)
{
    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (propose_valid_param_move(parID))
    {
        // Compute the new gobla loglikelihood
        // set loglike in households and store old values in buffer
        long double new_global_loglike = compute_global_loglike();

        //===== METROPOLIS-HASTINGS =====//
        long double logQ = setup_inv_T * (new_global_loglike - get_global_loglike());
        logQ += get_move_logratio(); // logratio priors - logratio proposal

        if (uniform() < exp(logQ))
        {
            accept_param_move(parID);
            set_global_loglike(new_global_loglike);
        }
        else
        {
            refuse_param_move(parID);
            // restore buffered loglikelihood in households
#pragma omp parallel for default(shared)
            for (int h = 0; h < num_househ; h++)
            {
                househ[h].refuse();
            }
        }
    }
}

/* Sweep over all model parameters. */
void MCMC::move_all_parameters()
{
    for (const std::string &parID : param_to_be_inferred)
    {
        move_single_param(parID);
    }
}

/* Sweep over all model parameters plus adaptive step. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / setup_num_adapt;

    for (const std::string &parID : param_to_be_inferred)
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

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : augmented move  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented variable for all individuals */
void MCMC::init_augmented_variables(bool is_init)
{
    if (is_init)
    {
        //===== random initialization on demand =====//
        if (!augm_rand_init_set.empty())
        {
            for (int h = 0; h < num_househ; h++)
            {
                for (int i = 0; i < househ[h].house_size; i++)
                {
                    for (const std::string &augmID : augm_rand_init_set)
                    {
                        househ[h].indiv[i].init_augm(augmID, get_prng());
                    }
                }

                // associate a first case (if any)
                update_first_case(h);
                // initialization make sense only if a case is present
                if (househ[h].first_case != -1)
                {
                    for (int i = 0; i < househ[h].house_size; i++)
                    {
                        augm_store_dependecies(h, i);
                    }
                }
            }
        }
    }
}

/* Update augmented variable `augmID` for individual `i` (). */
void MCMC::augm_update(int h, int i, std::string augmID, double new_val)
{
    househ[h].indiv[i].set_augm(augmID, new_val);
    update_first_case(h);         // Update information of first cases
    augm_store_dependecies(h, i); // Update infectivity profiles related to idx
}

/* Propose augmented variable `augmID` for individual `i`. */
double MCMC::augm_propose(int h, int i, std::string augmID)
{
    return househ[h].indiv[i].augm_propose(augmID, get_prng());
}

/* Move augmented variables for infectious members of a single household.
The move is performed one member at a time, one augmented variable at at time.
*/
long double MCMC::move_single_house_augm(int h)
{
    double old_value, new_value;
    long double logQ, old_house_loglike, new_house_loglike;
    bool is_accepted = false;

    long double delta_global_loglike = 0.0;

    // Loop over infectious residents of the househould
    // FIXME : only infectious are subject to data augmentation?
    for (const int &i : househ[h].infectious_set)
    {
        // Loop over augmented variables of infectious
        for (const std::string &augmID : augm_to_be_inferred)
        {
            is_accepted = false;
            old_value = househ[h].indiv[i].get_augm(augmID); // store old value
            old_house_loglike = househ[h].get_loglike();     // store old loglike
            new_value = augm_propose(h, i, augmID);          // propose new value

            // proceed only if the new value is consistent
            if (househ[h].indiv[i].augm_valid(augmID, new_value))
            {
                // update value to perform computations
                augm_update(h, i, augmID, new_value);
                new_house_loglike = compute_house_loglike(h);

                // MetropolisHasting contribution
                logQ = setup_inv_T * (new_house_loglike - old_house_loglike);
                logQ -= househ[h].indiv[i].augm_propose_logratio(augmID, old_value, new_value);

                if (uniform() < exp(logQ))
                {
                    is_accepted = true;
                    househ[h].set_loglike(new_house_loglike);
                    delta_global_loglike += new_house_loglike - old_house_loglike;
                }
                else
                {
                    augm_update(h, i, augmID, old_value); // restore old value
                }
            }
            househ[h].update_augm_accept_rate(augmID, is_accepted);
        }
    }

    return delta_global_loglike;
}

/* Sweep over all households to move augmented variables. */
void MCMC::move_all_augm_variables()
{
    if (__is_data_augm)
    {
        long double variation = 0.;

#pragma omp parallel for reduction(+ : variation)
        for (int h = 0; h < num_househ; h++)
        {
            variation += move_single_house_augm(h);
        }
        set_global_loglike(get_global_loglike() + variation);
    }
}

/* Return the value of the acceptance rate of the given augmented variable */
double MCMC::get_augm_accept_rate(std::string augmID)
{
    int accept_count = 0, propos_count = 0;

#pragma omp parallel for reduction(+ : accept_count, propos_count)
    for (int h = 0; h < num_househ; h++)
    {
        accept_count += househ[h].get_augm_accept_counts(augmID);
        propos_count += househ[h].get_augm_propos_counts(augmID);
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
long double MCMC::evaluate()
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
    // diagnostics for all parameters
    init_param_diagnostics(dgn_ss);

// diagnostics over households
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].init_diagnostics(dgn_ss);
    }

    dgn_avg = 0.; // (global) loglikelihood
    dgn_sqr = 0.; // mean squared (global) loglikelihood
}

/* Update Bayesian diagnostics */
void MCMC::update_diagnostics()
{
    if (!__is_power_post)
    {
        // diagnostics for all parameters
        update_param_diagnostics();

        // diagnostics over households
#pragma omp parallel for default(shared)
        for (int h = 0; h < num_househ; h++)
        {
            househ[h].update_diagnostics();
        }
    }

    // average (global) loglikelihood
    dgn_avg += global_loglike / dgn_ss;

    // average squared (global) loglikelihood
    dgn_sqr += std::pow(global_loglike, 2) / dgn_ss;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run markov chain  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The text to print in the Markov chian progress bar. */
std::string MCMC::get_progtxt(std::string txt)
{
    std::ostringstream oss;
    if (__is_power_post)
    {
        oss << "power post " << __rung << "/" << __num_rungs;
    }
    else if (setup_num_chains > 1)
    {
        oss << txt << " " << __chain << "/" << setup_num_chains;
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
    for (int h = 0; h < num_househ; h++)
    {
        for (int i = 0; i < househ[h].house_size; i++)
        {
            househ[h].indiv[i].disable_assert();
        }
    }
}

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed mcmc sweeps. */
void MCMC::run_monte_carlo(bool progbar, ulong sweep, bool is_init_augm)
{
    using namespace screen;

    //==========================//
    //===== Initialization =====//
    //==========================//

    //----- Initialize model -----//
    init_model_for_inference();

    //----- Initialize diagnostics -----//
    init_diagnostics();

    //===== Run Multiple MCMC Chains =====//
    for (__chain = 1; __chain <= setup_num_chains; __chain++)
    {
        //-----Initialize parameter inference-----//
        init_rand_params(setup_rnd_init_params);

        //-----Initialize augmented variables with model-----//
        init_augmented_variables(is_init_augm);

        //-----Initialize chain-----//
        init_single_chain();

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
            for (ulong t = 1; t <= setup_thinning; t++)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(i + t);
            }
            sweep += setup_thinning;
            progress_bar_run(get_progtxt("adapt"), sweep, __tot_num_sweeps, progbar);
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
            for (ulong t = 1; t <= setup_thinning; t++)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            sweep += setup_thinning;
            progress_bar_run(get_progtxt("burn"), sweep, __tot_num_sweeps, progbar);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        for (ulong i = setup_burn_in; i < setup_num_iter; i += setup_thinning)
        {
            // thinning mcmc moves
            for (ulong t = 1; t <= setup_thinning; t++)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            sweep += setup_thinning;

            // compute diagnostics for the previous iterations
            update_diagnostics();

            // print inference state
            inferstate_print();

            progress_bar_run(get_progtxt("sample"), sweep, __tot_num_sweeps, progbar);
        }
    }

    close_diagnostics();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : infer  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
// FIXME : output_file shouldn't be necessary
// FIXME : store diagnostics instead of canceling everything for each chain (~)
// FIXME : Gelman Rubin could be integrated
void MCMC::infer(std::string output_file, bool progbar)
{
    //----- Check if data is meaningful -----//
    if (report_data["n.miss.stat"] == num_indiv)
    {
        screen::raise_error("No infectious status has been provided in the data file.");
    }

    //----- total number of mcmc sweeps -----//
    __tot_num_sweeps = __mcmc_sweeps;

    //----- Initialize output file-----//
    inferstate_init(output_file);

    //----- Start progress bar and chrono -----//
    screen::progress_bar_init("init", progbar);

    //----- Run Markov Chain -----//
    run_monte_carlo(progbar);

    //----- Close progress bar -----//
    screen::progress_bar_end("end", progbar);

    //----- Close output file-----//
    inferstate_close();
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
        screen::raise_error("Number of rungs must be greater than 1.");
    }
    else
    {
        // fix internal number of rungs (for progress bar)
        __num_rungs = num_rungs;
    }

    //----- Backup user setup to restore them later -----//
    bool bkp_rnd_init_params = setup_rnd_init_params;
    long double bkp_inv_T = setup_inv_T;

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    std::vector<long double> inv_t = {};  // inverse temperature
    std::vector<long double> avg_ll = {}; // mean loglikelihood
    std::vector<long double> var_ll = {}; // variance loglikelihood

    //----- Fix internal states -----//
    __is_print_inferstate = false; // avoid printing
    __is_infer_complete = false;   // showing results becomes meaningless
    __is_power_post = true;

    //----- Total number of mcmc sweeps -----//
    __tot_num_sweeps = __mcmc_sweeps * (__num_rungs + 1);

    //----- Start progress bar and chrono -----//
    screen::progress_bar_init("init", progbar);

    //==================================//
    //===== INVERSE TEMPEREATURE 0 =====//
    //==================================//

    setup_rnd_init_params = true;
    __rung = 0;
    setup_inv_T = 0.;
    run_monte_carlo(progbar, 0, true);
    inv_t.push_back(setup_inv_T);
    avg_ll.push_back(dgn_avg);
    var_ll.push_back(dgn_var());
    setup_rnd_init_params = false;

    // thermodynamic ladder: use modified power fractions or adaptive scheme
    if (power_fractions)
    {
        //==================================//
        //===== POWER FRACTIONS LADDER =====//
        //==================================//

        //----- Inverse Temperatures (0, 1] -----//
        for (__rung = 1; __rung < __num_rungs; __rung++)
        {
            setup_inv_T = std::pow(static_cast<long double>(__rung) / (__num_rungs - 1), 5);
            run_monte_carlo(progbar, __rung * __mcmc_sweeps, false);
            inv_t.push_back(setup_inv_T);
            avg_ll.push_back(dgn_avg);
            var_ll.push_back(dgn_var());
        }
    }
    else
    {
        //===========================//
        //===== ADAPTIVE LADDER =====//
        //===========================//

        //----- Inverse Temperature 1 -----//
        __rung = 1;
        setup_inv_T = 1.;
        run_monte_carlo(progbar, __mcmc_sweeps, true);
        inv_t.push_back(setup_inv_T);
        avg_ll.push_back(dgn_avg);
        var_ll.push_back(dgn_var());

        //----- Inverse Temperatures (0, 1) -----//
        size_t k = 0; // index to interval lower bound
        while (__rung < __num_rungs)
        {
            __rung++;

            //----- adaptive choice of inverse temperature -----//
            // intersection between tangents
            setup_inv_T = (avg_ll[k + 1] - inv_t[k + 1] * var_ll[k + 1] - avg_ll[k] + inv_t[k] * var_ll[k]) / (var_ll[k] - var_ll[k + 1]);

            // otherwise weighted average
            if ((setup_inv_T <= inv_t[k]) || (setup_inv_T >= inv_t[k + 1]))
            {
                setup_inv_T = inv_t[k] + var_ll[k + 1] * (inv_t[k + 1] - inv_t[k]) / (var_ll[k] + var_ll[k + 1]);
            }

            //----- power posterior sampling -----//
            run_monte_carlo(progbar, __rung * __mcmc_sweeps, false);

            //----- find larger rectangular discretization error -----//
            is_right = (inv_t[k + 1] - setup_inv_T) * (avg_ll[k + 1] - dgn_avg) > (setup_inv_T - inv_t[k]) * (dgn_avg - avg_ll[k]);

            //----- insert in-between for all lists -----//
            inv_t.insert(inv_t.begin() + k + 1, setup_inv_T);
            avg_ll.insert(avg_ll.begin() + k + 1, dgn_avg);
            var_ll.insert(var_ll.begin() + k + 1, dgn_var());

            // if right interval is chosen, index points to inserted as lower bound
            k += (is_right);
        }
    }

    //----- trapezoidal rule to compute marginal likelihood -----//
    long double d_inv_t, output = 0.;
    for (int k = 0; k < num_rungs; k++)
    {
        d_inv_t = inv_t[k + 1] - inv_t[k];
        output += (avg_ll[k + 1] + avg_ll[k]) * d_inv_t / 2 - (var_ll[k + 1] - var_ll[k]) * std::pow(d_inv_t, 2) / 12;
    }

    //----- Close progress bar -----//
    screen::progress_bar_end("end", progbar);

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

/* Compute diagnostics at the end of inference. */
void MCMC::close_diagnostics()
{
    if (!__is_power_post)
    {
        //---- compute WAIC and implement warning -----//
        // Ref: Vhetari et al. (2016) Stat Comput.

        const long double p_waic_limit = 0.4;
        ulong count_p_waic_large = 0;

        // Loop over households
#pragma omp parallel for reduction(+ : dgn_lppd, dgn_eff_pars, count_p_waic_large)
        for (int h = 0; h < num_househ; h++)
        {
            long double p_waic = househ[h].dgn_var();

            dgn_lppd += std::log(househ[h].dgn_like);
            dgn_eff_pars += p_waic;
            count_p_waic_large += (p_waic > p_waic_limit);
        }

        //----- empirical WAIC warning -----//
        if (count_p_waic_large > 0)
        {
            std::ostringstream oss;
            oss << screen::percent(count_p_waic_large, num_househ)
                << " p_waic estimates greater than " << p_waic_limit << "."
                << " Unreliable WAIC (Vehtari-Gelman-Gabry criterion).";
            screen::raise_warning(oss.str(), Large_p_WAIC_Warning);
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

#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        simulate_single_house(h);
    }

    update_report_data();
}