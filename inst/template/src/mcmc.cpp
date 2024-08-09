// mcmc.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Inconmplete_Infer_Diagnostics_Warning = false; // wrong call to show diagnostics
bool Large_p_WAIC_Warning = false;                  // unreliable WAIC

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
void MCMC::augm_var_store_dependecies(int h, int i)
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
    for (std::string parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }

    // for augmented variables
    for (std::string augmID : augm_to_be_inferred)
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
    int new_case0 = -1;         // NOTE: init to impossible values

    for (int i : househ[h].infectious_set)
    {
        t_infec = get_infec_time(h, i);
        new_t0 = std::min(t_infec, new_t0);
        if (t_infec == new_t0)
        {
            new_case0 = i;
        }
    }
    househ[h].first_case = new_case0;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : methods (parameter move)  /////////////////////////////////
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
    for (std::string parID : param_to_be_inferred)
    {
        move_single_param(parID);
    }
}

/* Sweep over all model parameters plus adaptive step. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / setup_num_adapt;

    for (std::string parID : param_to_be_inferred)
    {
        move_single_param(parID);
        adapt(parID, adapt_ratio);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : methods (augmented move)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Propose augmented variable `augmID` for individual `i`. */
double MCMC::augm_var_propose(int h, int i, std::string augmID)
{
    return househ[h].indiv[i].augm_propose(augmID, get_prng());
}

/* Update augmented variable `augmID` for individual `i` (). */
void MCMC::augm_var_update(int h, int i, std::string augmID, double new_val)
{
    househ[h].indiv[i].set_augm(augmID, new_val);
    update_first_case(h);             // Update information of first cases
    augm_var_store_dependecies(h, i); // Update infectivity profiles related to idx
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
    for (int i : househ[h].infectious_set)
    {
        // Loop over augmented variables of infectious
        for (std::string augmID : augm_to_be_inferred)
        {
            is_accepted = false;
            old_value = househ[h].indiv[i].get_augm(augmID); // store old value
            old_house_loglike = househ[h].get_loglike();     // store old loglike
            new_value = augm_var_propose(h, i, augmID);      // propose new value

            // proceed only if the new value is consistent
            if (househ[h].indiv[i].augm_valid(augmID, new_value))
            {
                augm_var_update(h, i, augmID, new_value); // update value to perform computations
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
                    augm_var_update(h, i, augmID, old_value); // restore old value
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

    //-----random initialization of augmented variables with model-----//
    init_augmented_variables();

    //----- Initialize model -----//
    init_model_for_inference();

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
    init_param_diagnostics(dgn_sample);

// diagnostics over households
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].init_diagnostics(dgn_sample);
    }

    dgn_avg_ll = 0.; // (global) loglikelihood
    dgn_sqr_ll = 0.; // mean squared (global) loglikelihood
}

/* Update Bayesian diagnostics */
void MCMC::update_diagnostics()
{
    // diagnostics for all parameters
    update_param_diagnostics();

// diagnostics over households
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].update_diagnostics();
    }

    // (global) loglikelihood
    dgn_avg_ll += global_loglike / dgn_sample;

    // mean squared (global) loglikelihood
    dgn_sqr_ll += std::pow(global_loglike, 2) / dgn_sample;
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

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed mcmc sweeps. */
void MCMC::run_monte_carlo(bool progbar, ulong sweep, bool is_init_augm)
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    std::string progtxt;

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

        // avoid asserting parameters definitions after first call
        disable_assert_param();

        // initialize the acceptance counts
        init_move_counts();

        //=========================================//
        //======== (Finite) Adaptative MCMC =======//
        //=========================================//

        progtxt = get_progtxt("adapt");

        for (ulong i = 0; i < setup_num_adapt; i += setup_thinning)
        {
            for (ulong t = 1; t <= setup_thinning; t++)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(i + t);
            }
            sweep += setup_thinning;
            screen::progress_bar_run(progtxt, sweep, __tot_num_sweeps, progbar);
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

        progtxt = get_progtxt("burn");

        for (ulong i = 0; i < setup_burn_in; i += setup_thinning)
        {
            for (ulong t = 1; t <= setup_thinning; t++)
            {
                move_all_augm_variables();
                move_all_parameters();
            }
            sweep += setup_thinning;
            screen::progress_bar_run(progtxt, sweep, __tot_num_sweeps, progbar);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        progtxt = get_progtxt("sample");

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

            screen::progress_bar_run(progtxt, sweep, __tot_num_sweeps, progbar);
        }
    }
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

/* Use Metropolis-Hastings MCMC to conmpute marginal likelihood with power posteriors. */
double MCMC::power_posteriors(int num_rungs, bool progbar)
{
    // Ref 1. Friel, N., Hurn, M. & Wyse, J. Improving power posterior estimation of statistical evidence. Stat Comput 24, 709–723 (2014).

    // FLAG
    progbar = false;
    // FLAG-END

    if (num_rungs < 2)
    {
        screen::raise_error(
            "Number of rungs must be greater than 1.");
    }
    else
    {
        // fix internal number of rungs (for progress bar)
        __num_rungs = num_rungs;
    }

    //----- Backup user setup to restore them later -----//
    bool bkp_rnd_init_params = setup_rnd_init_params;
    long double bkp_inv_T = setup_inv_T;

    //----- Fix internal states -----//
    __is_print_inferstate = false; // avoid printing
    __is_infer_complete = false;   // showing results becomes meaningless
    __is_power_post = true;

    //----- Total number of mcmc sweeps -----//
    __tot_num_sweeps = __mcmc_sweeps * (__num_rungs + 1);

    //----- Thermodynamic ladder estimates -----//
    bool is_right = false;
    std::vector<long double> inv_t = {};  // inverse temperature
    std::vector<long double> avg_ll = {}; // mean loglikelihood
    std::vector<long double> var_ll = {}; // variance loglikelihood

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
    avg_ll.push_back(dgn_avg_ll);
    var_ll.push_back(dgn_var_ll());

    // FLAG
    std::cerr << "\n"
              << __rung << ")"
              << " it = " << inv_t[0]
              << " | avg_ll = " << avg_ll[0]
              << " | var_ll = " << var_ll[0]
              << "\n"
              << std::endl;
    // FLAG-END

    //=================================//
    //===== INVERSE TEMPERATURE 1 =====//
    //=================================//

    setup_rnd_init_params = true;
    __rung = 1;
    setup_inv_T = 1.;
    run_monte_carlo(progbar, __mcmc_sweeps, true);
    inv_t.push_back(setup_inv_T);
    avg_ll.push_back(dgn_avg_ll);
    var_ll.push_back(dgn_var_ll());

    // FLAG
    std::cerr << "\n"
              << __rung << ")"
              << " it = " << inv_t[1]
              << " | avg_ll = " << avg_ll[1]
              << " | var_ll = " << var_ll[1]
              << "\n"
              << std::endl;
    // FLAG-END

    //=================================//
    //===== INVERSE TEMPERATURE 1 =====//
    //=================================//

    setup_rnd_init_params = false;
    size_t k = 0;

    // trapezoidal rule to compute marginal likelihood with power posteriors
    while (__rung < __num_rungs)
    {
        __rung++;

        //----- adaptive choice of inverse temperature -----//
        setup_inv_T = (avg_ll[k + 1] - inv_t[k + 1] * var_ll[k + 1] - avg_ll[k] + inv_t[k] * var_ll[k]) / (var_ll[k] - var_ll[k + 1]);
        
        // FLAG
        std::cerr << "\n" << setup_inv_T
                  << std::endl;
        // FLAG-END

        if ((setup_inv_T <= inv_t[k]) || (setup_inv_T >= inv_t[k + 1]))
        {
            // FLAG
            std::cerr << "\n MEANY"
                      << std::endl;
            // FLAG-END
            setup_inv_T = inv_t[k] + var_ll[k + 1] * (inv_t[k + 1] - inv_t[k]) / (var_ll[k] + var_ll[k + 1]);
        }

        //----- power posterior sampling -----//
        run_monte_carlo(progbar, __rung * __mcmc_sweeps, false);

        //----- find larger rectangular discretization error -----//
        is_right = (inv_t[k + 1] - setup_inv_T) * (avg_ll[k + 1] - dgn_avg_ll) > (setup_inv_T - inv_t[k]) * (dgn_avg_ll - avg_ll[k]);

        // FLAG
        std::cerr << "\n"
                  << "is_right = " << is_right << " | "
                  << "(" << inv_t[k + 1] << " - " << setup_inv_T << ") * " << (avg_ll[k + 1] - dgn_avg_ll)
                  << " ?> "
                  << "(" << setup_inv_T << " - " << inv_t[k] << ") * " << (dgn_avg_ll - avg_ll[k])
                  << std::endl;
        // FLAG-END

        //----- insert in-between for all lists -----//
        inv_t.insert(inv_t.begin() + k + 1, setup_inv_T);
        avg_ll.insert(avg_ll.begin() + k + 1, dgn_avg_ll);
        var_ll.insert(var_ll.begin() + k + 1, dgn_var_ll());

        // FLAG
        std::cerr << "\n"
                  << "<< ";
        for (const auto &x : inv_t)
        {
            std::cerr << x << ", ";
        }
        std::cerr << " >>" << std::endl;
        // FLAG-END

        // move index to inserted if right interval is chosen
        k += (is_right);

        // FLAG
        if (is_right)
        {
            std::cerr << "\n"
                      << __rung << ")"
                      << " it = " << inv_t[k]
                      << " | avg_ll = " << avg_ll[k]
                      << " | var_ll = " << var_ll[k]
                      << " (RIGHT " << k << ")"
                      << "\n"
                      << std::endl;
        }
        else
        {
            std::cerr << "\n"
                      << __rung << ")"
                      << " it = " << inv_t[k + 1]
                      << " | avg_ll = " << avg_ll[k + 1]
                      << " | var_ll = " << var_ll[k + 1]
                      << " (LEFT " << k << ")"
                      << "\n"
                      << std::endl;
        }

        // FLAG-END
    }

    //----- trapezoidal rule to compute marginal likelihood -----//
    long double d_inv_t, marginal_loglike = 0.;
    for (int k = 0; k < num_rungs; k++)
    {
        d_inv_t = inv_t[k + 1] - inv_t[k];
        marginal_loglike += (avg_ll[k + 1] + avg_ll[k]) * d_inv_t / 2. - (var_ll[k + 1] - var_ll[k]) * std::pow(d_inv_t, 2) / 12.;
    }

    

    //----- Close progress bar -----//
    screen::progress_bar_end("end", progbar);

    //----- Restore user setup -----//
    setup_inv_T = bkp_inv_T;
    setup_rnd_init_params = bkp_rnd_init_params;

    return marginal_loglike;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Display results as sampled during iterations. */
void MCMC::show_results()
{
    using namespace screen;

    if (!__is_infer_complete)
    {
        raise_warning("Ignored attempt to show results before inference.",
                      Inconmplete_Infer_Diagnostics_Warning);
    }
    else
    {
        //===============================//
        //===== Compute Diagnostics =====//
        //===============================//

        // compute WAIC
        long double sample_correc = static_cast<long double>(dgn_sample) / (dgn_sample - 1);
        long double p_waic = 0., log_point_pred_dens = 0., effec_num_params_waic = 0.;

        //---- implement WAIC warning -----//
        // Ref: Vhetari et al. (2016) Stat Comput.

        Large_p_WAIC_Warning = false; // note :print warning at every call
        bool is_p_waic_large = false;

        const long double p_waic_limit = 0.4;
        unsigned p_waic_large_count = 0;

        // Loop over households
        for (int h = 0; h < num_househ; h++)
        {
            log_point_pred_dens += std::log(househ[h].dgn_likelihood);
            p_waic = (househ[h].dgn_sqr_ll - std::pow(househ[h].dgn_avg_ll, 2)) * sample_correc;
            if (p_waic > p_waic_limit)
            {
                is_p_waic_large = true;
                p_waic_large_count++;
            }
            effec_num_params_waic += p_waic;
        }

        // diagnostics (sampled every thinning iteration)
        std::ostringstream oss;
        oss << std::setprecision(2) << std::fixed
            << "loglikelihood_mean" << C_DELIM
            << dgn_avg_ll
            << "\nloglikelihood_sd" << C_DELIM
            << std::sqrt((dgn_sqr_ll - std::pow(dgn_avg_ll, 2)) * sample_correc)
            << "\nlogpoint_predictive_density" << C_DELIM
            << log_point_pred_dens
            << "\neffective_parameters" << C_DELIM
            << effec_num_params_waic
            << "\nWAIC" << C_DELIM
            << 2 * (effec_num_params_waic - log_point_pred_dens);
        print_stdout(oss.str() + "\n");

        // FIXME : this should be computed alone
        // practical WAIC warning
        if (is_p_waic_large == true)
        {
            std::ostringstream oss;
            oss << percent(p_waic_large_count, num_househ)
                << " p_waic estimates greater than " << p_waic_limit << "."
                << " Unreliable WAIC (Vehtari-Gelman-Gabry criterion).";
            raise_warning(oss.str(), Large_p_WAIC_Warning);
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

    //-----random initialization of augmented variables with model-----//
    init_augmented_variables();

    //----- Initialize model -----//
    init_model_for_simulation();

#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        simulate_single_house(h);
    }

    update_report_data();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : methods (display) ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* display the user model parameters */
void MCMC::show_parameters()
{
    using namespace screen;

    std::ostringstream oss;
    oss << param_file_headers[0];
    for (size_t k = 1; k < param_file_headers.size(); k++)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    for (std::string parID : param_to_be_inferred)
    {
        oss << "\n"
            << parID << C_DELIM << get_message_to_display(parID);
    }
    for (std::string parID : param_to_be_static)
    {
        oss << "\n"
            << parID << C_DELIM << get_message_to_display(parID);
    }

    print_stdout(oss.str() + "\n");
}

/* display the diagnostics */
void MCMC::show_param_diagnostics()
{
    using namespace screen;

    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";
    for (size_t k = 2; k < param_file_headers.size(); k++)
    {
        oss << C_DELIM << param_file_headers[k];
    }
    for (std::string parID : param_to_be_inferred)
    {
        oss << "\n"
            << parID << C_DELIM << get_message_to_diagnostics(parID);
    }

    print_stdout(oss.str() + "\n");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : methods (Rcpp) ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef USE_RCPP
Rcpp::List MCMC::return_setup()
{
    using namespace Rcpp;
    List output;

    output = List::create(
        _["seed"] = setup_seed,
        _["num_threads"] = get_num_threads(),
        _["num_chains"] = setup_num_chains,
        _["num_adapt"] = setup_num_adapt,
        _["burn_in"] = setup_burn_in,
        _["num_iter"] = setup_num_iter,
        _["thinning"] = setup_thinning,
        _["rnd_init_params"] = screen::truth(setup_rnd_init_params));

    if (setup_inv_T != 1.)
    {
        output["inv_temper"] = setup_inv_T;
    }

    return output;
}

#endif