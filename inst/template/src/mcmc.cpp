// mcmc.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Inconmplete_Infer_Diagnostics_Warning = false; // Warning for unreliable WAIC
bool Large_p_WAIC_Warning = false;                  // Warning for unreliable WAIC

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (user defined)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign group rate to an individual  */
void MCMC::assign_rate_per_group_to_individual(Individual &this_indiv)
{
    this_indiv.assign_rate_per_group("Infectivity", REFERENCE);
    this_indiv.assign_rate_per_group("Susceptibility", REFERENCE);
}

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_rate_per_contact_to_household(int h, int from, int to)
{
    househ[h].assign_rate_per_contact(from, to, REFERENCE);
}

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
double MCMC::compute_single_house_loglike(int h)
{
    UNUSED(h);
    return 0.;
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : methods (general)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute global loglikelihood */
double MCMC::compute_global_loglike()
{
    double output = 0.0;

#pragma omp parallel for reduction(+ : output)
    for (int h = 0; h < num_househ; h++)
    {
        double tmp = compute_single_house_loglike(h);
        househ[h].set_loglike(tmp);
        output += tmp;
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

    for (int i : househ[h].infectious_list)
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

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : methods (parameter move)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
void MCMC::single_param_move(std::string parID)
{
    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (propose_valid_param_move(parID))
    {
        // Compute the new gobla loglikelihood
        // set loglike in households and store old values in buffer
        double new_global_loglike = compute_global_loglike();

        //=====Metropolis-Hastings=====//
        double logQ = inverse_temperature * (new_global_loglike - get_global_loglike());
        logQ += get_move_logratio(); // logratio priors - logratio proposal

        if (log(uniform()) < logQ)
        {
            accept_param_move(parID);
            set_global_loglike(new_global_loglike);
        }
        else
        {
            refuse_param_move(parID);
            // restore loglike in households
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
        single_param_move(parID);
    }
}

/* Sweep over all model parameters plus adaptive step. */
void MCMC::move_and_adapt_all_parameters(ulong step)
{
    double adapt_ratio = static_cast<double>(step) / setup_num_adapt_iter;

    for (std::string parID : param_to_be_inferred)
    {
        single_param_move(parID);
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
double MCMC::single_household_augm_move(int h)
{
    double old_value, new_value, logQ;
    double old_house_loglike = 0.0, new_house_loglike = 0.0;
    double delta_global_loglike = 0.0;
    bool is_accepted = false;

    // Loop over infectious residents of the househould
    // FIXME : only infectious are subject to data augmentation?
    for (int i : househ[h].infectious_list)
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
                new_house_loglike = compute_single_house_loglike(h);

                // metropolis-hasting contribution
                logQ = inverse_temperature * (new_house_loglike - old_house_loglike);
                logQ -= househ[h].indiv[i].augm_propose_logratio(augmID, old_value, new_value);

                if (log(uniform()) < logQ)
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
    if (setup_data_augm)
    {
        double variation = 0.;

#pragma omp parallel for reduction(+ : variation)
        for (int h = 0; h < num_househ; h++)
        {
            variation += single_household_augm_move(h);
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
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : display likelihood  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method to only compute the loglikelihood and display it. */
double MCMC::run_compute_loglikelihood(bool verbose)
{
    //===== Initialization =====//

    //-----random initialization of augmented variables with model-----//
    init_augmented_variables();

    //----- Initialize model -----//
    initialize_model();

    //===== Execution =====//
    double ll = compute_global_loglike();

    //----- Print output -----//
    if (verbose)
    {
        std::ostringstream oss;
        oss << "Loglikelihood" << screen::C_DELIM << ll << "\n";
        screen::print_stdout(oss.str());
    }
#ifndef USE_RCPP
    else
    {
        screen::print_stdout(std::to_string(ll) + "\n");
    }
#endif

    return ll;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MCMC : inference  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings with power posteriors. */
void MCMC::run_inference(std::string output_file, bool verbose, bool progbar)
{
    std::string adap_prog_txt, burn_prog_txt, exec_prog_txt;

    //==========================//
    //===== Initialization =====//
    //==========================//

    screen::display_decorator("Inference", screen::C_DECOR, UNDEFINED, verbose);
    intern_is_infer_complete = false;

    //-----Initialize model-----//
    initialize_model();

    //-----Initialize output file-----//
    init_param_file(output_file);

    //-----Initialize diagnostics-----//
    init_model_diagnostics();

    //-----Multiple Chains-----//
    for (intern_chain = 1; intern_chain <= setup_num_chains; intern_chain++)
    {
        // Set progress bars text
        if (setup_num_chains == 1)
        {
            adap_prog_txt = "adaptive";
            burn_prog_txt = "burn-in";
            exec_prog_txt = "execution";
        }
        else 
        {
            adap_prog_txt = "adap " + std::to_string(intern_chain) + "/" + std::to_string(setup_num_chains);
            burn_prog_txt = "burn " + std::to_string(intern_chain) + "/" + std::to_string(setup_num_chains);
            exec_prog_txt = "exec " + std::to_string(intern_chain) + "/" + std::to_string(setup_num_chains);
        }

        //-----Initialize parameter inference-----//
        init_rand_params(setup_rnd_init_params);

        //-----Initialize augmented variables with model-----//
        init_augmented_variables();

        //-----Initialize chain-----//
        initialize_chain();

        // set the value of the global loglikelihood
        set_global_loglike(compute_global_loglike());

        // avoid asserting parameters definitions after first call
        disable_assert_param();

        // reset the acceptance counts
        init_move_counts();

        //=======================================//
        //=====Adaptative MCMC on parameters=====//
        //=======================================//

        if (setup_num_adapt_iter > 0)
        {
            screen::display_progress_bar_start(adap_prog_txt, progbar);
            for (ulong i = 1; i <= setup_num_adapt_iter; i++)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(i);
                screen::display_progress_bar(adap_prog_txt, i, setup_num_adapt_iter, progbar);
            }
            screen::display_progress_bar_end(adap_prog_txt, progbar);

            // reset the acceptance counts
            init_move_counts();
        }

        //============================================//
        //========= MCMC Burn-in Iterations ==========//
        //============================================//

        if (setup_burnin > 0)
        {
            screen::display_progress_bar_start(burn_prog_txt, progbar);
            for (ulong i = 1; i <= setup_burnin; i++)
            {
                move_all_augm_variables();
                move_all_parameters();
                screen::display_progress_bar(burn_prog_txt, i, setup_num_prints, progbar);
            }
            screen::display_progress_bar_end(burn_prog_txt, progbar);
        }

        //============================================//
        //========= MCMC Posterior Sampling ==========//
        //============================================//

        screen::display_progress_bar_start(exec_prog_txt, progbar);
        for (ulong i = 1; i <= setup_num_prints; i++)
        {
            for (ulong j = 1; j <= setup_thinning; j++)
            {
                move_all_augm_variables();
                move_all_parameters();
            }

            // compute diagnostics for the previous iterations
            update_model_diagnostics();

            print_parameters();
            screen::display_progress_bar(exec_prog_txt, i, setup_num_prints, progbar);
        }
        screen::display_progress_bar_end(exec_prog_txt, progbar);

        // FIXME : 
        // missing a storing diagnostics : everything is canceled for each chain
        // a Gelman Rubin could be integrated then
    }

    //==========================//
    //=====end of inference=====//
    //==========================//

    intern_is_infer_complete = true;
    close_param_file();

    if (verbose)
    {
        screen::display_decorator("Inference::results", '-', UNDEFINED);
        display_param_diagnostics();
        screen::display_decorator("Inference::diagnostics", '-', UNDEFINED);
        display_thin_results();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void MCMC::init_model_diagnostics()
{
    ulong sample_size = setup_num_chains * setup_num_prints; 

    // diagnostics for all parameters
    init_param_diagnostics(sample_size);

// diagnostics over households
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].init_diagnostics(sample_size);
    }

    dgn_loglike = 0;     // (global) loglikelihood
    dgn_sqr_loglike = 0; // mean squared (global) loglikelihood
}

/* Update Bayesian diagnostics */
void MCMC::update_model_diagnostics()
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
    dgn_loglike += global_loglike / (setup_num_chains * setup_num_prints);    
    
    // mean squared (global) loglikelihood
    dgn_sqr_loglike += std::pow(global_loglike, 2) / (setup_num_chains * setup_num_prints);
}

/* Display results as sampled during thin iterations.
// FIXME : this output needs to be user specified probably
*/
void MCMC::display_thin_results()
{
    using namespace screen; 

    if (!intern_is_infer_complete)
    {
        raise_warning("Ignored attempt to show results before inference.", Inconmplete_Infer_Diagnostics_Warning);
    }
    else
    {
        //===========================//
        //===== Compute Results =====//
        //===========================//

        // compute WAIC
        double sample_size = setup_num_chains * setup_num_prints;
        double samp_var_corr = (sample_size) / (sample_size - 1.);
        double p_waic = 0., log_point_pred_dens = 0., effec_num_params_waic = 0.;

        // implement WAIC warning
        // Ref: Vhetari et al. (2016) Stat Comput.
        Large_p_WAIC_Warning = false; // print warning at every call
        bool is_p_waic_large = false;
        int p_waic_large_count = 0;
        const double p_waic_limit = 0.4;

        // Loop over households
        for (int h = 0; h < num_househ; h++)
        {
            log_point_pred_dens += std::log(househ[h].dgn_likelihood);
            p_waic = (househ[h].dgn_sqr_loglike - std::pow(househ[h].dgn_loglike, 2)) * samp_var_corr;
            if (p_waic > p_waic_limit)
            {
                is_p_waic_large = true;
                p_waic_large_count++;
            }
            effec_num_params_waic += p_waic;
        }
        // diagnostics (sampled every thinning iterration)
        std::ostringstream oss;
        oss << "loglikelihood.mean" << C_DELIM << std::setprecision(2) << std::fixed << dgn_loglike
            << "\nloglikelihood.sd" << C_DELIM << (dgn_sqr_loglike - std::pow(dgn_loglike, 2)) * samp_var_corr
            << "\nlogpoint.predictive.density" << C_DELIM << log_point_pred_dens
            << "\neffective.parameters" << C_DELIM << effec_num_params_waic
            << "\nWAIC" << C_DELIM << 2 * (effec_num_params_waic - log_point_pred_dens);
        print_stdout(oss.str() + "\n");

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

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : methods (simulation)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the time of infection of a patient. */
void MCMC::simulate_single_house(int h)
{
    UNUSED(h);
}

/* Perform simulation as provided by the user */
void MCMC::run_simulation(std::string output_filename, bool verbose)
{
    if (output_filename == UNDEFINED)
    {
        screen::raise_error("Missing output file for simulation.");
    }
    else
    {
        //=====INITIALIZATION=====//
        screen::display_decorator("Simulation", screen::C_DECOR, UNDEFINED, verbose);

        init_augmented_variables();
        initialize_model();

#pragma omp parallel for default(shared)
        for (int h = 0; h < num_househ; h++)
        {
            simulate_single_house(h);
        }

        // print results as augmented variable file
        print_augm_variables(output_filename, verbose);
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : methods (display) ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* display the user model parameters */
void MCMC::display_parameters()
{
    using namespace screen;

    std::ostringstream oss;
    oss << param_file_required_headers[0];
    for (size_t k = 1; k < param_file_required_headers.size(); k++)
    {
        oss << C_DELIM << param_file_required_headers[k];
    }
    for (std::string parID : param_to_be_inferred)
    {
        oss << "\n"
            << parID << C_DELIM << get_message_to_display(parID);
    }
    for (std::string parID : param_static)
    {
        oss << "\n"
            << parID << C_DELIM << get_message_to_display(parID);
    }
    
    print_stdout(oss.str() + "\n");
}

/* display the diagnostics */
void MCMC::display_param_diagnostics()
{
    using namespace screen;

    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";
    for (size_t k = 2; k < param_file_required_headers.size(); k++)
    {
        oss << C_DELIM << param_file_required_headers[k];
    }
    for (std::string parID : param_to_be_inferred)
    {
        oss << "\n"
            << parID << C_DELIM << get_message_to_diagnostics(parID);
    }
    
    print_stdout(oss.str() + "\n");
}