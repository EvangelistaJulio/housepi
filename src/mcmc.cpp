// mcmc.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

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
void MCMC::move_and_adapt_all_parameters(int step)
{
    double adapt_ratio = (double)step / setup_num_adapt_iter;

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
    if (is_data_augmentation)
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
/////////////////////////////////  MCMC : methods (inference)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings with power posteriors. */
void MCMC::inference(std::string output_file, bool is_rand_params, double inv_temper)
{
    //=====Initialization=====//
    display_decorator("INFERENCE");
    display_decorator("MCMC::initialization", 1, '-');

    //-----Initialize output file-----//
    init_param_file(output_file);

    //-----Initialize parameter inference-----//
    init_rand_params(is_rand_params);

    //-----random initialization of augmented variables with model-----//
    init_augmented_variables();

    //-----Initialize model-----//
    initialize_model();
    set_global_loglike(compute_global_loglike());
    disable_assert_param(); // avoid asserting definitions of parameters

    //-----Assign inverse temperature-----//
    // UPDATE : is someone interested in the case inv_temp > 1?
    if ((inv_temper < 0.) || (inv_temper > 1.))
    {
        display_error_message("Inverse temperature " + std::to_string(inv_temper) + " out of range.");
        std::exit(EXIT_FAILURE);
    }
    else
    {
        inverse_temperature = inv_temper;
        if (inv_temper != 1.)
        {
            display_stdout("Inverse temperature " + std::to_string(inv_temper));
        }
    }

    //=======================================//
    //=====Adaptative MCMC on parameters=====//
    //=======================================//

    if (setup_num_adapt_iter > 0)
    {
        int num_adapt_iter_thin = setup_num_adapt_iter / setup_thinning;

        display_decorator("MCMC::adaptive", 1, '-');
        display_init_progress_bar("adaptive");

        init_move_counts();
        for (int i = 1; i <= num_adapt_iter_thin; i++)
        {
            for (int j = 1; j <= setup_thinning; j++)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(j + (i - 1) * setup_thinning);
            }
#ifdef PRINT_ADAPT_ITER
            print_parameters(i * setup_thinning);
#endif
            display_progress_bar("adaptive", i, num_adapt_iter_thin);
        }
        display_close_progress_bar("adaptive");
    }

    //============================================//
    //=====Monte Carlo Markov Chain execution=====//
    //============================================//

    if (setup_num_prints > 0)
    {
        display_decorator("MCMC::execution", 1 + (setup_num_adapt_iter > 0), '-');
        display_init_progress_bar("execution");

        init_move_counts();
        init_model_diagnostics();
        for (int i = 1; i <= setup_num_prints; i++)
        {
            for (int j = 1; j <= setup_thinning; j++)
            {
                move_all_augm_variables();
                move_all_parameters();
            }

            // compute diagnostics for the previous iterations
            update_model_diagnostics();

#ifdef PRINT_ADAPT_ITER
            print_parameters(setup_num_adapt_iter + i * setup_thinning);
#else
            print_parameters(i * setup_thinning);
#endif
            display_progress_bar("execution", i, setup_num_prints);
        }
        display_close_progress_bar("execution");
    }

    //==========================//
    //=====end of inference=====//
    //==========================//

    close_param_file();
    display_thin_results();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (diagnostics) ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init Bayesian diagnostics. */
void MCMC::init_model_diagnostics()
{
    // diagnostics for all parameters
    init_param_diagnostics(setup_num_prints);

// diagnostics over households
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].init_diagnostics(setup_num_prints);
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

    dgn_loglike += global_loglike / setup_num_prints;                  // (global) loglikelihood
    dgn_sqr_loglike += std::pow(global_loglike, 2) / setup_num_prints; // mean squared (global) loglikelihood
}

/* Display results as sampled during thin iterations.
// FIXME : this output needs to be user specified probably
*/
void MCMC::display_thin_results()
{
    if (setup_num_iter >= setup_thinning)
    {
        double sample_var_correc = setup_num_prints / (setup_num_prints - 1.);
        double log_point_pred_dens = 0.;
        double effec_num_params_waic = 0.;
        std::ostringstream oss;

        // Loop over households
#pragma omp parallel for reduction(+ : log_point_pred_dens, effec_num_params_waic)
        for (int h = 0; h < num_househ; h++)
        {
            log_point_pred_dens += std::log(househ[h].dgn_likelihood);
            effec_num_params_waic += (househ[h].dgn_sqr_loglike - std::pow(househ[h].dgn_loglike, 2)) * sample_var_correc;
        }

        // parameters (sampled every thinning iterration)
        display_decorator("MCMC::results", 2, '-');
        display_param_diagnostics();

        // diagnostics (sampled every thinning iterration)
        display_decorator("MCMC::diagnostics", 1, '-');
        oss << "loglikelihood_mean" << C_DELIM << std::setprecision(2) << std::fixed << dgn_loglike
            << "\nloglikelihood_sd" << C_DELIM << (dgn_sqr_loglike - std::pow(dgn_loglike, 2)) * sample_var_correc
            << "\nlogpoint_predict_dens" << C_DELIM << log_point_pred_dens
            << "\neffective_parameters" << C_DELIM << effec_num_params_waic
            << "\nWAIC" << C_DELIM << log_point_pred_dens - effec_num_params_waic;
        display_stdout(oss);
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
void MCMC::simulation()
{
    display_decorator("SIMULATION");

    //=====INITIALIZATION=====//
    init_augmented_variables();
    initialize_model();

#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        simulate_single_house(h);
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//