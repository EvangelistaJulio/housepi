// mcmc.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

/* Set the number of threads for parrell computations. */
void MCMC::set_num_threads(int n)
{
    omp_set_num_threads(n);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////  MCMC : ???? ////////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the time of infection of a patient. */
double MCMC::get_infec_time(int i) const
{
    UNUSED(i);
    return INFINITY;
}

/* Assign the contact patterns between a couple of individuals. */
std::string MCMC::assign_rate_per_contact(int i, int k)
{
    UNUSED(i);
    UNUSED(k);
    return REFERENCE;
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void MCMC::store_augm_var_dependecies(int i)
{
    UNUSED(i);
}

/* Update household first case and first infection time.
Note: this structure assumes a single first infector.
*/
void MCMC::update_first_case(int h)
{
    double t_infec;
    double new_t0 = INFINITY; // NOTE: init to impossible values
    int new_case0 = -1;       // NOTE: init to impossible values

    for (int i : househ[h].infectious_list)
    {
        t_infec = get_infec_time(i);
        new_t0 = std::min(t_infec, new_t0);
        if (t_infec == new_t0)
        {
            new_case0 = i;
        }
    }
    househ[h].first_case = new_case0;
}

/* Compute the loglikelihood for a single household. */
double MCMC::compute_single_house_loglike(int h)
{
    UNUSED(h);
    return 0.;
};

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

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////  MCMC : move parameters /////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set accept rates to zero. */
void MCMC::reset_moves()
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

/* Attempt to move a single parameter for the given inverse temperature. */
void MCMC::single_param_move(std::string parID)
{
    double logQ = -INFINITY;
    double new_global_loglike = 0.;

    // FLAG
    /*
    std::ostringstream oss;
    oss << parID << ": " << get_value(parID);
    bool is_valid = propose_valid_param_move(parID);
    oss << " -> " << get_value(parID)
        << ((is_valid) ? " (valid)" : " (not valid)");
    */
   // END-of-FLAG

    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (propose_valid_param_move(parID))
    //if (is_valid)
    {
        // Compute the new gobla loglikelihood
        // set loglike in households and store old values in buffer
        new_global_loglike = compute_global_loglike();

        //=====Metropolis-Hastings=====//
        logQ = inverse_temperature * (new_global_loglike - get_global_loglike());
        logQ += get_move_logratio(); // logratio priors - logratio proposal

        if (log(uniform()) < logQ)
        {
            accept_param_move(parID);
            set_global_loglike(new_global_loglike);
            
            // FLAG
            //oss << " | accept (" << exp(logQ) << ")";
            // END-of-FLAG
        }
        else
        {
            refuse_param_move(parID);
            // restore loglike in households
            // FLAG
#pragma omp parallel for default(shared)
            for (int h = 0; h < num_househ; h++)
            {
                househ[h].refuse();
            }

            // FLAG
            //oss << " | refuse (" << exp(logQ) << ")";
            // END-of-FLAG
        }
    }

    // FLAG
    //display_stdout(oss.str());
    // END-of-FLAG
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
void MCMC::move_and_adapt_all_parameters(int step, double target)
{
    double adapt_ratio = (double)step / setup_num_adapt_iter;

    for (std::string parID : param_to_be_inferred)
    {
        single_param_move(parID);
        adapt(parID, adapt_ratio, target);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////  MCMC : move augmented variables //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::update_augm_var(int i, std::string augmID, double new_val)
{
    indiv[i].set_augm(augmID, new_val);
    update_first_case(indiv[i].home); // Update information of first cases
    store_augm_var_dependecies(i);    // Update infectivity profiles related to idx
}

/* Move augmented variables for infectious members of a single household.
The move is performed one member at a time, one augmented variable at at time.
*/
double MCMC::single_household_augm_move(int h)
{
    double old_value, new_value, logQ;
    double old_house_loglike = 0.0, new_house_loglike = 0.0;
    double augm_model_logratio = 0., delta_global_loglike = 0;
    bool is_accepted = false;

    // Loop over infectious residents of the househould
    // FIXME : only infectious are subject to data augmentation?
    for (int i : househ[h].infectious_list)
    {
        // Loop over augmented variables of infectious
        for (std::string augmID : augm_to_be_inferred)
        {
            is_accepted = false;
            old_value = indiv[i].get_augm(augmID);           // store old value
            old_house_loglike = househ[h].get_loglike();     // store old loglike
            new_value = indiv[i].augm_propose(augmID, prng); // propose new value

            // proceed only if the new value is consistent
            if (indiv[i].augm_valid(augmID, new_value))
            {
                update_augm_var(i, augmID, new_value); // update value to perform computations
                new_house_loglike = compute_single_house_loglike(h);
                augm_model_logratio = indiv[i].augm_model_logratio(augmID, old_value, new_value);
                // metropolis-hasting contribution
                logQ = inverse_temperature * (augm_model_logratio + new_house_loglike - old_house_loglike);
                logQ -= indiv[i].augm_propose_logratio(augmID, old_value, new_value);

                if (log(uniform()) < logQ)
                {
                    is_accepted = true;
                    househ[h].set_loglike(new_house_loglike);
                    delta_global_loglike += new_house_loglike - old_house_loglike;
                }
                else
                {
                    update_augm_var(i, augmID, old_value); // restore old value
                }
            }
            househ[h].update_augm_accept_rate(augmID, is_accepted);

            // FLAG
            /*
            std::ostringstream oss;
            oss << i << ", " << augmID << ") " << (is_accepted? ": true" : ": false")
                << " | " << old_val << "
                (" << old_house_loglike << "), "
                << new_val << " (" << new_house_loglike << ")"
                << " | " << indiv[i].augm_model_logratio(augmID, old_val, new_val)
                << " | " << indiv[i].augm_propose_logratio(augmID, old_val, new_val)
                << " | " << logQ;
            display_stdout(oss.str());
            */
            // END-of-FLAG
        }
    }

    return delta_global_loglike;
}

/* Sweep over all households to move augmented variables. */
void MCMC::move_all_augm_variables()
{
    if (is_data_augmentation)
    {
        double output = 0.;

#pragma omp parallel for reduction(+ : output)
        for (int h = 0; h < num_househ; h++)
        {
            output += single_household_augm_move(h);
        }
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////  MCMC : diagnostics /////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update Bayesian diagnostics */
void MCMC::update_model_diagnostics()
{
    // diagnostics for all parameters
    update_diagnostics();

// diagnostics over households
#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].update_diagnostics();
    }

    // update the mean of the global loglikelihood
    dgn_sum_loglike += global_loglike;
    // update the mean of the squared global loglikelihood
    dgn_sum_sqr_loglike += std::pow(global_loglike, 2);
    // update the diagnostics sample size
    dgn_ss++;
}

/* Display Bayesian diagnostics */
void MCMC::display_model_diagnostics()
{
    double log_point_pred_dens = 0.;
    double effec_num_params_waic = 0.;

    // Loop over households
#pragma omp parallel for reduction(+ : log_point_pred_dens, effec_num_params_waic)
    for (int h = 0; h < num_househ; h++)
    {
        log_point_pred_dens += std::log(househ[h].dgn_sum_likelihood / househ[h].dgn_ss);
        effec_num_params_waic += est_sample_variance(househ[h].dgn_sum_loglike, househ[h].dgn_sum_sqr_loglike, househ[h].dgn_ss);
    }

    // display diagnostics of parameters
    display_diagnostics();

    std::cout << "loglikelihood_mean" << DELIM << std::setprecision(2) << std::fixed << dgn_sum_loglike / dgn_ss << "\n"
              << "loglikelihood_sd" << DELIM << est_sample_variance(dgn_sum_loglike, dgn_sum_sqr_loglike, dgn_ss) << "\n"
              << "logpoint_predict_dens" << DELIM << log_point_pred_dens << "\n"
              << "effective_parameters" << DELIM << effec_num_params_waic << "\n"
              << "WAIC" << DELIM << log_point_pred_dens - effec_num_params_waic;
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////  MCMC : inference /////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings with power posteriors. */
void MCMC::inference(std::string output_file, bool is_rand_params, double inv_temper)
{
    // use the thin number of iteration (euclidean division)
    int num_adapt_iter_thin = setup_num_adapt_iter / setup_thinning;
    int num_iter_thin = setup_num_iter / setup_thinning;

    display_decorator("INFERENCE");
    display_decorator("MCMC::initialization", '-', 1);

    //=====Initialize output file=====//
    init_param_file(output_file);

    //=====Initialize parameter inference====//
    if (is_rand_params)
    {
        init_rand_params();
    }
    initialize_model();
    set_global_loglike(compute_global_loglike());
    reset_moves();

    //=====Assign inverse temperature=====//
    // UPDATE : is someone interested in the case inv_temp > 1?
    if ((inv_temper < 0.) || (inv_temper > 1.))
    {
        std::string err_text = "Inverse temperature " + std::to_string(inv_temper) + " out of range.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
    else
    {
        inverse_temperature = inv_temper;
        if (inv_temper != 1.)
        {
            std::string text = "Inverse temperature " + std::to_string(inv_temper);
            display_stdout(text);
        }
    }

    //=======================================//
    //=====Adaptative MCMC on parameters=====//
    //=======================================//

    if (setup_num_adapt_iter > 0)
    {
        display_decorator("MCMC::adaptive", '-', 1);
        display_progress_bar("adaptive", 0, 100);
        for (int i = 1; i <= num_adapt_iter_thin; i++)
        {
            for (int j = 0; j < setup_thinning; j++)
            {
                move_all_augm_variables();
                move_and_adapt_all_parameters(j + i * setup_thinning);
            }
#ifdef PRINT_ADAPT_ITER
            print_parameters(i * setup_thinning);
#endif
            display_progress_bar("adaptive", i, num_adapt_iter_thin);
        }
        display_progress_bar("adaptive", 100, 100);
    }

    //============================================//
    //=====Monte Carlo Markov Chain execution=====//
    //============================================//

    if (num_iter_thin > 0)
    {
        display_decorator("MCMC::execution", '-', 1 + (setup_num_adapt_iter > 0));
        display_progress_bar("execution", 0, 100);
        for (int i = 1; i <= num_iter_thin; i++)
        {
            for (int j = 0; j < setup_thinning; j++)
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
            display_progress_bar("execution", i, num_iter_thin);
        }
        display_progress_bar("execution", 100, 100);
    }
    //==========================//
    //=====end of inference=====//
    //==========================//

    close_param_file();
    display_decorator("MCMC::end", '-', 2);
    if (num_iter_thin > 1)
    {
        // FIXME : this output needs to be user specified probably
        display_model_diagnostics();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////  MCMC : simulation /////////////////////////////////////// //
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
    initialize_model();

#pragma omp parallel for default(shared)
    for (int h = 0; h < num_househ; h++)
    {
        simulate_single_house(h);
    }
    format_simulation();
}
