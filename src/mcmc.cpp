// mcmc.cpp //

#include "mcmc.h"

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
void MCMC::assign_rate_per_contact(int i, int k)
{
    UNUSED(i);
    UNUSED(k);
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void MCMC::store_augm_var_dependecies(int i)
{
    UNUSED(i);
}

/* Update household first case and first infection time.
WARNING: this strucutre assumes a single first infector.
FIXME : this should be a function of the time-models maybe
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
    househ[h].set_data("start_followup", new_t0);
}

/* */ // FIXME : is then study_period a standard name? In general, how to load household data?
// FIXME : this should be a function of the time-models maybe
void MCMC::assign_end_followup(int h)
{
    // find maximum study period within the household
    double end_followup = 0;
    for (int i : househ[h].infectious_list)
    {
        end_followup = std::max(end_followup, indiv[i].get_data("study_period"));
    }
    househ[h].set_data("end_followup", end_followup);
}

double MCMC::compute_single_house_loglike(int h)
{
    UNUSED(h);
    return 0.;
};

/* Compute global loglikelihood */
double MCMC::compute_global_loglike()
{
    double output = 0.0;
    double tmp;

    for (int h = 0; h < num_househ; h++)
    {
        tmp = compute_single_house_loglike(h);
        househ[h].set_loglike(tmp);
        output += tmp;
    }

    return output;
}

/* Set accept rates to zero. */
void MCMC::reset_moves()
{
    for (std::string augmID : augm_to_infer_set)
    {
        m_num_propos_data[augmID] = 0;
        m_num_accept_data[augmID] = 0;
    }

    for (std::string parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////  MCMC : move parameters /////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter */
void MCMC::_single_param_move_annealing(std::string parID)
{
    double Q, new_global_loglike;

    // proceed only if the new value is consistent, otherwise automatically execute refuse_move(parID)
    if (propose_valid_move(parID, rng))
    {
        // Compute the new log likelihood
        new_global_loglike = compute_global_loglike();

        //=====Metropolis-Hastings=====//
        Q = inverse_temperature * (new_global_loglike - get_global_loglike());
        Q += get_move_log_ratio();
        if (log(uniform(rng)) < Q)
        {
            accept_param_move(parID);
            set_global_loglike(new_global_loglike);
        }
        else
        {
            refuse_param_move(parID);
            // restore loglike in households
            for (int h = 0; h < num_househ; h++)
            {
                househ[h].refuse();
            }
        }
    }
}

/* Attempt to move a single parameter at zero inverse temperature */
void MCMC::_single_param_move_zero_inv_temp_(std::string parID)
{
    // proceed only if the new value is consistent
    // otherwise automatically execute refuse_move(parID)
    if (propose_valid_move(parID, rng))
    {
        //=====Metropolis-Hastings (only with priors)=====//
        if (log(uniform(rng)) < get_move_log_ratio())
        {
            accept_param_move(parID);
            set_global_loglike(compute_global_loglike());
        }
        else
        {
            refuse_param_move(parID);
            // restore loglike in households
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
        _single_param_move_annealing(parID);
    }
}

/* Sweep over all model parameters plus adaptive step. */
void MCMC::move_and_adapt_all_parameters(int step, double target)
{
    double adapt_ratio = (double)step / setup_num_adapt_iter;

    for (std::string parID : param_to_be_inferred)
    {
        _single_param_move_annealing(parID);
        adapt(parID, adapt_ratio, target);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////  MCMC : move augmented variables //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/*
Sweep over augmented variables.
FIXME : data agumentation is not always time dependent
*/
void MCMC::move_all_augm_variables()
{
    if (augm_to_infer_set.size() > 0)
    {
        double old_val, new_val, Q;
        double log_ratio_proposal = 0.0;
        double old_loglike_of_HH = 0.0, new_loglike_of_HH = 0.0;

        // Loop over households
        for (int h = 0; h < num_househ; h++)
        {
            // Loop over infectious residents of the househould
            for (int i : househ[h].infectious_list)
            {
                for (std::string augmID : augm_to_infer_set)
                {
                    old_val = indiv[i].get_augm(augmID);          // store old value
                    old_loglike_of_HH = househ[h].get_loglike();  // store old loglike
                    new_val = indiv[i].augm_propose(augmID, rng); // propose new value

                    // proceed only if the new value is consistent
                    // UPDATE : always valid if prior proposer
                    if (indiv[i].augm_valid(augmID, new_val))
                    {                                       // Update with proposed value to perform computations
                        indiv[i].set_augm(augmID, new_val); // UPDATE : should be a method of Model
                        update_first_case(h);               // Update information of first cases
                        store_augm_var_dependecies(i);      // Update infectivity profiles related to idx

                        new_loglike_of_HH = compute_single_house_loglike(h);
                        log_ratio_proposal = indiv[i].augm_propose_log_ratio(augmID, old_val, new_val);

                        Q = inverse_temperature * (new_loglike_of_HH - old_loglike_of_HH);
                        Q -= log_ratio_proposal;

                        if (log(uniform(rng)) < Q)
                        {
                            househ[h].set_loglike(new_loglike_of_HH);
                            set_global_loglike(get_global_loglike() + new_loglike_of_HH - old_loglike_of_HH);
                            m_num_accept_data[augmID]++;
                        }
                        else
                        {
                            indiv[i].set_augm(augmID, old_val); // UPDATE: should be method Model
                            update_first_case(h);               // Update information of first cases
                            store_augm_var_dependecies(i);      // Update infectivity profiles related to idx
                        }
                    }
                    m_num_propos_data[augmID]++;
                }
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////  MCMC : update_diagnostics //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update Bayesian diagnostics */
void MCMC::update_diagnostics()
{
    // Loop over households
    for (int h = 0; h < num_househ; h++)
    {
        househ[h].update_dgn();
    }

    dgn_mean_ll = (global_loglike + dgn_mean_ll * dgn_sample_size) / (dgn_sample_size + 1);
    dgn_mean_sqr_ll = (std::pow(global_loglike, 2) + dgn_mean_sqr_ll * dgn_sample_size) / (dgn_sample_size + 1);
    dgn_sample_size++;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////  MCMC : inference /////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::inference(
    std::string output_file,
    bool is_rand_params,
    double inv_temper,
    bool return_diagnostics)
{
    // use the thin number of iteration (euclidean division)
    int num_iter_thin = setup_num_iter / setup_thinning;

    // FIXME : temporaru
    UNUSED(return_diagnostics);

    if (Verbose)
    {
        std::cout << "\n================MCMC==================" << std::endl;
        std::cout << "\n---------MCMC:initialization----------" << std::endl;
    }

    //=====Initialize output file=====//
    init_param_file(output_file);

    //=====Initialize parameter inference====//
    if (is_rand_params)
    {
        init_rand_params(rng);
    }
    initialize_model();
    set_global_loglike(compute_global_loglike());
    reset_moves();

    //=====Assign inverse temperature=====//
    if ((inv_temper < 0.) || (inv_temper > 1.))
    {
        std::string err_text = "Inverse temperature " + std::to_string(inv_temper) + " out of range.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
    else
    {
        inverse_temperature = inv_temper;
        if ((Verbose) && (inv_temper != 1.))
        {
            std::cout << "Inverse_temperature : " << inverse_temperature << std::endl;
        }
    }

    //=======================================//
    //=====Adaptative MCMC on parameters=====//
    //=======================================//

    // UPDATE : probably can be modularized

    // UPDATE: add option for zero inverse temperature

    if (Verbose && (setup_num_adapt_iter > 0))
    {
        std::cout << "\n------------MCMC:adaptive-------------" << std::endl;
    }

    for (int step = 1; step <= setup_num_adapt_iter; step++)
    {
        move_and_adapt_all_parameters(step);
        move_all_augm_variables();

        if (Verbose)
        {
            print_progress_bar(step, setup_num_adapt_iter);
        }
    }

    //============================================//
    //=====Monte Carlo Markov Chain execution=====//
    //============================================//

    if (Verbose)
    {
        std::cout << "\n-----------MCMC:execution-------------" << std::endl;
    }

    for (int i = 1; i <= num_iter_thin; i++)
    {
        for (int j = 0; j < setup_thinning; j++)
        {
            move_all_parameters();
            move_all_augm_variables();
        }
        print_parameters(i * setup_thinning);

        if (Verbose)
        {
            print_progress_bar(i, num_iter_thin);
        }
    }

    //==========================//
    //=====end of inference=====//
    //==========================//

    close_param_file();
    if (Verbose)
    {
        std::cout << "\n--------------MCMC:end----------------" << std::endl;
        display_parameters();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////////  MCMC : simulation /////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::simulation(std::string output_filename)
{
    if (Verbose)
    {
        std::cout << "\n=============SIMULATION===============" << std::endl;
    }

    initialize_model();
    simulate_infection();
    print_augm_variables(output_filename);
}
