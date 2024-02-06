// mcmc.cpp //

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////  MCMC : setter ////////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////


/* Initialize parameters and augmented variables */
void MCMC::initialize_model(bool is_rand_params)
{
    if (is_rand_params)
    {
        init_rand_params(rng);
    }

    assign_infec_profile_function(); // initialize infectivity profiles

    // compute infectivity profiles for each of the symptomatic and relative infectees
    for (size_t i = 0; i < num_indiv; i++)
    {
        store_augm_var_dependecies(i);
    }

    // find the first case for each household
    for (size_t h = 0; h < num_househ; h++)
    {
        househ[h].assign_end_followup();
        househ[h].update_first_case();
    }
}

/* Compute global loglikelihood */
// UPDATE : emberassing parallel
double MCMC::compute_global_loglike()
{
    double output = 0.0;
    double tmp;

    for (size_t h = 0; h < num_househ; h++)
    {
        tmp = compute_loglike(h);
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
void MCMC::move_single_param(std::string parID)
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
        if (log(unif(rng)) < Q)
        {
            accept_param_move(parID);
            set_global_loglike(new_global_loglike);
        }
        else
        {
            refuse_param_move(parID);
            // restore loglike in households
            for (size_t h = 0; h < num_househ; h++)
            {
                househ[h].refuse();
            }
        }
    }
}

/* Sweep over all model parameters. */
void MCMC::move_parameters()
{
    std::string parID = UNDEFINED;

    for (size_t z = 0; z < num_param_to_infer; z++)
    {
        parID = param_to_be_inferred[param_picker(rng)];
        move_single_param(parID);
    }
}

/* Sweep over all model parameters plus adaptive step. */
void MCMC::move_and_adapt_parameters(unsigned step, double target)
{
    std::string parID = UNDEFINED;
    double adapt_ratio = (double)(1.0 + step) / setup.num_adapt_iter;

    for (size_t z = 0; z < num_param_to_infer; z++)
    {
        parID = param_to_be_inferred[param_picker(rng)];
        move_single_param(parID);
        adapt(parID, adapt_ratio, target);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////  MCMC : move augmented variables //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Sweep over augmented variables. */
void MCMC::move_augm_variables()
{
    double old_val, new_val, Q;
    double log_ratio_proposal = 0.0;
    double old_loglike_of_HH = 0.0, new_loglike_of_HH = 0.0;

    // Loop over households
    for (size_t h = 0; h < num_househ; h++)
    {
        // Loop over infectious residents of the househould
        for (size_t i : househ[h].infectors_list)
        {
            for (std::string augmID : augm_to_infer_set)
            {
                old_val = indiv[i].get_augm(augmID);          // store old value
                old_loglike_of_HH = househ[h].get_loglike();  // store old loglike
                new_val = indiv[i].augm_propose(augmID, rng); // propose new value

                // FIXME : always valid if prior proposer
                // proceed only if the new value is consistent
                if (indiv[i].augm_valid(augmID, new_val))
                {                                       // Update with proposed value to perform computations
                    indiv[i].set_augm(augmID, new_val); // UPDATE: should be method Model
                    househ[h].update_first_case();      // Update information of first cases
                    store_augm_var_dependecies(i);      // Update infectivity profiles related to idx

                    new_loglike_of_HH = compute_loglike(h);
                    log_ratio_proposal = indiv[i].augm_propose_log_ratio(augmID, old_val, new_val);

                    Q = inverse_temperature * (new_loglike_of_HH - old_loglike_of_HH);
                    Q -= log_ratio_proposal;

                    if (log(unif(rng)) < Q)
                    {
                        househ[h].set_loglike(new_loglike_of_HH);
                        set_global_loglike(get_global_loglike() + new_loglike_of_HH - old_loglike_of_HH);
                        m_num_accept_data[augmID]++;
                    }
                    else
                    {
                        indiv[i].set_augm(augmID, old_val); // UPDATE: should be method Model
                        househ[h].update_first_case();      // Update information of first cases
                        store_augm_var_dependecies(i);      // Update infectivity profiles related to idx
                    }
                }
                m_num_propos_data[augmID]++;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////  MCMC : inference /////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::inference(std::string output_file, bool is_rand_params, double inv_temper)
{
    if (Verbose)
    {
        std::cout << "\n================MCMC==================" << std::endl;
        std::cout << "\n---------MCMC:initialization----------" << std::endl;
    }

    initialize_model(is_rand_params);
    set_global_loglike(compute_global_loglike());
    reset_moves();

    //=====Initialize output file=====//
    init_param_file(output_file);

    //=====Assign ivnverse temperature=====//
    // UPDATE : code can be optimized to avoid * 1 or * 0
    if (inv_temper == INFINITY)
    {
        inverse_temperature = 1.0;
    }
    else if ((inverse_temperature < 0) || (inverse_temperature > 1))
    {
        std::string err_text = "Inverse temperature " + std::to_string(inverse_temperature) + " out of range.";
        exit_with_failure_error(err_text);
    }
    else
    {
        inverse_temperature = inv_temper;
        if (Verbose)
        {
            std::cout << "Inverse_temperature : " << inverse_temperature << std::endl;
        }
    }

    //=====Adaptative MCMC on parameters=====//
    if (Verbose && (setup.num_adapt_iter > 0))
    {
        std::cout << "\n------------MCMC:adaptive-------------" << std::endl;
    }
    for (unsigned step = 0; step < setup.num_adapt_iter; step++)
    {
        move_and_adapt_parameters(step);
        move_augm_variables();
    }

    //=====MCMC execution=====//
    if (Verbose)
    {
        std::cout << "\n-----------MCMC:execution-------------" << std::endl;
    }
    for (unsigned i = 0; i < setup.num_iter; i += setup.thinning)
    {
        for (unsigned j = 0; j < setup.thinning; j++)
        {
            move_parameters();
            move_augm_variables();
        }
        print_parameters(i + setup.thinning);
    }

    //=====MCMC ending=====//
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

    initialize_model(false);
    simulate_infection(setup.delta_t, rng);
    print_augm_variables(output_filename);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////  MCMC : add individuals and context ///////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Add individual */
void MCMC::add_individual(Individual new_indiv)
{
    int status = new_indiv.get_status();
    std::string parID = UNDEFINED;

    //=====check for the validity of the status (FIXME: temporary)=====//
    if ((status != HEALTHY) && (status != SYMPTOMATIC))
    {
        std::string err_text = "Unkown status: " + std::to_string(status) + ".\nCheck your model.";
        exit_with_failure_error(err_text);
    }
    else
    {
        //-----Assign Model Rates-----//
        new_indiv._allocate_rates(group_rate_list.size());
        assign_rate_per_group_to_individual(&new_indiv);
        for (size_t r = 0; r < group_rate_list.size(); r++)
        {
            if (new_indiv._which_rate[r] == UNDEFINED)
            {
                std::string err_text = "Rate `" + group_rate_list[r] + "` is not provided for individual " + std::to_string(new_indiv.get_index()) + ".";
                exit_with_failure_error(err_text);
            }
            else if (new_indiv._which_rate[r] == REFERENCE)
            {
                new_indiv._set_rate(r, get_reference());
                model_report[r].at("[reference]") += 1;
            }
            else
            {
                parID = new_indiv._which_rate[r];
                new_indiv._set_rate(r, get_address(parID));
                if (model_report[r].find(parID) == model_report[r].end())
                {
                    model_report[r].emplace(parID, 1);
                }
                else
                {
                    model_report[r].at(parID) += 1;
                }
            }
        }

        //-----Store New Individual-----//
        indiv[num_indiv] = new_indiv;
        num_indiv++;
    }
}

/* Assign community structure to individuals. I returns a report. */
void MCMC::assign_community_structure(std::string HH_index_column)
{
    size_t h, size;
    std::unordered_map<int, Household> load_househ = {};

    //=====safe load households=====//
    for (size_t i = 0; i < num_indiv; i++)
    {
        h = indiv[i].get_data(HH_index_column);
        // if househ_idx doesn't exist, append new household
        if (load_househ.find(h) == load_househ.end())
        {
            load_househ.emplace(h, Household());
            load_househ.at(h).set_index(h);
        }
        load_househ.at(h).link_individual(i, indiv[i].get_status());
    }

    //=====produce list of households=====//
    num_househ = load_househ.size();    // get number of unique household index
    househ = new Household[num_househ]; // initialize the array of household
    h = 0;
    for (const auto &pair : load_househ)
    {
        househ[h] = pair.second;
        househ[h].set_indiv_ptr(&indiv); // associate the array of individuals to the Structure
        h++;
    }
    load_househ.clear(); // clear the dictionary for the loading

    //=====assign home to individuals=====//
    for (size_t h = 0; h < num_househ; h++)
    {
        size = 0;
        for (size_t i : househ[h].infectors_list)
        {
            indiv[i].home = h;
            size++;
        }
        for (size_t i : househ[h].susceptibles_list)
        {
            indiv[i].home = h;
            size++;
        }
        househ[h].size = size;
    }
}