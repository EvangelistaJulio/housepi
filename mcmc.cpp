// mcmc.cpp //

#include "mcmc.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : Setter ///////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables */
void MCMC::initialize_model(bool random_params)
{
    if (random_params)
    {
        init_rand_params(rng);
    }
    initialize_infectivity_profile(); 
    assign_augm_data(rng);
}

/* Initialize global loglikelihood */
// UPDATE : this initialization is used because changes are not stored for params update
void MCMC::init_global_loglike()
{
    double tmp;
    // #pragma omp parallel num_threads(NUM_THREAD)
    {
        global_loglike = 0.0;
        // #pragma omp for reduction(+ : sum)
        for (size_t h = 0; h < num_househ; h++)
        {
            tmp = compute_loglike(h);
            househ[h].set_loglike(tmp);
            global_loglike += tmp;
        }
    }
}

/* Set accept rates to zero */
void MCMC::reset_moves()
{
    m_num_propos_data = 0;
    m_num_accept_data = 0;
    for (std::string parID : param_to_be_inferred)
    {
        reset_accept_rate(parID);
    }
}

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : moves ////////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Sweep over model parameters. */
void MCMC::move_parameters(double inv_temper)
{
    double Q, new_global_loglike;

    // loop over all parameters to be inferred
    for (std::string parID : param_to_be_inferred)
    {
        // proceed only if the new value is consistent
        if (propose_valid_move(parID, rng))
        {
            // Compute the new log likelihood
            // UPDATE : emberassing parallel // FIXME : buffer
            new_global_loglike = 0.0;
            for (size_t h = 0; h < num_househ; h++)
            {
                new_global_loglike += compute_loglike(h);
            }

            //=====Metropolis-Hastings=====//
            Q = inv_temper * (new_global_loglike - global_loglike);
            Q += get_move_log_ratio();
            if (std::log(unif(rng)) < Q)
            {
                accept_move(parID);
                global_loglike = new_global_loglike;
            }
            else
            {
                refuse_move(parID);
            }
        }
    }
}

/* Sweep over augmented variables. */
void MCMC::move_augm_variables(double inv_temper)
{
    double old_val, new_val, Q;
    double log_ratio_proposal = 0.0;
    double old_loglike_of_HH = 0.0, new_loglike_of_HH = 0.0;

    // #pragma omp parallel num_threads(NUM_THREAD)
    {
        // #pragma omp for reduction(+ : sum)
        // Loop over households
        for (size_t h = 0; h < num_househ; h++)
        {
            // Loop over infectious residents of the househould
            for (size_t res : househ[h].symptomatic_list)
            {
                for (std::string augmID : augm_to_be_inferred)
                {
                    old_val = indiv[res].get_augm(augmID);          // store old value
                    old_loglike_of_HH = househ[h].get_loglike();    // store old loglike
                    new_val = indiv[res].augm_propose(augmID, rng); // propose new value

                    // FIXME : always valid if prior proposer
                    // proceed only if the new value is consistent
                    if (indiv[res].augm_valid(augmID, new_val))
                    {                                         // Update with proposed value to perform computations
                        indiv[res].set_augm(augmID, new_val); // UPDATE: should be method Model
                        update_first_case(h);                 // Update information of first cases
                        update_infectivity_profiles(res);     // Update infectivity profiles related to idx

                        new_loglike_of_HH = compute_loglike(h);
                        log_ratio_proposal = indiv[res].augm_propose_log_ratio(augmID, old_val, new_val);

                        Q = inv_temper * (new_loglike_of_HH - old_loglike_of_HH);
                        Q -= log_ratio_proposal;

                        if (std::log(unif(rng)) < Q)
                        {
                            global_loglike += new_loglike_of_HH - old_loglike_of_HH;
                            househ[h].set_loglike(new_loglike_of_HH);
                            m_num_accept_data++; // NOTE : problematic for normalization
                        }
                        else
                        {
                            indiv[res].set_augm(augmID, old_val); // UPDATE: should be method Model
                            update_first_case(h);                 // Update information of first cases
                            update_infectivity_profiles(res);     // Update infectivity profiles related to idx
                        }
                    }
                    m_num_propos_data++;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////
// /////////////////////////// MCMC : run /////////////////////////// //
////////////////////////////////////////////////////////////////////////

void MCMC::run(std::string output_file, double inv_temper, bool random_params)
{
    unsigned step = 0;
    double adapt_ratio;

    if (setup.verbose)
    {
        std::cout << "\n=============MCMC===============" << std::endl;
        std::cout << "\n======MCMC:initialization=======" << std::endl;
    }

    initialize_model(random_params);
    init_global_loglike();
    reset_moves();

    //=====Initialize output file=====//
    init_param_file(output_file);

    //=====MCMC execution=====//
    if (setup.verbose)
    {
        std::cout << "\n=========MCMC:execution=========" << std::endl;
    }
    for (unsigned i = 0; i < setup.num_iter; i += setup.thinning)
    {
        //=====MC sweeps over alla parameters to be inferred=====//
        for (unsigned j = 0; j < setup.thinning; j++)
        {
            //=====Change all parameters=====//
            step = 1 + j + i;
            move_parameters(inv_temper);
            init_global_loglike(); // FIXME: I reassign loglike to each household

            //=====Data augmentation=====//
            move_augm_variables(inv_temper);

            //=====Adaptative MCMC on parameters=====//
            if (step < setup.num_adapt_iter)
            {
                adapt_ratio = (double)step / setup.num_adapt_iter;
                for (std::string parID : param_to_be_inferred)
                {
                    adapt(parID, adapt_ratio);
                }
            }
        }
        print_parameters(step);
    }

    //=====MCMC ending=====//
    close_param_file();
}

////////////////////////////////////////////////////////////////////////
// /////////////////////////// MCMC : run /////////////////////////// //
////////////////////////////////////////////////////////////////////////

void MCMC::simulate(std::string output_filename)
{
    if (setup.verbose)
    {
        std::cout << "\n=============Simulation===============" << std::endl;
    }

    initialize_model(false);
    open_augm_file(output_filename);
    simulate_infection(setup.delta_t, rng);
    print_augm_variables();
    close_augm_file();
}

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : print ////////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Initialize the output param file with column names and initial values */
void MCMC::init_param_file(std::string output_filename)
{
    std::cout << "Model parameter file: " << output_filename << std::endl;
    param_ofstream.open(output_filename.c_str());

    if (param_ofstream.is_open())
    {
        param_ofstream << "iter" << setup.delim << "loglike";
        param_ofstream << setup.delim << "augm__a%";
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << setup.delim << parID << setup.delim << parID << "__a%";
        }
        param_ofstream << std::endl;
        print_parameters(0);
    }
    else
    {
        std::cerr << "\nERROR: Unable to open file `" << output_filename << "`" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }
}

/* Print the model parameters at the given iteration time */
void MCMC::print_parameters(unsigned step)
{
    if (param_ofstream.is_open())
    {
        param_ofstream << step << setup.delim << get_global_loglike();
        param_ofstream << setup.delim << 100.0 * ((double)m_num_accept_data / m_num_propos_data);
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << setup.delim << get_value(parID) << setup.delim << 100.0 * get_accept_rate(parID);
        }
        param_ofstream << std::endl;
    }
    else
    {
        std::cerr << "\nERROR: Output file closed unexpectedly." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }
}

/* Open the file for the augmented variables */
void MCMC::open_augm_file(std::string output_filename)
{
    augm_ofstream.open(output_filename.c_str());
    if (augm_ofstream.is_open())
    {
        std::cout << "Augmented variable file: " << output_filename << std::endl;

        //=====Print column headers=====//
        augm_ofstream << "index" << setup.delim << "status";
        for (std::string dataID : data_headers)
        {
            augm_ofstream << setup.delim << dataID;
        }
        for (std::string augmID : augm_list)
        {
            augm_ofstream << setup.delim << augmID;
        }
        augm_ofstream << std::endl;
    }
    else
    {
        std::cerr << "\nERROR: Unable to open file `" << output_filename << "`" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }
}

/* Print the ensemble of augmented variables */
void MCMC::print_augm_variables()
{
    for (size_t i = 0; i < num_indiv; i++)
    {
        augm_ofstream << indiv[i].get_index() << setup.delim << indiv[i].get_status();
        for (std::string dataID : data_headers)
        {
            augm_ofstream << setup.delim << indiv[i].get(dataID);
        }
        for (std::string augmID : augm_list)
        {
            augm_ofstream << setup.delim << indiv[i].get_augm(augmID);
        }
        augm_ofstream << std::endl;
    }
}