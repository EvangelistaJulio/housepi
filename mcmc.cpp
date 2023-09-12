// mcmc.cpp //

#include "mcmc.hpp"

using namespace std;
using namespace std::chrono;

/////////////////////////////////////////////////////
//  MCMC_Class : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Methods : Setter
// ----------------------------

// Modification of attributes
void MCMC_Class::reset_moves()
{
    m_num_propos_data = 0;
    m_num_accept_data = 0;
    for (std::string parID : glob_params.to_be_inferred)
    {
        glob_params.reset_accept_rate(parID);
    }
}

/* Initialize augmented variables */
void MCMC_Class::init_augm_data(std::mt19937_64 &rng)
{
    for (size_t h = 0; h < num_households(); h++)
    {
        m_data[h].assign_infectivity_profiles(&infectivity_profile);
        m_data[h].init_augm_data(rng);
    }
}

/* Initialize global loglikelihood */
// UPDATE : this initialization is used because changes are not stored for params update
void MCMC_Class::update_global_loglike()
{
    // #pragma omp parallel num_threads(NUM_THREAD)
    {
        m_global_loglike = 0.0;
        // #pragma omp for reduction(+ : sum)
        for (size_t h = 0; h < num_households(); h++)
        {
            m_data[h].h_loglike = m_data[h].compute_loglike();
            m_global_loglike += m_data[h].h_loglike;
        }
    }
}

// ----------------------------
//  Methods : Computations
// ----------------------------

/* Sweep over model parameters. */
void MCMC_Class::update_parameters()
{
    double old_value, new_value;
    double Q, log_ratio_proposal, log_ratio_prior, new_global_loglike;
    bool move_accepted;

    // loop over all parameters to be inferred
    for (std::string parID : glob_params.to_be_inferred)
    {
        old_value = glob_params.get_value(parID);      // store old value
        new_value = glob_params.propose(parID, m_rng); // propose new value

        // proceed only if the new value is consistent
        if (glob_params.invalid_param(parID, new_value))
        {
            move_accepted = false;
        }
        else
        {
            log_ratio_proposal = glob_params.log_ratio_proposal(parID, old_value, new_value);
            log_ratio_prior = glob_params.log_ratio_prior(parID, old_value, new_value);

            // Compute the new log likelihood
            glob_params.set_value(parID, new_value);

            // FIXME : emberassing parallel
            new_global_loglike = 0.0;
            for (size_t h = 0; h < num_households(); h++)
            {
                new_global_loglike += m_data[h].compute_loglike();
            }

            // Metropolis-Hastings
            Q = new_global_loglike - m_global_loglike + log_ratio_prior - log_ratio_proposal;
            if (std::log(m_unif(m_rng)) < Q)
            {
                move_accepted = true;
                m_global_loglike = new_global_loglike; // Update log likelihood
            }
            else
            {
                move_accepted = false;
                glob_params.set_value(parID, old_value);
            }
        }
        glob_params.update_accept_rate(parID, move_accepted);
    }
}

/* Sweep over augmented variables. */
void MCMC_Class::update_augm_variables()
{
    bool move_accepted = false;
    double old_value, new_value, Q;
    double log_ratio_proposal = 0.0;
    double old_loglike_of_HH = 0.0, new_loglike_of_HH = 0.0;

    // #pragma omp parallel num_threads(NUM_THREAD)
    {
        // #pragma omp for reduction(+ : sum)
        // Loop over households
        for (size_t h = 0; h < m_data.size(); h++)
        {
            // Loop over infectious residents of the househould
            for (size_t res : m_data[h].h_symptomatic_list)
            {
                for (std::string augmID : m_data[h].get_augm_to_be_inferred(res))
                {
                    old_value = m_data[h].get_augm(augmID, res);            // store old value
                    old_loglike_of_HH = m_data[h].h_loglike;                // store old loglike
                    new_value = m_data[h].augm_propose(augmID, res, m_rng); // propose new value

                    // FIXME (FUTURE) : is this improvable? this is basically always true
                    // FIXME : always valid if prior proposer
                    // proceed only if the new value is consistent
                    if (m_data[h].augm_invalid(augmID, res, new_value))
                    {
                        move_accepted = false;
                    }
                    else
                    { // Update with proposed value to perform computations
                        m_data[h].set_augm(augmID, res, new_value);
                        new_loglike_of_HH = m_data[h].compute_loglike();
                        log_ratio_proposal = m_data[h].augm_propose_log_ratio(augmID, res, old_value, new_value);
                        Q = new_loglike_of_HH - old_loglike_of_HH - log_ratio_proposal;

                        if (std::log(m_unif(m_rng)) < Q)
                        {
                            move_accepted = true;
                            m_global_loglike += new_loglike_of_HH - old_loglike_of_HH;
                            m_data[h].h_loglike = new_loglike_of_HH;
                            m_num_accept_data++; // NOTE : problematic for normalization
                        }
                        else
                        {
                            move_accepted = false;
                            m_data[h].set_augm(augmID, res, old_value);
                        }
                    }
                    m_num_propos_data++;
                }
            }
        }
    }
}

//////////////////////////////////////////////////
// /////////////// MCMC : infer /////////////// //
//////////////////////////////////////////////////

void MCMC_Class::infer(std::string output_file)
{
    unsigned this_iter;
    double adapt_ratio;

    if (m_setup.verbose)
    {
        std::cout << "\n=============MCMC===============" << std::endl;
        std::cout << "\n======MCMC:initialization=======" << std::endl;
    }

    glob_params.random_init_values(m_rng);
    init_augm_data(m_rng);
    update_global_loglike();
    reset_moves();

    // --- Initialize output file --- //

    // PRINT : column names and initial values to parameters output file
    init_param_file(output_file);

    // --- MCMC execution --- //
    if (m_setup.verbose)
    {
        std::cout << "\n=========MCMC:execution=========" << std::endl;
    }
    auto start_exec_time = high_resolution_clock::now(); // Get starting timepoint
    for (unsigned i = 0; i < m_setup.num_iter; i += m_setup.thinning)
    {
        // FLAG
        /*
        if (m_setup.verbose)
        {
            progress_bar(i, m_setup.num_iter);
        }
        */

        // --- MC sweeps over alla parameters to be inferred --- /
        for (unsigned j = 0; j < m_setup.thinning; j++)
        {
            // --- Change all parameters --- //
            this_iter = 1 + j + i;
            update_parameters();
            update_global_loglike(); // FIXME: I reassign loglike to each household

            // --- Data augmentation --- //
            for (unsigned k = 0; k < m_setup.num_augm_sweeps; k++)
            {
                update_augm_variables();
            }

            // --- Adaptative MCMC on parameters --- //
            if (this_iter < m_setup.num_adapt_iter)
            {
                adapt_ratio = (double)this_iter / m_setup.num_adapt_iter;
                for (std::string parID : glob_params.to_be_inferred)
                {
                    glob_params.adapt(parID, adapt_ratio);
                }
            }
        }
        // PRINT : parameters on parameter output file
        print_parameters(this_iter);
    }

    // close parameter output file
    close_param_file();

    if (m_setup.verbose)
    {
        std::cout << "\n===========MCMC:end=============" << std::endl;
        // print execution time
        auto stop_exec_time = high_resolution_clock::now();
        auto exec_duration = duration_cast<seconds>(stop_exec_time - start_exec_time);
        long int millisec_time = exec_duration.count();
        std::cout << "Execution time : " << millisec_time << " ms" << endl;
    }
}

//////////////////////////////////////////////////
// /////////////// MCMC : print /////////////// //
//////////////////////////////////////////////////

/* initialize the output param file */
void MCMC_Class::init_param_file(std::string output_filename)
{
    param_ofstream.open(output_filename.c_str());
    if (param_ofstream.is_open())
    {
        param_ofstream << "iter" << m_setup.delim << "loglike";
        if (m_setup.is_data_augm)
        {
            param_ofstream << m_setup.delim << "augm__a%";
        }
        for (std::string parID : glob_params.to_be_inferred)
        {
            param_ofstream << m_setup.delim << parID << m_setup.delim << parID << "__a%";
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

/* print the model parameters at the given iteration time */
void MCMC_Class::print_parameters(unsigned step)
{
    param_ofstream << step << m_setup.delim << get_global_loglike();
    if (m_setup.is_data_augm)
    {
        param_ofstream << m_setup.delim << 100.0 * ((double)m_num_accept_data / m_num_propos_data);
    }
    for (std::string parID : glob_params.to_be_inferred)
    {
        param_ofstream << m_setup.delim << glob_params.get_value(parID) << m_setup.delim << 100.0 * glob_params.get_accept_rate(parID);
    }
    param_ofstream << std::endl;
}

/* print the ensemble of augmented variables */
void MCMC_Class::print_augm_variables()
{
    for (size_t h = 0; h < num_households(); h++)
    {
        for (size_t res = 0; res < m_data[h].h_size; res++)
        {
            augm_ofstream << m_data[h].get_index() << m_setup.delim << m_data[h].h_size;
            augm_ofstream << m_setup.delim << m_data[h].get_index(res) << m_setup.delim << m_data[h].get_status(res);
            for (std::string augmID : m_data[h].get_augm_list(res))
            {
                augm_ofstream << m_setup.delim << m_data[h].get_augm(augmID, res);
            }
            augm_ofstream << std::endl;
        }
    }
}

/* open the file for the augmented variables */
void MCMC_Class::open_augm_file(std::string output_filename)
{
    augm_ofstream.open(output_filename.c_str());
    if (!augm_ofstream.is_open())
    {
        std::cerr << "\nERROR: Unable to open file `" << output_filename << "`" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }
}

/* progress bar */
void MCMC_Class::progress_bar(size_t iter, size_t n_iter)
{
    double progress = 100.0 * (1.0 + iter) / n_iter;

    // Write percentage of work done
    std::cout << "\rRunning... " << std::floor(progress) << "%" << std::flush;
}

//////////////////////////////////////////////////
// ////////////// MCMC : generate ///////////// //
//////////////////////////////////////////////////

/*
void MCMC_Class::generate_from(std::string output_file)

*/