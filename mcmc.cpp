// mcmc.cpp //

#include "mcmc.hpp"

using namespace std;
using namespace std::chrono;

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : Setter ///////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Set accept rates to zero */
void MCMC::reset_moves()
{
    m_num_propos_data = 0;
    m_num_accept_data = 0;
    for (std::string parID : to_be_inferred)
    {
        reset_accept_rate(parID);
    }
}

/* Initialize augmented variables */
void MCMC::declare_augm_data(std::mt19937_64 &rng)
{
    for (size_t h = 0; h < num_households(); h++)
    {
        init_infectivity_profiles(&m_data[h]);
        assign_augm_data(&m_data[h], rng);
    }
}

/* Initialize global loglikelihood */
// UPDATE : this initialization is used because changes are not stored for params update
void MCMC::update_global_loglike()
{
    // #pragma omp parallel num_threads(NUM_THREAD)
    {
        m_global_loglike = 0.0;
        // #pragma omp for reduction(+ : sum)
        for (size_t h = 0; h < num_households(); h++)
        {
            m_data[h].h_loglike = compute_loglike(&m_data[h]);
            m_global_loglike += m_data[h].h_loglike;
        }
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
    for (std::string parID : to_be_inferred)
    {
        // proceed only if the new value is consistent
        if (propose_valid_move(parID, m_rng))
        {
            // Compute the new log likelihood
            // UPDATE : emberassing parallel
            new_global_loglike = 0.0;
            for (size_t h = 0; h < num_households(); h++)
            {
                new_global_loglike += compute_loglike(&m_data[h]);
            }

            //=====Metropolis-Hastings=====//
            Q = inv_temper * (new_global_loglike - m_global_loglike);
            Q += get_log_ratio();
            if (std::log(m_unif(m_rng)) < Q)
            {
                accept_move(parID);
                m_global_loglike = new_global_loglike;
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
    bool move_accepted = false;
    double old_val, new_val, Q;
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
                    old_val = m_data[h].get_augm(augmID, res);            // store old value
                    old_loglike_of_HH = m_data[h].h_loglike;              // store old loglike
                    new_val = m_data[h].augm_propose(augmID, res, m_rng); // propose new value

                    // FIXME : always valid if prior proposer
                    // proceed only if the new value is consistent
                    move_accepted = false;
                    if (m_data[h].augm_valid(augmID, res, new_val))
                    {                                                 // Update with proposed value to perform computations
                        m_data[h].set_augm(augmID, res, new_val);     // UPDATE: should be method Model
                        update_first_case(&m_data[h]);                // Update information of first cases
                        update_infectivity_profiles(&m_data[h], res); // Update infectivity profiles related to idx

                        new_loglike_of_HH = compute_loglike(&m_data[h]);
                        log_ratio_proposal = m_data[h].augm_propose_log_ratio(augmID, res, old_val, new_val);

                        Q = inv_temper * (new_loglike_of_HH - old_loglike_of_HH);
                        Q -= log_ratio_proposal;

                        if (std::log(m_unif(m_rng)) < Q)
                        {
                            move_accepted = true;
                            m_global_loglike += new_loglike_of_HH - old_loglike_of_HH;
                            m_data[h].h_loglike = new_loglike_of_HH;
                            m_num_accept_data++; // NOTE : problematic for normalization
                        }
                        else
                        {
                            m_data[h].set_augm(augmID, res, old_val);     // UPDATE: should be method Model
                            update_first_case(&m_data[h]);                // Update information of first cases
                            update_infectivity_profiles(&m_data[h], res); // Update infectivity profiles related to idx
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

void MCMC::run(std::string output_file, double inv_temper)
{
    unsigned step;
    double adapt_ratio;

    if (m_setup.verbose)
    {
        std::cout << "\n=============MCMC===============" << std::endl;
        std::cout << "\n======MCMC:initialization=======" << std::endl;
    }

    init_rand_params(m_rng);
    declare_augm_data(m_rng);
    update_global_loglike();
    reset_moves();

    //=====Initialize output file=====//
    init_param_file(output_file);

    //=====MCMC execution=====//
    if (m_setup.verbose)
    {
        std::cout << "\n=========MCMC:execution=========" << std::endl;
    }
    auto start_exec_time = high_resolution_clock::now(); // Get starting timepoint
    for (unsigned i = 0; i < m_setup.num_iter; i += m_setup.thinning)
    {
        //=====MC sweeps over alla parameters to be inferred=====//
        for (unsigned j = 0; j < m_setup.thinning; j++)
        {
            //=====Change all parameters=====//
            step = 1 + j + i;
            move_parameters(inv_temper);
            update_global_loglike(); // FIXME: I reassign loglike to each household

            //=====Data augmentation=====// 
            if (m_setup.is_data_augm)
            {
                move_augm_variables(inv_temper);
            }

            //=====Adaptative MCMC on parameters=====//
            if (step < m_setup.num_adapt_iter)
            {
                adapt_ratio = (double)step / m_setup.num_adapt_iter;
                for (std::string parID : to_be_inferred)
                {
                    adapt(parID, adapt_ratio);
                }
            }
        }
        print_parameters(step);
    }
    
    //=====MCMC ending=====//
    close_param_file();

    if (m_setup.verbose)
    {
        std::cout << "\n===========MCMC:end=============" << std::endl;
        // print execution time
        auto stop_exec_time = high_resolution_clock::now();
        auto exec_duration = duration_cast<seconds>(stop_exec_time - start_exec_time);
        long int millisec_time = exec_duration.count();
        std::cout << "Execution time : " << millisec_time << " ms\n"
                  << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : print ////////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Initialize the output param file with column names and initial values */
void MCMC::init_param_file(std::string output_filename)
{
    param_ofstream.open(output_filename.c_str());
    if (param_ofstream.is_open())
    {
        param_ofstream << "iter" << m_setup.delim << "loglike";
        if (m_setup.is_data_augm)
        {
            param_ofstream << m_setup.delim << "augm__a%";
        }
        for (std::string parID : to_be_inferred)
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

/* Print the model parameters at the given iteration time */
void MCMC::print_parameters(unsigned step)
{
    param_ofstream << step << m_setup.delim << get_global_loglike();
    if (m_setup.is_data_augm)
    {
        param_ofstream << m_setup.delim << 100.0 * ((double)m_num_accept_data / m_num_propos_data);
    }
    for (std::string parID : to_be_inferred)
    {
        param_ofstream << m_setup.delim << get_value(parID) << m_setup.delim << 100.0 * get_accept_rate(parID);
    }
    param_ofstream << std::endl;
}

/* Print the ensemble of augmented variables */
void MCMC::print_augm_variables()
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

/* Open the file for the augmented variables */
void MCMC::open_augm_file(std::string output_filename)
{
    augm_ofstream.open(output_filename.c_str());
    if (!augm_ofstream.is_open())
    {
        std::cerr << "\nERROR: Unable to open file `" << output_filename << "`" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }
}