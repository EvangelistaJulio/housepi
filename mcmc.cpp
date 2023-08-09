// mcmc.cpp //

#include "mcmc.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  MCMC_Class : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

// default
MCMC_Class::MCMC_Class()
{
    m_rand_gen.seed(DEFAULT_SEED);
    m_global_loglike = 0.0;
    m_data = {};

    // init_param_to_infer(); // NOTE: parameters are initialized to given value
    reset_moves();
}

// user
MCMC_Class::MCMC_Class(
    unsigned seed,
    std::vector<Household> data)
{
    m_rand_gen.seed(seed);
    m_data = data;

    // init_param_to_infer(); // NOTE: parameters are initialized to given value
    reset_moves();
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Initialize param values generating them from the prior
void MCMC_Class::init_param_to_infer()
{
    for (std::string parID : glob_params.to_be_inferred)
    {
        glob_params.init_value(parID, m_rand_gen);
    }
}

// Modification of attributes
void MCMC_Class::reset_moves()
{
    m_num_propos_data = 0;
    m_num_accept_data = 0;
    for (std::string parID : glob_params.to_be_inferred)
    {
        glob_params.reset_moves(parID);
    }
}

/* Initialize augmented variables */
void MCMC_Class::init_augm_data()
{
    if (Verbose == true)
    {
        std::cout << "-> init augmented data;" << std::endl;
    }

    // loop over houses
    for (size_t h = 0; h < num_households(); h++)
    {
        m_data[h].init_augm_data(m_rand_gen);
    }
}

/* Initialize global loglikelihood */
void MCMC_Class::init_global_loglike()
{
    if (Verbose == true)
    {
        std::cout << "-> init loglikelihood;" << std::endl;
    }

    m_global_loglike = 0.0;
    // loop over houses
    for (size_t h = 0; h < num_households(); h++)
    {
        m_global_loglike += m_data[h].compute_loglike();
    }
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// FIXME : infectivity profiles change if tost params change!!!
// Update param value in the mcmc chain
bool MCMC_Class::update_param(std::string parID)
{
    // Random walk
    double old_value = glob_params.value(parID); // store old value
    double new_value, Q;
    double log_ratio_proposal, log_ratio_prior, new_global_loglike;
    bool move_accepted = false;

    // propose new value
    new_value = glob_params.propose(parID, m_rand_gen);

    // proceed only if the new value is within range, i. e. no infinities
    if (glob_params.within_range(parID, new_value) == true)
    {
        log_ratio_proposal = glob_params.log_ratio_proposal(parID);
        log_ratio_prior = glob_params.log_ratio_prior(parID, new_value);

        // Compute new log likelihood
        glob_params.set_value(parID, new_value);
        new_global_loglike = 0.0;
        for (size_t h = 0; h < num_households(); h++)
        {
            new_global_loglike += m_data[h].compute_loglike();
        }

        // Metropolis-Hastings
        Q = new_global_loglike - m_global_loglike + log_ratio_proposal + log_ratio_prior;

        if (unif_rand(m_rand_gen) < std::exp(Q))
        { // move accepted
            move_accepted = true;
            glob_params.set_num_accept_pp(parID);
            // Update log likelihood
            m_global_loglike = new_global_loglike;
        }
        else
        { // move refused
            glob_params.set_value(parID, old_value);
        }
    }
    glob_params.set_num_propos_pp(parID);

    return move_accepted;
}

// FIXME no consistency check
// Update augmented variables
void MCMC_Class::update_augm_variables()
{
    double old_value, new_value, Q;
    double delta_proposal = 0.0;
    double delta_loglike_of_HH = 0.0;

    // Loop over households
    for (size_t h = 0; h < m_data.size(); h++)
    {
        // Loop over infectious residents of the househould
        for (size_t res : m_data[h].h_symptomatic_list)
        {
            // UPDATE : these two chunks to update augmented should be integrated...
            // --- onset time --- //
            old_value = m_data[h].get_augm_onset_time(res); // store old value
            delta_proposal = 0.0;                           // unif
            delta_loglike_of_HH = -m_data[h].compute_loglike();

            // Update onset time with proposed value
            new_value = m_data[h].propose_augm_onset_time(res, m_rand_gen);
            m_data[h].set_augm_onset_time(res, new_value);
            delta_proposal += 0.0; // unif
            delta_loglike_of_HH += m_data[h].compute_loglike();

            Q = delta_loglike_of_HH - delta_proposal;

            if (unif_rand(m_rand_gen) < std::exp(Q))
            {
                m_global_loglike += delta_loglike_of_HH;
                m_num_accept_data++;
            }
            else
            {
                m_data[h].set_augm_onset_time(res, old_value);
            }
            m_num_propos_data++;

            // --- incub period --- //
            old_value = m_data[h].get_augm_incub_period(res); // store old value
            delta_proposal = -m_data[h].log_incub_period_pdf(res);
            delta_loglike_of_HH = -m_data[h].compute_loglike();

            // Update incubation period time with proposed value
            new_value = m_data[h].propose_augm_incub_period(res, m_rand_gen);
            m_data[h].set_augm_incub_period(res, new_value);
            delta_proposal += m_data[h].log_incub_period_pdf(res);
            delta_loglike_of_HH += m_data[h].compute_loglike();

            Q = delta_loglike_of_HH - delta_proposal;

            if (unif_rand(m_rand_gen) < std::exp(Q))
            {
                m_global_loglike += delta_loglike_of_HH;
                m_num_accept_data++;
            }
            else
            {
                m_data[h].set_augm_incub_period(res, old_value);
            }
            m_num_propos_data++;
        }
    }
}

//////////////////////////////////////////////////
// /////////////// MCMC : infer /////////////// //
//////////////////////////////////////////////////

void MCMC_Class::infer(
    std::string output_file,
    int num_iter,
    int num_sweeps,
    int num_augm_sweeps,
    int num_adapt_iter,
    std::string tab)
{
    ofstream output(output_file.c_str());
    int step;
    bool move_accepted;

    // Column names of the output file
    std::string colNames = "iteration" + tab + "loglike";
    colNames += tab + "num_p" + tab + "num_a";
    for (std::string parID : glob_params.to_be_inferred)
    {
        colNames += tab + parID + tab + parID + "_p" + tab + parID + "_a";
    }
    output << colNames << std::endl;

    std::cout << "\n=============MCMC===============" << std::endl;
    std::cout << "\n======MCMC:initialization=======" << std::endl;

    init_augm_data();
    init_global_loglike();    

    output << 0 << tab << get_global_loglike(); // Log likelihood
    output << tab << 0 << tab << 0;             // n. of proposed/accepted move
    for (std::string parID : glob_params.to_be_inferred)
    {
        output << tab << glob_params.value(parID) << tab << 0 << tab << 0;
    }
    output << std::endl;

    std::cout << "\n=========MCMC:execution=========" << std::endl;
    for (int i = 0; i < num_iter; i++)
    {
        std::cout << "-> running iteration : " << (i + 1) << "/" << num_iter << endl;
        reset_moves();

        for (int j = 0; j < num_sweeps; j++)
        {
            for (std::string parID : glob_params.to_be_inferred)
            {
                step = i * num_sweeps + j + 1;
                move_accepted = update_param(parID);

                // Adaptative MCMC
                glob_params.adapt_rate(parID, move_accepted, step, num_adapt_iter);
            }

            // Augmented data
            for (int k = 0; k < num_augm_sweeps; k++)
            {
                update_augm_variables(); // Update loglik at the household level
            }
        }

        // Print output file
        output << (i + 1) << tab << get_global_loglike();               // Log likelihood
        output << tab << m_num_propos_data << tab << m_num_accept_data; // n. of proposed/accepted move
        for (std::string parID : glob_params.to_be_inferred)
        {
            output << tab << glob_params.value(parID) << tab << glob_params.num_propos(parID) << tab << glob_params.num_accept(parID);
        }
        output << std::endl;
    }

    output.close();
    std::cout << "\n===========MCMC:end=============\n"
              << std::endl;
}

//////////////////////////////////////////////////
// ////////////// MCMC : generate ///////////// //
//////////////////////////////////////////////////

/*
void MCMC_Class::generate(std::string output_file, double temperature, int num_sweeps, std::string tab)
*/