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

    reset_moves();
}

// user
MCMC_Class::MCMC_Class(
    unsigned seed,
    std::vector<Household> data)
{
    m_rand_gen.seed(seed);
    m_data = data;

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
void MCMC_Class::init_augm_data(std::mt19937_64 &rand_gen)
{
    if (Verbose == true)
    {
        std::cout << "-> init augmented data;" << std::endl;
    }

    // loop over houses
    for (size_t h = 0; h < num_households(); h++)
    {
        m_data[h].init_augm_data(rand_gen);
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

// Update param value in the mcmc chain
bool MCMC_Class::update_param(std::string parID)
{
    // Random walk
    double old_value = glob_params.get_value(parID); // store old value
    double new_value, Q = 0.0;
    double log_ratio_proposal, log_ratio_prior, new_global_loglike;
    bool move_accepted = false;

    // propose new value
    new_value = glob_params.propose(parID, m_rand_gen, old_value);

    // proceed only if the new value is consistent
    if (glob_params.valid_param(parID, new_value) == true)
    {
        log_ratio_proposal = glob_params.log_ratio_proposal(parID, old_value, new_value);
        log_ratio_prior = glob_params.log_ratio_prior(parID, old_value, new_value);

        // Compute new log likelihood
        glob_params.set_value(parID, new_value);
        // FIXME (FUTURE) : infectivity profiles change if tost params change!!!

        new_global_loglike = 0.0;
        for (size_t h = 0; h < num_households(); h++)
        {
            new_global_loglike += m_data[h].compute_loglike();
        }

        // Metropolis-Hastings
        Q = new_global_loglike - m_global_loglike + log_ratio_prior + log_ratio_proposal;

        // FLAG

        if (parID == "alpha")
        {
            std::cerr << "alpha old : " << old_value << "| new : " << new_value << endl;
            std::cerr << "DLL : " << new_global_loglike - m_global_loglike << " | log_ratio_prior : " << log_ratio_prior << " | log_ratio_proposal : " << log_ratio_proposal << endl;
            std::cerr << endl;
        }
        // END-of-FLAG*/

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

// Update augmented variables
void MCMC_Class::update_augm_variables()
{
    double old_value, new_value, Q;
    double log_ratio_proposal = 0.0, log_ratio_prior = 0.0;
    double old_loglike_of_HH = 0.0, new_loglike_of_HH = 0.0;

    // Loop over households
    for (size_t h = 0; h < m_data.size(); h++)
    {
        // Loop over infectious residents of the househould
        for (size_t res : m_data[h].h_symptomatic_list)
        {
            for (std::string augmID : m_data[h].get_augm_list(res))
            {
                old_value = m_data[h].get_augm(augmID, res); // store old value
                old_loglike_of_HH = m_data[h].compute_loglike();

                new_value = m_data[h].propose_augm(augmID, res, m_rand_gen);

                // proceed only if the new value is consistent
                // FIXME (FUTURE) : is this necessary? this is basically always true
                if (m_data[h].valid_augm(augmID, res, new_value) == true)
                {
                    // Update onset time with proposed value
                    m_data[h].set_augm(augmID, res, new_value);
                    new_loglike_of_HH = m_data[h].compute_loglike();
                    log_ratio_prior = m_data[h].log_ratio_prior_augm(augmID, res, new_value);
                    log_ratio_proposal = m_data[h].log_ratio_proposal_augm(augmID, res, new_value);
                    Q = new_loglike_of_HH - old_loglike_of_HH + log_ratio_proposal;

                    // FLAG
                    std::cerr << augmID << " (" << res << ") | old : " << old_value << " | new : " << new_value << endl;
                    std::cerr << "DLL : " << new_loglike_of_HH - old_loglike_of_HH << " | log_ratio_proposal : " << log_ratio_proposal << endl;
                    std::cerr << endl;
                    // END-of-FLAG

                    if (unif_rand(m_rand_gen) < std::exp(Q))
                    { // move accepted
                        m_global_loglike += new_loglike_of_HH - old_loglike_of_HH;
                        m_num_accept_data++;
                    }
                    else
                    { // move refused
                        m_data[h].set_augm(augmID, res, old_value);
                    }
                }
                m_num_propos_data++;
            }
        }
    }
}

//////////////////////////////////////////////////
// /////////////// MCMC : infer /////////////// //
//////////////////////////////////////////////////

void MCMC_Class::infer(
    std::string output_file,
    int num_iter,
    int num_param_sweeps,
    int num_augm_sweeps,
    int num_adapt_steps,
    std::string tab)
{
    ofstream output(output_file.c_str());
    int step;
    bool move_accepted;
    bool is_adaptive_mcmc = num_adapt_steps > 0;
    bool is_data_augmentation = num_augm_sweeps > 0;

    // Column names of the output file
    std::string colNames = "iter" + tab + "loglike";
    if (is_data_augmentation == true)
    {
        colNames += tab + "augm__a%";
    }

    for (std::string parID : glob_params.to_be_inferred)
    {
        colNames += tab + parID + tab + parID + "__a%";
    }
    output << colNames << std::endl;

    std::cout << "\n=============MCMC===============" << std::endl;
    std::cout << "\n======MCMC:initialization=======" << std::endl;

    init_param_to_infer(); // NOTE: shouldn't parameters be initialized to given value?
    init_augm_data(m_rand_gen);
    init_global_loglike();

    output << 0.0 << tab << get_global_loglike();
    if (is_data_augmentation == true)
    {
        output << tab << 0.0; // n. of proposed/accepted move
    }
    for (std::string parID : glob_params.to_be_inferred)
    {
        output << tab << glob_params.get_value(parID) << tab << 0.0;
    }
    output << std::endl;

    std::cout << "\n=========MCMC:execution=========" << std::endl;
    for (int i = 0; i < num_iter; i++)
    {
        std::cout << "-> running iteration : " << (i + 1) << "/" << num_iter << endl;
        reset_moves();

        // Data augmentation
        for (int k = 0; k < num_augm_sweeps; k++)
        { // Update loglik at the household level
            update_augm_variables();
        }

        for (int j = 0; j < num_param_sweeps; j++)
        {
            step = 1 + j + i * num_param_sweeps;

            // Parameter inference
            for (std::string parID : glob_params.to_be_inferred)
            {
                move_accepted = update_param(parID);
                glob_params.update_accept_rate(parID, move_accepted, step);
                // Adaptative MCMC
                if (is_adaptive_mcmc == true && step <= num_adapt_steps)
                {
                    glob_params.adapt(parID, step, num_adapt_steps);
                }
            }
        }

        // Print output file
        output << (i + 1) << tab << get_global_loglike();
        if (is_data_augmentation == true)
        {
            output << tab << 100.0 * ((double)m_num_accept_data / m_num_propos_data);
        }
        for (std::string parID : glob_params.to_be_inferred)
        {
            output << tab << glob_params.get_value(parID) << tab << 100.0 * glob_params.get_accept_rate(parID);
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
void MCMC_Class::generate(std::string output_file, double temperature, int num_param_sweeps, std::string tab)
*/