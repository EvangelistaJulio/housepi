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
    m_num_iter = 0;
    m_num_augm_iter = 0;
    m_hh_loglike = {};
    m_data = {};
    m_num_adapt_iter = 0;

    init_param_to_infer();
    reset_moves();
}

// user
MCMC_Class::MCMC_Class(
    unsigned seed,
    int num_iter,
    int num_adapt_iter,
    std::vector<Household> data,
    int num_augm_iter)
{

    // Copy input vectors
    m_rand_gen.seed(seed);
    m_hh_loglike.resize(data.size());
    m_num_iter = num_iter;
    m_num_adapt_iter = num_adapt_iter;
    m_num_augm_iter = num_augm_iter;
    m_data = data;

    init_param_to_infer();
    reset_moves();
}

// ----------------------------
//  Methods : Setter
// ----------------------------

// Initialize param values

// FIXME : initial value should be the given one, no?
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

// Initialize infection time
void MCMC_Class::init_augm_data()
{
    if (Verbose == true)
    {
        std::cout << "-> init augmented data;" << std::endl;
    }

    // loop over houses
    for (size_t h = 0; h < num_households(); h++)
    {
        m_data[h].model_init_augm_data(m_rand_gen);
    }
}

// Initialize infection time
void MCMC_Class::init_loglike()
{
    if (Verbose == true)
    {
        std::cout << "-> init loglikelihood;" << std::endl;
    }

    for (size_t house = 0; house < num_households(); house++)
    {
        m_hh_loglike[house] = m_data[house].compute_loglike();
        m_global_loglike += m_hh_loglike[house];
    }
}

// ----------------------------
//  Methods : Computations
// ----------------------------

// Update param value in the mcmc chain
bool MCMC_Class::update_param(std::string parID)
{
    // Random walk
    double old_value = glob_params.value(parID); // store old value
    double new_value, Q;
    double log_ratio_proposal, log_ratio_prior, new_global_loglike;
    std::vector<double> new_loglike_HH(num_households(), 0.0);
    bool move_accepted = false;

    // propose new value
    new_value = glob_params.propose_move(parID, m_rand_gen);

    // proceed only if the new value is within range, i. e. no infinities
    if (glob_params.within_range(parID, new_value) == true)
    {
        log_ratio_proposal = glob_params.log_ratio_proposal(parID);
        log_ratio_prior = glob_params.log_ratio_prior(parID, new_value);

        // Compute new log likelihood
        glob_params.set_value(parID, new_value);
        new_global_loglike = 0.0;
        for (size_t house = 0; house < num_households(); house++)
        {
            new_loglike_HH[house] = m_data[house].compute_loglike();
            new_global_loglike += new_loglike_HH[house];
        }

        // FLAG
        std::cout << "(old LL) " << m_global_loglike << " | (new LL) " << new_global_loglike << endl;
        // END-of-FLAG

        // Metropolis
        Q = new_global_loglike - m_global_loglike + log_ratio_proposal + log_ratio_prior;

        if (unif_rand(m_rand_gen) < std::exp(Q))
        { // move accepted
            move_accepted = true;
            glob_params.set_num_accept_pp(parID);
            // Update log likelihood
            m_global_loglike = new_global_loglike;
            m_hh_loglike = new_loglike_HH;
        }
        else
        { // move refused
            glob_params.set_value(parID, old_value);
        }
    }
    glob_params.set_num_propos_pp(parID);

    return move_accepted;
}

// Update augmented infection times
void MCMC_Class::update_augm_infec_times()
{
    double old_value, new_value, Q;
    double old_proposal, new_proposal;
    double old_loglike_of_HH, new_loglike_of_HH;

    // Loop over households
    for (size_t h = 0; h < m_data.size(); h++)
    {
        // Loop over infected residents within house
        for (size_t i = 0; i < m_data[h].h_size; i++)
        {
            old_value = m_data[h].get_augm_infec_time(i);
            old_proposal = m_data[h].log_incub_period_pdf(i);
            old_loglike_of_HH = m_hh_loglike[h];

            // Update infection time with proposed value
            new_value = m_data[h].propose_augm_infec_time(i, m_rand_gen);
            m_data[h].set_augm_infec_time(i, new_value);
            new_proposal = m_data[h].log_incub_period_pdf(i);
            new_loglike_of_HH = m_data[h].compute_loglike();

            Q = new_loglike_of_HH - old_loglike_of_HH + old_proposal - new_proposal;

            if (unif_rand(m_rand_gen) < std::exp(Q))
            {
                m_global_loglike += new_loglike_of_HH - old_loglike_of_HH;
                m_hh_loglike[h] = new_loglike_of_HH;
                m_num_accept_data++;
            }
            else
            {
                m_data[h].set_augm_infec_time(i, old_value);
            }

            m_num_propos_data++;
        }
    }
}

// Update augmented symptom onsets
// UPDATE : these two loops for update augmented should be integrated...
void MCMC_Class::update_augm_symptom_onset()
{
    double old_value, new_value, Q;
    double old_proposal, new_proposal;
    double old_loglike_of_HH, new_loglike_of_HH;

    // Loop over households
    for (size_t h = 0; h < m_data.size(); h++)
    {
        // Loop over infected residents of the househould
        for (size_t i = 0; i < m_data[h].h_size; i++)
        {
            old_value = m_data[h].get_augm_onset_time(i);
            old_proposal = 0.0; // unif
            old_loglike_of_HH = m_hh_loglike[h];

            // Update onset time with proposed value
            new_value = m_data[h].propose_augm_onset_time(i, m_rand_gen);
            m_data[h].set_augm_onset_time(i, new_value);
            new_proposal = 0.0; // unif
            new_loglike_of_HH = m_data[h].compute_loglike();

            Q = new_loglike_of_HH - old_loglike_of_HH + old_proposal - new_proposal;

            if (unif_rand(m_rand_gen) < std::exp(Q))
            {
                m_global_loglike += new_loglike_of_HH - old_loglike_of_HH;
                m_hh_loglike[h] = new_loglike_of_HH;
                m_num_accept_data++;
            }
            else
            {
                m_data[h].set_augm_onset_time(i, old_value);
            }

            m_num_propos_data++;
        }
    }
}

//----------------------------------------------------------------------
// MCMC : infer
//----------------------------------------------------------------------

void MCMC_Class::infer(std::string output_file, int num_sweeps, std::string tab)
{

    ofstream output(output_file.c_str());

    int step;
    bool move_accepted;
    int num_augm_iter = get_num_augm_iter();

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
    init_loglike();

    output << 0 << tab << get_global_loglike(); // Log likelihood
    output << tab << 0 << tab << 0;             // n. of proposed/accepted move
    for (std::string parID : glob_params.to_be_inferred)
    {
        output << tab << glob_params.value(parID) << tab << 0 << tab << 0;
    }
    output << std::endl;

    std::cout << "\n=========MCMC:execution=========" << std::endl;
    for (int i = 0; i < get_num_iter(); i++)
    {
        std::cout << "-> running iteration : " << (i + 1) << "/" << get_num_iter() << endl;
        reset_moves();

        for (int j = 0; j < num_sweeps; j++)
        {
            for (std::string parID : glob_params.to_be_inferred)
            {
                step = i * num_sweeps + j + 1;
                move_accepted = update_param(parID);

                // Adaptative MCMC
                glob_params.adapt_rate(parID, move_accepted, step, m_num_adapt_iter);
            }

            // Augmented data
            for (int k = 0; k < num_augm_iter; k++)
            {
                update_augm_infec_times(); // Update loglik at the household level
                update_augm_symptom_onset();
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