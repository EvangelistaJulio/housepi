// mcmc.cpp //

#include "mcmc.h"

using namespace std;

extern bool Verbose;

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : Setter ///////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables */
void MCMC::initialize_model(bool is_rand_params)
{
    if (is_rand_params)
    {
        init_rand_params(rng);
    }

    assign_augm_data(rng);
    initialize_infectivity_profile();
}

/* Compute global loglikelihood */
// UPDATE : emberassing parallel
double MCMC::compute_global_loglike()
{
    double output = 0.0, tmp;
    {
        for (size_t h = 0; h < num_househ; h++)
        {
            tmp = compute_loglike(h);
            househ[h].set_loglike(tmp);
            output += tmp;
        }
    }
    return output;
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
            new_global_loglike = compute_global_loglike();

            //=====Metropolis-Hastings=====//
            Q = inv_temper * (new_global_loglike - get_global_loglike());
            Q += get_move_log_ratio();
            if (std::log(unif(rng)) < Q)
            {
                accept_move(parID);
                set_global_loglike(new_global_loglike);
            }
            else
            {
                refuse_move(parID);
                // restore loglike in households
                for (size_t h = 0; h < num_househ; h++)
                {
                    househ[h].refuse();
                }
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

    {
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
                            househ[h].set_loglike(new_loglike_of_HH);
                            set_global_loglike(get_global_loglike() + new_loglike_of_HH - old_loglike_of_HH);
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
// ////////////////////// MCMC : inference ////////////////////////// //
////////////////////////////////////////////////////////////////////////

void MCMC::inference(std::string output_file, bool is_rand_params, double inv_temper)
{
    unsigned step = 0;
    double adapt_ratio;

    if (Verbose)
    {
        std::cout << "\n=============MCMC===============" << std::endl;
        std::cout << "\n======MCMC:initialization=======" << std::endl;
    }

    initialize_model(is_rand_params);
    set_global_loglike(compute_global_loglike());
    reset_moves();

    //=====Initialize output file=====//
    init_param_file(output_file);

    //=====MCMC execution=====//
    if (Verbose)
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

            //=====Data augmentation=====//
            // FLAG : tests with frozen augmented variables
            // move_augm_variables(inv_temper);

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
// ///////////////////// MCMC : simulation ////////////////////////// //
////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : print ////////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Initialize the output param file with column names and initial values */
void MCMC::init_param_file(std::string output_filename)
{
    param_ofstream.open(output_filename.c_str());

    if (param_ofstream.is_open())
    {
        if (Verbose)
        {
            std::cout << "Inferred parameter file: " << output_filename << std::endl;
        }

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

/* Print the ensemble of augmented variables */
void MCMC::print_augm_variables(std::string output_filename)
{
    std::ofstream augm_ofstream;
    augm_ofstream.open(output_filename.c_str());

    if (augm_ofstream.is_open())
    {
        if (Verbose)
        {
            std::cout << "Augmented variable file: " << output_filename << std::endl;
        }

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

        for (size_t i = 0; i < num_indiv; i++)
        {
            augm_ofstream << indiv[i].get_index() << setup.delim << indiv[i].get_status();
            for (std::string dataID : data_headers)
            {
                augm_ofstream << setup.delim << indiv[i].get_data(dataID);
            }
            for (std::string augmID : augm_list)
            {
                augm_ofstream << setup.delim << indiv[i].get_augm(augmID);
            }
            augm_ofstream << std::endl;
        }

        augm_ofstream.close();
    }
    else
    {
        std::cerr << "\nERROR: Unable to open file `" << output_filename << "`" << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
    }
}

/* Add individual */ // FIXME: base model
void MCMC::add_individual(Individual new_indiv)
{
    int status = new_indiv.get_status();

    if ((status != HEALTHY) && (status != SYMPTOMATIC))
    {
        std::cerr << "ERROR: status not in the model: " << status << endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    indiv.push_back(new_indiv);
    num_indiv++;
}

/* Assign community structure to individuals */
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
        load_househ.at(h).add_individual(i, indiv[i].get_status());
    }

    //=====produce list of households=====//
    househ.clear();
    for (const auto &pair : load_househ)
    {
        househ.push_back(pair.second);
    }
    num_househ = househ.size();
    load_househ.clear();

    //=====assign home to individuals=====//
    num_symptomatic = 0;
    for (size_t h = 0; h < num_househ; h++)
    {
        size = 0;
        for (size_t i : househ[h].symptomatic_list)
        {
            indiv[i].home = h;
            num_symptomatic += 1;
            size++;
        }
        for (size_t i : househ[h].healthy_list)
        {
            indiv[i].home = h;
            size++;
        }
        househ[h].size = size;
    }
}