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
    for (std::string augmID : augm_to_be_inferred)
    {
        m_num_propos_data[augmID] = 0;
        m_num_accept_data[augmID] = 0;
    }

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
                            m_num_accept_data[augmID]++;
                        }
                        else
                        {
                            indiv[res].set_augm(augmID, old_val); // UPDATE: should be method Model
                            update_first_case(h);                 // Update information of first cases
                            update_infectivity_profiles(res);     // Update infectivity profiles related to idx
                        }
                    }
                    m_num_propos_data[augmID]++;
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
        for (std::string augmID : augm_to_be_inferred)
        {
            param_ofstream << setup.delim << "augm_" << augmID << "__a%";
        }
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
        for (std::string augmID : augm_to_be_inferred)
        {
            param_ofstream << setup.delim << 100.0 * ((double)m_num_accept_data[augmID] / m_num_propos_data[augmID]);
        }
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

/////////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : load  data ////////////////////////// //
/////////////////////////////////////////////////////////////////////////////

/* Add individual */
void MCMC::add_individual(Individual new_indiv)
{
    int status = new_indiv.get_status();

    //=====check for the validity of the status=====//
    if ((status != HEALTHY) && (status != SYMPTOMATIC))
    {
        std::cerr << "ERROR: status not in the model: " << status << endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    else
    {
        indiv[num_indiv] = new_indiv;
        num_indiv++;
    }
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
    num_househ = load_househ.size();    // get number of unique household index
    househ = new Household[num_househ]; // initialize the array of household
    h = 0;
    for (const auto &pair : load_househ)
    {
        househ[h] = pair.second;
        h++;
    }
    load_househ.clear(); // clear the dictionary for the loading

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

/* Read Data file */
void MCMC::build_data(std::string data_pathtofile)
{
    ifstream input(data_pathtofile.c_str());
    std::vector<Individual> aux_individuals = {};

    if (input.is_open())
    {
        if (Verbose)
        {
            cout << "\n=============DATA===============" << endl;
            cout << "Data file: " << data_pathtofile << endl;
        }

        //=====Read column names for variables in the model=====//
        std::string line;
        getline(input, line);
        istringstream iss(line);
        std::string token;
        unordered_map<string, double> data = {};

        // Tokenize the line using the specified delimiter
        while (getline(iss, token, setup.delim))
        {
            string_cleaner(&token);
            data_headers.push_back(token);
            data.emplace(token, 0);
        }

        //=====Load individuals from the raw with values=====//
        for (std::string line; getline(input, line);)
        {
            istringstream iss(line);

            // Store information in variables
            size_t n_col = 0;
            while (getline(iss, token, setup.delim))
            {
                data.at(data_headers[n_col]) = stod(token);
                n_col++;
            }

            // Create an onbject Individual and add it to the model
            Individual this_indiv(data.at("status"), data.at("index"), data);
            aux_individuals.push_back(this_indiv);
        }

        //=====assign individuals to mcmc=====//
        indiv = new Individual[aux_individuals.size()]; // initialize vector of individuals
        for (Individual this_indiv : aux_individuals)
        {
            add_individual(this_indiv);
        }
        aux_individuals.clear();

        //=====assign structure of the community=====//
        assign_community_structure("HH_idx");

        //=====Display=====//
        if (Verbose)
        {
            cout << "N. of individuals : " << num_indiv << endl;
            cout << "N. of households : " << num_househ << endl;
            cout << "N. of cases : " << num_symptomatic << endl;
        }
    }
    else
    {
        cerr << "\nERROR: Impossible to open: " << data_pathtofile << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
}

/////////////////////////////////////////////////////////////////////////////
// ////////////////////// MCMC : load  parameters //////////////////////// //
/////////////////////////////////////////////////////////////////////////////

/* Read Model Parameters file */
void MCMC::build_parameters(std::string data_pathtofile)
{
    char delim = setup.delim;
    char dict_delim = setup.dict_delim;
    std::ifstream input(data_pathtofile.c_str());

    if (input.is_open())
    {
        if (Verbose)
        {
            cout << "\n=============PARAMETERS==============" << endl;
            cout << "Model file: " << data_pathtofile << endl;
        }

        // Create a stringstream of the current line
        for (std::string line; getline(input, line);)
        {
            std::string name, value_str, infer, strategy, rate_str, range_str, prior_name, prior_args_str;
            double value;
            kwargs prior_args = {}, rate = {}, range = {};
            std::istringstream iss(line);

            //=====NAME=====//
            if (!getline(iss, name, delim))
            {
                build_column_error("name", line);
            }

            //=====VALUE=====//
            if (!getline(iss, value_str, delim))
            {
                build_column_error("value", line);
            }
            value = stod(value_str);

            //=====INFER=====//
            if (!getline(iss, infer, delim))
            {
                build_column_error("infer", line);
            }

            if (stobool(infer) == true)
            {
                //=====PRIOR=====//
                if (!getline(iss, prior_name, delim))
                {
                    build_column_error("prior", line);
                }

                //=====ARGS=====//
                if (!getline(iss, prior_args_str, delim))
                {
                    build_column_error("args", line);
                }
                prior_args = parse_dict_string(prior_args_str, dict_delim);

                //=====STRATEGY=====//
                if (!getline(iss, strategy, delim))
                {
                    build_column_error("strategy", line);
                }

                //=====STRATEGY_ARGS=====//
                if (!getline(iss, rate_str, delim))
                {
                    build_column_error("strat_args", line);
                }
                rate = parse_dict_string(rate_str, dict_delim);

                //=====RANGE=====//
                if (!getline(iss, range_str, delim))
                {
                    build_column_error("range", line);
                }
                range = parse_dict_string(range_str, dict_delim);

                // add entry to shared parameter dictionary
                add_param(name, value, prior_name, prior_args, strategy, rate, range);
            }
            else
            {
                // add entry to shared parameter dictionary
                add_param(name, value);
            }
        }

        if (Verbose)
        {
            display_parameters();
        }
    }
    else
    {
        cerr << "\nERROR: Impossible to read `" << data_pathtofile << "`" << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
}