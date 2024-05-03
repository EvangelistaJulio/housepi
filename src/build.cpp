// build.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

// Warning for status
bool Unknown_Status_Warning = true;

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// MCMC : build extra ////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Extra File (default unused) */
void MCMC::build_extra(std::string data_pathtofile)
{
    UNUSED(data_pathtofile);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// MCMC : build setup ////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Setup file */
void MCMC::build_setup(std::string setup_pathtofile)
{
    if (setup_pathtofile != UNDEFINED)
    {
        char delim = ',';
        std::ifstream input(setup_pathtofile.c_str());
        ssargs safe_load = {};

        if (input.is_open())
        {
            display_decorator("SETUP");
            display_stdout("Loaing setup from: " + setup_pathtofile);

            // Create a std::stringstream of the current line
            for (std::string line; getline(input, line);)
            {
                std::string name, choice;
                std::istringstream iss(line);

                //=====NAME=====//
                if (!getline(iss, name, delim))
                {
                    std::string err_text = "Unable to extract setup name.\nCheck file line: " + line;
                    raise_error_message(err_text);
                    std::exit(EXIT_FAILURE);
                }

                //=====CHOICE=====//
                if (!getline(iss, choice))
                {
                    std::string err_text = "Unable to extract setup choice.\nCheck file line: " + line;
                    raise_error_message(err_text);
                    std::exit(EXIT_FAILURE);
                }

                // add the new setup
                safe_load.emplace(string_cleaner(name), string_cleaner(choice));
            }

            // MCMC number of iterations
            setup_num_iter = stoi(safe_load.at("num_iter"));
            setup_num_adapt_iter = stoi(safe_load.at("num_adapt_iter"));
            setup_thinning = stoi(safe_load.at("thinning"));

            // seed for pseudo-random number generator (default is time)
            std::string seed_str = safe_load.at("seed");
            if (seed_str != UNDEFINED)
            {
                setup_seed = stoi(seed_str);
            }
        }
        else
        {
            std::string err_text = "Impossible to read file: " + setup_pathtofile + ".";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
    }

    //=====Pseudo-Random Number Generator Seeding=====//
    prng_seed(setup_seed);

    //=====Display=====//
    display_decorator("SETUP::definitions", '-', 1);
    std::ostringstream oss;
    oss << "seed" << DELIM << setup_seed
        << "\nnum_iter" << DELIM << setup_num_iter
        << "\nnum_adapt_iter" << DELIM << setup_num_adapt_iter
        << "\nthinning" << DELIM << setup_thinning;
    if (setup_num_omp_threads > 1)
    {
        oss << "\nnum_threads" << DELIM << setup_num_omp_threads;
    }
    display_stdout(&oss);
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// MCMC : build parameters //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Parameters file */
void MCMC::build_parameters(std::string data_pathtofile)
{
    std::ifstream input(data_pathtofile.c_str());

    if (input.is_open())
    {
        display_decorator("PARAMETERS");
        display_stdout("Loading parameters definition from: " + data_pathtofile);

        //=====Load column headers=====//
        std::vector<std::string> headers = read_column_headers(input, DELIM);
        // check the validity of the headers
        for (std::string req : param_file_required_headers)
        {
            if (std::find(headers.begin(), headers.end(), req) == headers.end())
            {
                std::string err_text = "Missing column `" + req + "`. Please check file : " + data_pathtofile + ".";
                raise_error_message(err_text);
                std::exit(EXIT_FAILURE);
            }
        }

        //=====Load parameters from the row with options=====//
        for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row, DELIM); this_row.clear())
        {
            // add entry to shared parameter dictionary
            if (stobool(this_row.at("infer")))
            {
                add_parameter(
                    this_row["name"],
                    std::stod(this_row["value"]),
                    this_row["prior"],
                    parse_dict_args(this_row["prior_args"], DICT_DELIM),
                    this_row["strategy"],
                    parse_dict_args(this_row["strategy_args"], DICT_DELIM),
                    parse_dict_args(this_row["range"], DICT_DELIM));
            }
            else
            {
                add_parameter(this_row.at("name"), stod(this_row.at("value")));
            }
        }

        //=====display a summary of the used parameters=====//
        display_decorator("PARAMETERS::definitions", '-', 1);
        display_parameters();
    }
    else
    {
        std::string err_text = "Impossible to read file:" + data_pathtofile + ".";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

void MCMC::display_model_report()
{
    if (Verbose)
    {
        // check if user defined groups other than reference
        bool has_non_reference_groups = false;
        if (!model_report.empty())
        {
            for (const auto &entry : model_report)
            {
                if (entry.second.size() > 1)
                {
                    has_non_reference_groups = true;
                    break;
                }
            }
        }
        if (has_non_reference_groups)
        {
            std::ostringstream oss;
            display_decorator("DATABASE::groups", '-', 1);
            for (const auto &which : model_report)
            {
                oss << which.first << ":";
                for (auto const &x : which.second)
                {
                    oss << "\n    * " << x.first << ": " << x.second
                        << std::fixed << std::setprecision(2) << " (" << (100.0 * x.second / num_indiv) << " %)";
                }
            }
            display_stdout(&oss);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////// MCMC : build_data //////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read data file and pre-process the data*/
void MCMC::build_data(std::string data_pathtofile)
{
    std::ifstream input(data_pathtofile.c_str());
    std::vector<std::string> required_headers = {"index", "status", "HH_idx"};

    if (input.is_open())
    {
        display_decorator("DATABASE");
        display_stdout("Loading data from: " + data_pathtofile);

        //=====Load column headers=====//
        data_headers = read_column_headers(input, DELIM);
        // check the validity of the headers
        for (std::string req : required_headers)
        {
            if (std::find(data_headers.begin(), data_headers.end(), req) == data_headers.end())
            {
                std::string err_text = "Missing column `" + req + "`. Please check file : " + data_pathtofile + ".";
                raise_error_message(err_text);
                std::exit(EXIT_FAILURE);
            }
        }

        //=====Load individuals from the row with values=====//
        std::vector<Individual> aux_individuals = {};
        for (ssargs this_row = {}; read_row_by_headers(input, data_headers, this_row, DELIM); this_row.clear())
        {
            Individual this_indiv(row_to_sdargs(this_row));
            aux_individuals.push_back(this_indiv);
        }

        //=====assign individuals to mcmc=====//
        int num_cases = 0;
        int num_unknown_status = 0;
        indiv = new Individual[aux_individuals.size()];
        for (Individual this_indiv : aux_individuals)
        {
            add_individual(this_indiv);
            if (this_indiv.get_status() == INFECTIOUS)
            {
                num_cases++;
            }
            else if (this_indiv.get_status() != SUSCEPTIBLE)
            {
                num_unknown_status++;
            }
        }

        //=====assign structure of the community=====//
        assign_community_structure("HH_idx"); // FIXME : HH_idx is found where??

        // database analysis
        int min_house_size = 1000, max_house_size = 0;
        int no_house_cases = 0, atleast_one_case = 0, atleast_two_cases = 0;
        for (int h = 0; h < num_househ; h++)
        {
            min_house_size = std::min(househ[h].house_size, min_house_size);
            max_house_size = std::max(househ[h].house_size, max_house_size);
            int aux = househ[h].infectious_list.size();
            if (aux == 0)
            {
                no_house_cases++;
            }
            else
            {
                atleast_one_case++;
                if (aux > 1)
                {
                    atleast_two_cases++;
                }
            }
        }

        //=====initialize contact patterns=====//
        // UPDATE: this scale like individuls, but could scale as househodls
        // NOTE : however the indexes are tricky to be mapped back
        contact_patterns = new const double *[num_indiv * num_indiv];
        std::fill_n(contact_patterns, num_indiv * num_indiv, nullptr);

        // FIXME : missing report and control
        for (int h = 0; h < num_househ; h++)
        {
            for (int i : househ[h].get_house_members())
            {
                for (int k : househ[h].get_house_members())
                {
                    if (i != k)
                    {
                        std::string parID = assign_rate_per_contact(i, k);
                        if (parID == REFERENCE)
                        {
                            contact_patterns[from_to(i, k)] = get_reference();
                        }
                        else
                        {
                            contact_patterns[from_to(i, k)] = get_address(parID);
                        }
                    }
                }
            }
        }

        //=====Display General Information=====//
        if (Verbose)
        {
            std::ostringstream oss;

            //-----Display Database Summary-----//
            display_decorator("DATABASE::summary", '-', 1);
            oss << "N. of individuals: " << num_indiv
                << "\n    * N. of cases: " << num_cases
                << std::fixed << std::setprecision(2) << " (" << (100.0 * num_cases / num_indiv) << " %)";
            if (num_unknown_status > 0)
            {
                oss << "\nN. of missing status : " << num_unknown_status;
            }
            oss << "\nN. of households: " << num_househ;
            if (no_house_cases > 0)
            {
                oss << "\n    * with no cases: " << no_house_cases
                    << std::fixed << std::setprecision(2) << " (" << (100.0 * no_house_cases / num_househ) << " %)";
            }
            oss << "\n    * with at least 1 case: " << atleast_one_case
                << std::fixed << std::setprecision(2) << " (" << (100.0 * atleast_one_case / num_househ) << " %)"
                << "\n    * with at least 2 cases: " << atleast_two_cases
                << std::fixed << std::setprecision(2) << " (" << (100.0 * atleast_two_cases / num_househ) << " %)"
                << "\nMin n. house members: " << min_house_size
                << "\nMax n. house members: " << max_house_size;
            display_stdout(&oss);
            oss.str("");
            oss.clear();
        }
        display_model_report();
        model_report.clear();
    }
    else
    {
        std::string err_text = "Impossible to open file: `" + data_pathtofile + "`";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

//---------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////  MCMC : add individuals and context ///////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Add individual */
void MCMC::add_individual(Individual new_indiv)
{
    int status = new_indiv.get_status();
    int i = num_indiv; // index of the individual
    std::string rateID = UNDEFINED;
    std::string parID = UNDEFINED;

    //=====Store New Individual=====//
    indiv[i] = new_indiv;
    num_indiv++;

    //=====check for the validity of the status=====//
    if ((status != SUSCEPTIBLE) && (status != INFECTIOUS))
    {
        std::string warn_text = "Unknown status `" + std::to_string(status) + "`.";
        warn_text += " Use " + std::to_string(SUSCEPTIBLE) + "/" + std::to_string(INFECTIOUS) + " for healthy/non-healthy state.";
        warn_text += "\n If it was intentional, ignore this message.";
        raise_warning_message(warn_text, &Unknown_Status_Warning);
    }

    //=====Assign Model Rates=====//
    assign_rate_per_group_to_individual(i);

    for (const auto &which : indiv[i]._which_rate)
    {
        rateID = which.first;
        parID = which.second;

        if (num_indiv == 1)
        {
            // if first individual, initialize rateID-report
            model_report.emplace(rateID, std::map<std::string, int>());
            model_report.at(rateID).emplace("[reference]", 0);
        }
        else
        {
            // all individuals must have the same set of rates of first individual
            if (model_report.find(rateID) == model_report.end())
            {
                std::string err_text = "Rate `" + rateID + "` is not provided for individuals.";
                raise_error_message(err_text);
                std::exit(EXIT_FAILURE);
            }
        }

        // set rates and update report
        if (parID == UNDEFINED)
        {
            // error : undefined rate
            std::string err_text = "Rate `" + rateID + "` is not provided for individual " + std::to_string(indiv[i].get_index()) + ".";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
        else if (parID == REFERENCE)
        {
            // rate set to reference
            indiv[i]._set_rate(rateID, get_reference());
            model_report.at(rateID).at("[reference]") += 1;
        }
        else
        {
            // rate associated to parameter
            indiv[i]._set_rate(rateID, get_address(parID));
            if (model_report.at(rateID).find(parID) == model_report.at(rateID).end())
            {
                model_report.at(rateID).emplace(parID, 1);
            }
            else
            {
                model_report.at(rateID).at(parID) += 1;
            }
        }
    }
}

/* Assign group rate to an individual */
void MCMC::assign_rate_per_group_to_individual(int i)
{
    indiv[i].assign_rate_per_group("Infectivity", REFERENCE);
    indiv[i].assign_rate_per_group("Susceptibility", REFERENCE);
}

/* Assign community structure to individuals. I returns a report. */
void MCMC::assign_community_structure(std::string HH_index_column)
{
    int h;
    std::unordered_map<int, Household> load_househ = {};

    //=====safe load households=====//
    for (int i = 0; i < num_indiv; i++)
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
    num_househ = load_househ.size(); // get number of unique household indexes

    //=====produce list of households=====//
    househ = new Household[num_househ]; // initialize the MCMC array of households

    h = 0;
    for (const auto &pair : load_househ)
    {
        househ[h] = pair.second;
        househ[h].set_indiv_ptr(&indiv); // associate MCMC array of individuals to the Structure
        h++;
    }
    load_househ.clear(); // clear the loading dictionary

    //=====assign home to individuals=====//
    for (int h = 0; h < num_househ; h++)
    {
        for (int i : househ[h].get_house_members())
        {
            indiv[i].home = h;
        }
    }
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////  MCMC : build augmented //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file */
void MCMC::build_augmented(std::string data_pathtofile)
{
    if (data_pathtofile != UNDEFINED)
    {
        std::ifstream input(data_pathtofile.c_str());
        std::vector<std::string> required_headers = {
            "name", "value", "model", "model_args", "strategy", "strategy_args", "move_from", "range", "init_rand"};

        if (input.is_open())
        {
            display_decorator("AUGMENTED");
            display_stdout("Loading augmented definition from: " + data_pathtofile);

            //=====Load column headers=====//
            std::vector<std::string> headers = read_column_headers(input, DELIM);
            // check the validity of the headers
            for (std::string req : required_headers)
            {
                if (std::find(headers.begin(), headers.end(), req) == headers.end())
                {
                    std::string err_text = "Missing column `" + req + "`. Please check file : " + data_pathtofile + ".";
                    raise_error_message(err_text);
                    std::exit(EXIT_FAILURE);
                }
            }

            //=====Load augmented from the row with options=====//
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row, DELIM); this_row.clear())
            {
                load_augm_definition.push_back(this_row);
            }

            //=====assign augmented variables to all individuals=====//
            assign_augmented_variables();
        }
        else
        {
            std::string err_text = "Impossible to read file:" + data_pathtofile + ".";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////  MCMC : assign augmented variables /////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign augmented variables to all individuals. */
void MCMC::assign_augmented_variables()
{
    std::string augmID;

    // load the definitions of agumetned variablews and define lists of usage
    for (ssargs loader : load_augm_definition)
    {
        augmID = loader.at("name");
        augm_full_set.emplace(augmID);

        // list of augmented variables to be inferred
        if (loader.at("strategy") != UNDEFINED)
        {
            augm_to_be_inferred.emplace(augmID);

            // list of augmented variables with a model
            // note: no contribution of constant variables to Metropolis
            if (loader.at("model") != UNDEFINED)
            {
                augm_with_model_set.emplace(augmID);
                if (stobool(loader.at("init_rand")) == true)
                {
                    augm_rand_init_set.emplace(augmID);
                }
            }
        }
    }
    is_data_augmentation = !augm_to_be_inferred.empty();

    // assign augmented variables to each individual
    for (int i = 0; i < num_indiv; i++)
    {
        for (ssargs loader : load_augm_definition)
        {
            Augmented augm_var(
                type_to_double(loader.at("value"), i),
                loader.at("model"),
                typedict_to_sdargs_ext(loader.at("model_args")),
                loader.at("strategy"), 
                typedict_to_sdargs(loader.at("strategy_args"), i),
                type_to_double(loader.at("move_from"), i),
                typedict_to_sdargs(loader.at("range"), i));
            indiv[i].add_augm(loader.at("name"), augm_var);
        }
    }

    // free memory
    load_augm_definition.clear();
}

/* Initialize augmented variable for all individuals */
void MCMC::init_augmented_variables()
{
    for (std::string augmID : augm_rand_init_set)
    {
        for (int i = 0; i < num_indiv; i++)
        {
            indiv[i].init_augm(augmID, prng);
        }
    }
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////  MCMC : init param file //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void MCMC::init_param_file(std::string output_filename)
{
    param_ofstream.open(output_filename.c_str());

    if (param_ofstream.is_open())
    {
        display_stdout("Saving MCMC iterations to: " + output_filename);

        param_ofstream << "iter" << DELIM << "loglike";
        for (std::string augmID : augm_to_be_inferred)
        {
            param_ofstream << DELIM << "augm_" << augmID << "__a%";
        }
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << DELIM << parID << DELIM << parID << "__a%";
        }
        param_ofstream << std::endl;
        // print_parameters(0); // FIXME: meaningless if augmented steps are not saved
    }
    else
    {
        std::string err_text = "Unable to open file `" + output_filename + "`";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// MCMC : print parameters //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the model parameters at the given iteration time */
void MCMC::print_parameters(int step)
{
    if (param_ofstream.is_open())
    {
        param_ofstream << step << DELIM << std::setprecision(2) << std::fixed << get_global_loglike();
        for (std::string augmID : augm_to_be_inferred)
        {
            param_ofstream << DELIM << std::setprecision(2) << std::fixed << 100.0 * get_augm_accept_rate(augmID);
        }
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << DELIM << std::setprecision(5) << std::fixed << get_value(parID);               
            param_ofstream << DELIM << std::setprecision(2) << std::fixed << 100.0 * get_accept_rate(parID);
        }
        param_ofstream << std::endl;
    }
    else
    {
        std::string err_text = "Output file closed unexpectedly.";
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// MCMC : print augm variables ///////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print all augmented variables of each individual */
void MCMC::print_augm_variables(std::string save_augm_pathtofile)
{
    if (save_augm_pathtofile != UNDEFINED)
    {
        std::ofstream augm_ofstream;
        augm_ofstream.open(save_augm_pathtofile.c_str());

        if (augm_ofstream.is_open())
        {
            display_decorator("OUTPUT::augmented database", '-', 1);
            display_stdout("Saving augmented database to: " + save_augm_pathtofile);

            // change format of augmented data
            format_augmented_data();

            //=====Clean up list of headers for output files=====//
            data_headers.erase(find(data_headers.begin(), data_headers.end(), "index"));
            data_headers.erase(find(data_headers.begin(), data_headers.end(), "status"));
            for (std::string augmID : augm_full_set)
            {
                auto it = find(data_headers.begin(), data_headers.end(), augmID);
                if (it != data_headers.end())
                {
                    data_headers.erase(it);
                }
            }

            //=====Print column headers=====//
            augm_ofstream << "index" << DELIM << "status";
            for (std::string dataID : data_headers)
            {
                augm_ofstream << DELIM << dataID;
            }
            for (std::string augmID : augm_full_set)
            {
                augm_ofstream << DELIM << augmID;
            }
            augm_ofstream << std::endl;

            for (int i = 0; i < num_indiv; i++)
            {
                augm_ofstream << indiv[i].get_index() << DELIM << indiv[i].get_status();
                for (std::string dataID : data_headers)
                {
                    augm_ofstream << DELIM << indiv[i].get_data(dataID);
                }
                for (std::string augmID : augm_full_set)
                {
                    augm_ofstream << DELIM << indiv[i].get_augm(augmID);
                }
                augm_ofstream << std::endl;
            }

            augm_ofstream.close();
        }
        else
        {
            std::string err_text = "Unable to open file `" + save_augm_pathtofile + "`";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
    }
}