// build.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Unknown_Status_Warning = false; // Warning for unknown status

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MCMC : build  ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::build(bool is_infer, bool is_simulate,
                 std::string params_file, std::string data_file,
                 std::string augm_file, std::string setup_file,
                 std::string extra_file)
{

    //=====OPTION VALIDITY=====//

    // UPDATE : complete check for options validity

    if (is_infer && is_simulate)
    {
        raise_error_message("Incompatible commands: choose one between -infer and -simulate.");
    }

    // 1) Load Parameters and Check (required)
    if (params_file == UNDEFINED)
    {
        raise_error_message("Required argument: parameter path/to/file not found.");
    }
    else
    {
        build_parameters(params_file);
    }

    // 2) Load Data and Check (required)
    if (data_file == UNDEFINED)
    {
        raise_error_message("Required argument: database path/to/file not found.");
    }
    else
    {
        build_data(data_file);
    }

    // 3) Load Augmented Variables (optional)
    build_augmented(augm_file);

    // 4) Load Setup (optional)
    build_setup(setup_file, is_infer);

    // 5) Load Extra File (optional, user defined)
    // UPDATE: there could be plenty extra files to be store somewhere
    build_extra(extra_file);
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : build parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read model parameter file and build paramerter dictionary. */
void MCMC::build_parameters(std::string pathtofile)
{
    std::ifstream input(pathtofile.c_str());

    if (input.is_open())
    {
        display_decorator("Parameters");
        display_stdout("import parameters definition from:\n" + pathtofile);

        //=====Load column headers=====//
        std::vector<std::string> headers = read_column_headers(input);
        // check the validity of the headers
        for (std::string req : param_file_required_headers)
        {
            if (std::find(headers.begin(), headers.end(), req) == headers.end())
            {
                raise_error_message("Missing column `" + req + "`. Please check file: " + pathtofile);
            }
        }

        //=====Load parameters from the row with options=====//
        for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row); this_row.clear())
        {
            // add entry to shared parameter dictionary
            if (stobool(this_row.at("infer")))
            {
                assign_parameter(
                    this_row["name"],
                    std::stod(this_row["value"]),
                    this_row["prior"],
                    parse_dict_args(this_row["prior_args"]),
                    this_row["strategy"],
                    parse_dict_args(this_row["strategy_args"]),
                    parse_dict_args(this_row["range"]));
            }
            else
            {
                assign_parameter(this_row.at("name"), stod(this_row.at("value")));
            }
        }

        //=====display a summary of the used parameters=====//
        display_decorator("Parameters::definition", 1, '-');
        display_parameters();
    }
    else
    {
        raise_error_message("Impossible to read file: " + pathtofile);
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : build data  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read data file and pre-process the data. */
void MCMC::build_data(std::string pathtofile)
{
    std::ifstream input(pathtofile.c_str());
    std::vector<std::string> required_headers = {"index", "status", "house_index"};

    if (input.is_open())
    {
        display_decorator("Database");
        display_stdout("import dataset from:\n" + pathtofile);

        //=====Load column headers=====//
        data_headers = read_column_headers(input);
        // check the validity of the headers
        for (std::string req : required_headers)
        {
            if (std::find(data_headers.begin(), data_headers.end(), req) == data_headers.end())
            {
                raise_error_message("Missing column `" + req + "`. Please check file: " + pathtofile);
            }
        }

        //=====Load individuals from the row with values=====//
        std::vector<Individual> aux_individuals = {};
        for (ssargs this_row = {}; read_row_by_headers(input, data_headers, this_row); this_row.clear())
        {
            Individual this_indiv(row_to_sdargs(this_row));
            assign_rate_per_group_to_individual(this_indiv);
            set_rate_per_group(this_indiv);
            aux_individuals.push_back(this_indiv);
            report_data["n.indiv"]++;

            //=====check for the validity of the status=====//
            if (this_indiv.get_status() == INFECTIOUS)
            {
                report_data["n.cases"]++;
            }
            else if (this_indiv.get_status() != SUSCEPTIBLE)
            {
                report_data["n.unk.stat"]++;

                std::ostringstream oss;
                oss << "Status `" << this_indiv.get_status() << "` interpreted as missing."
                    << "\nUse " << SUSCEPTIBLE << "/" << INFECTIOUS << " for susceptible/infectious states."
                    << "\nIf it was intentional, ignore this message.";
                raise_warning_message(oss, Unknown_Status_Warning);
            }
        }

        //=====safe load households=====//
        int hh_idx;
        std::unordered_map<int, int> aux_house_idx = {};
        for (Individual this_indiv : aux_individuals)
        {
            hh_idx = this_indiv.get_data("house_index");

            // if household index doesn't exist, add new household
            if (aux_house_idx.find(hh_idx) == aux_house_idx.end())
            {
                aux_house_idx.emplace(hh_idx, num_househ);
                househ.emplace_back();
                househ[num_househ].set_index(hh_idx);
                num_househ++;
            }
            househ[aux_house_idx.at(hh_idx)].assign_individual(this_indiv);
        }
        for (int h = 0; h < num_househ; h++)
        {
            set_rate_per_contact(h);
        }

        //=====Produce Household List and Set Contact Rates=====//

        //=====Display Database and Model Reports=====//
        display_report();
    }
    else
    {
        raise_error_message("Impossible to open file: " + pathtofile);
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : set rate per group  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update model report and set group rates associate to new individual. */
void MCMC::set_rate_per_group(Individual &this_indiv)
{
    for (const auto &which : this_indiv.which_rate)
    {
        std::string rateID = which.first;
        if (rateID == CONTACT_RATE)
        {
            std::ostringstream oss;
            oss << "The group rate name `" << CONTACT_RATE << "` is reserved."
                << "Please choose a different one.";
            raise_error_message(oss);
        }

        std::string parID = which.second;

        //=====initialize model report=====//
        if (report_data["n.indiv"] == 0)
        {
            // if first individual, initialize rateID-report
            report_model.emplace(rateID, std::map<std::string, int>());
            report_model.at(rateID).emplace("[reference]", 0);
        }
        else
        {
            // all individuals must have the same set of rates of first individual
            if (report_model.find(rateID) == report_model.end())
            {
                raise_error_message("Group rate `" + rateID + "` has not been provided for all individuals.");
            }
        }

        //=====set rates and update report=====//
        if (parID == UNDEFINED)
        {
            // error : undefined rate
            std::ostringstream oss;
            oss << "Rate `" << rateID << "` is not provided for individual " << this_indiv.get_index() << ".";
            raise_error_message(oss);
        }
        else if (parID == REFERENCE)
        {
            // rate set to reference
            this_indiv.set_rate(rateID, get_reference());
            report_model.at(rateID).at("[reference]") += 1;
        }
        else
        {
            // rate associated to parameter
            this_indiv.set_rate(rateID, get_address(parID));
            if (report_model.at(rateID).find(parID) == report_model.at(rateID).end())
            {
                report_model.at(rateID).emplace(parID, 1);
            }
            else
            {
                report_model.at(rateID).at(parID) += 1;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : set rate per contact  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set contact rate according to user assignement. */
void MCMC::set_rate_per_contact(int h)
{
    std::string rateID = CONTACT_RATE;

    // allocate contact patterns
    int s = househ[h].house_size;
    househ[h].which_rate.resize(s, std::vector<std::string>(s, UNDEFINED));
    househ[h].contact_rates.resize(s, std::vector<const double *>(s, &REFERENCE_VALUE));

    if (h == 0)
    {
        // if first household, initialize contact rate report
        report_model.emplace(rateID, std::map<std::string, int>());
        report_model.at(rateID).emplace("[reference]", 0);
    }
    else
    {
        // all household must have a contact rate
        if (report_model.find(rateID) == report_model.end())
        {
            raise_error_message("Contact rate has not been provided for all households.");
        }
    }

    for (int i = 0; i < househ[h].house_size; i++)
    {
        for (int k = 0; k < househ[h].house_size; k++)
        {
            // note : from `i` --> to `k`
            if (i != k)
            {
                assign_rate_per_contact_to_household(h, i, k);
                std::string parID = househ[h].which_rate[i][k];

                if (parID == UNDEFINED)
                {
                    std::ostringstream oss;
                    oss << "Contact rate not provided for individuals " << househ[h].indiv[i].get_index()
                        << " and " << househ[h].indiv[k].get_index() << " within household " << househ[h].get_index() << ".";
                    raise_error_message(oss);
                }
                else if (parID == REFERENCE)
                {
                    househ[h].set_rate(i, k, get_reference()); // reference rate
                    report_model.at(rateID).at("[reference]") += 1;
                }
                else
                {
                    househ[h].set_rate(i, k, get_address(parID)); // parameter rate
                    if (report_model.at(rateID).find(parID) == report_model.at(rateID).end())
                    {
                        report_model.at(rateID).emplace(parID, 1);
                    }
                    else
                    {
                        report_model.at(rateID).at(parID) += 1;
                    }
                }

                // count the total number of pairs
                report_data.at("_n.contacts_")++;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : display report  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void MCMC::display_report()
{
    if (VERBOSE)
    {
        //=====Perform Database Analysis=====//
        for (int h = 0; h < num_househ; h++)
        {
            report_data["min.HH.size"] = std::min(househ[h].house_size, report_data["min.HH.size"]);
            report_data["max.HH.size"] = std::max(househ[h].house_size, report_data["max.HH.size"]);
            int aux = househ[h].infectious_list.size();
            report_data["no.HH.case"] += (aux == 0);
            report_data["lst.1.case"] += (aux > 0);
            report_data["lst.2.case"] += (aux > 1);
        }

        //=================================//
        //=====Display Database Report=====//
        //=================================//

        display_decorator("Database::summary", 1, '-');
        std::ostringstream oss;
        oss << "n. of individuals: " << report_data["n.indiv"]
            << "\n    * n. of cases: " << report_data["n.cases"]
            << std::fixed << std::setprecision(2)
            << " (" << (100.0 * report_data["n.cases"] / report_data["n.indiv"]) << " %)";
        if (report_data["n.unk.stat"] > 0)
        {
            oss << "\nn. of missing status : " << report_data["n.unk.stat"];
        }
        oss << "\nn. of households: " << num_househ;
        if (report_data["no.HH.case"] > 0)
        {
            oss << "\n    * with no cases: " << report_data["no.HH.case"]
                << std::fixed << std::setprecision(2)
                << " (" << (100.0 * report_data["no.HH.case"] / num_househ) << " %)";
        }
        oss << "\n    * with at least 1 case: " << report_data["lst.1.case"]
            << std::fixed << std::setprecision(2)
            << " (" << (100.0 * report_data["lst.1.case"] / num_househ) << " %)"
            << "\n    * with at least 2 cases: " << report_data["lst.2.case"]
            << std::fixed << std::setprecision(2)
            << " (" << (100.0 * report_data["lst.2.case"] / num_househ) << " %)"
            << "\nmin n. house members: " << report_data["min.HH.size"]
            << "\nmax n. house members: " << report_data["max.HH.size"];
        display_stdout(oss);
        oss.str("");
        oss.clear();

        //==============================//
        //=====Display Model Report=====//
        //==============================//

        // check if user defined groups other than reference
        bool has_non_reference_groups = false;
        if (!report_model.empty())
        {
            for (const auto &entry : report_model)
            {
                if (entry.second.size() > 1)
                {
                    has_non_reference_groups = true;
                    break;
                }
            }
        }

        // in such case display the model report
        if (has_non_reference_groups)
        {
            std::ostringstream oss;
            display_decorator("Database::rates", 1, '-');

            for (const auto &entry : report_model)
            {
                if (entry.first != CONTACT_RATE)
                {
                    oss << entry.first << ":";
                    for (auto const &x : entry.second)
                    {
                        oss << "\n    * " << x.first << ": " << x.second
                            << std::fixed << std::setprecision(2)
                            << " (" << (100.0 * x.second / report_data["n.indiv"]) << " %)";
                    }
                    oss << "\n";
                }
            }
            oss << CONTACT_RATE << ":";
            for (auto const &x : report_model.at(CONTACT_RATE))
            {
                oss << "\n    * " << x.first << ": " << x.second
                    << std::fixed << std::setprecision(2)
                    << " (" << (100.0 * x.second / report_data["_n.contacts_"]) << " %)";
            }
            display_stdout(oss);
        }
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : build augmented  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file */
void MCMC::build_augmented(std::string pathtofile)
{
    if (pathtofile != UNDEFINED)
    {
        std::ifstream input(pathtofile.c_str());
        std::vector<std::string> required_headers = {
            "name", "value", "model", "model_args", "strategy", "strategy_args", "move_from", "range", "init_rand"};

        if (input.is_open())
        {
            display_decorator("Augmented");
            display_stdout("import augmented variables definition from:\n" + pathtofile);

            //=====Load column headers=====//
            std::vector<std::string> headers = read_column_headers(input);
            // check the validity of the headers
            for (std::string req : required_headers)
            {
                if (std::find(headers.begin(), headers.end(), req) == headers.end())
                {
                    raise_error_message("Missing column `" + req + "`. Please check file: " + pathtofile);
                }
            }

            //=====Load augmented from the row with options=====//
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row); this_row.clear())
            {
                load_augm_definition.push_back(this_row);
            }

            //=====assign augmented variables to all individuals=====//
            assign_augmented_variables();
        }
        else
        {
            raise_error_message("Impossible to read file: " + pathtofile);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : assign augmented variables  /////////////////////////////////
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
        }

        // list of augmented variables with a model
        // note: no contribution of constant variables to Metropolis
        if (loader.at("model") == UNDEFINED)
        {
            if (stobool(loader.at("init_rand")) == true)
            {
                raise_error_message("Unknown random initialization without a model for `" + augmID + "`.");
            }
        }
        else
        {
            augm_with_model_set.emplace(augmID);

            if (stobool(loader.at("init_rand")) == true)
            {
                augm_rand_init_set.emplace(augmID);
            }
        }
    }
    setup_data_augm = !augm_to_be_inferred.empty();

    // assign augmented variables to each individual
    for (int h = 0; h < num_househ; h++)
    {
        for (int i = 0; i < househ[h].house_size; i++)
        {
            for (ssargs loader : load_augm_definition)
            {
                Augmented augm_var(
                    type_to_double(loader.at("value"), h, i),
                    loader.at("model"),
                    typedict_to_sdargs_ext(loader.at("model_args")),
                    loader.at("strategy"),
                    typedict_to_sdargs(loader.at("strategy_args"), h, i),
                    type_to_double(loader.at("move_from"), h, i),
                    typedict_to_sdargs(loader.at("range"), h, i));
                househ[h].indiv[i].assign_augm(loader.at("name"), augm_var);
            }
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
        for (int h = 0; h < num_househ; h++)
        {
            for (int i = 0; i < househ[h].house_size; i++)
            {
                househ[h].indiv[i].init_augm(augmID, get_prng());
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : build setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Setup file */
void MCMC::build_setup(std::string pathtofile, bool is_infer)
{
    int num_threads = 1;

    if (pathtofile != UNDEFINED)
    {
        std::ifstream input(pathtofile.c_str());
        ssargs safe_load = {};

        if (input.is_open())
        {
            display_decorator("Setup");
            display_stdout("import setup definition from:\n" + pathtofile);

            // Create a std::stringstream of the current line
            for (std::string line; getline(input, line);)
            {
                std::string name, choice;
                std::istringstream iss(line);

                //=====NAME=====//
                if (!getline(iss, name, C_DELIM))
                {
                    raise_error_message("Unable to extract setup name.\nCheck file line: " + line);
                }

                //=====CHOICE=====//
                if (!getline(iss, choice))
                {
                    raise_error_message("Unable to extract setup choice.\nCheck file line: " + line);
                }

                // add the new setup
                safe_load.emplace(string_cleaner(name), string_cleaner(choice));
            }

            //-----verbose-----//
            if (safe_load.find("verbose") != safe_load.end())
            {
                VERBOSE = stobool(safe_load.at("verbose"));
            }

            //-----number of threads-----//
            if (safe_load.find("num_threads") != safe_load.end())
            {
                num_threads = std::stoi(safe_load.at("num_threads"));
            }

            //-----prng seed-----//
            if ((safe_load.find("seed") == safe_load.end()) || (safe_load.at("seed") == UNDEFINED))
            {
                auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();
                setup_seed = static_cast<long long>(time_now);
            }
            else
            {
                setup_seed = std::stoi(safe_load.at("seed"));
                if (setup_seed < 1)
                {
                    raise_error_message("Setup value of `seed` must be " + UNDEFINED + " (default) or greater than 0.");
                }
            }

            if (is_infer)
            {
                //-----number of iterations-----//
                if (safe_load.find("num_iter") != safe_load.end())
                {
                    setup_num_iter = std::stoi(safe_load.at("num_iter"));
                    if (setup_num_iter < 1)
                    {
                        raise_error_message("Setup value of `num_iter` must be greater than 1.");
                    }
                }

                //-----number of finite adaptive iterations-----//
                if (safe_load.find("num_adapt_iter") != safe_load.end())
                {
                    setup_num_adapt_iter = std::stoi(safe_load.at("num_adapt_iter"));
                    if (setup_num_adapt_iter < 0)
                    {
                        raise_error_message("Setup value of `num_adapt_iter` must be greater than 0.");
                    }
                }

                //-----thinning-----//
                if (safe_load.find("thinning") != safe_load.end())
                {
                    setup_thinning = std::stoi(safe_load.at("thinning"));
                    if (setup_thinning < 1)
                    {
                        raise_error_message("Setup value of `thinning` must be greater than 1.");
                    }
                }

                //-----random init parameters-----//
                if (safe_load.find("rnd_init_params") != safe_load.end())
                {
                    setup_rnd_init_params = stobool(safe_load.at("rnd_init_params"));
                }

                //-----inverse temperature-----//
                if (safe_load.find("inv_tmprt") != safe_load.end())
                {
                    setup_inv_tmprt = std::stod(safe_load.at("inv_tmprt"));
                    if ((setup_inv_tmprt < 0.) || (setup_inv_tmprt > 1.))
                    {
                        raise_error_message("Inverse temperature " + std::to_string(setup_inv_tmprt) + " out of range.");
                    }
                }
            }
        }
        else
        {
            raise_error_message("Impossible to read file: " + pathtofile);
        }
    }

    //=====Define Number of Threads (Monte_Carlo_Core)=====//
    set_num_threads(num_threads);

    //=====Pseudo-Random Number Generator Seeding=====//
    set_prng_seed(setup_seed);

    //=====Define MCMC Number of Prints=====//
    setup_num_prints = setup_num_iter / setup_thinning;

    //=====Display=====//
    display_decorator("Setup::definitions", 1, '-');
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nnum_threads" << C_DELIM << get_num_threads()
        << "\nverbose" << C_DELIM << ((VERBOSE) ? "true" : "false");
    if (is_infer == true)
    {
        oss << "\nnum_iter" << C_DELIM << setup_num_iter
            << "\nnum_adapt_iter" << C_DELIM << setup_num_adapt_iter
            << "\nthinning" << C_DELIM << setup_thinning
            << "\nrnd_init_params" << C_DELIM << ((setup_rnd_init_params) ? "true" : "false");
        if (setup_inv_tmprt != 1.)
        {
            oss << "\ninv_tmprt" << C_DELIM << setup_inv_tmprt;
        }
    }

    display_stdout(oss);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : build extra  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Extra File (default unused) */
void MCMC::build_extra(std::string pathtofile)
{
    UNUSED(pathtofile);
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : init param file  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void MCMC::init_param_file(std::string output_filename)
{
    param_ofstream.open(output_filename.c_str());

    if (param_ofstream.is_open())
    {
        display_stdout("save MCMC iterations to:\n" + output_filename);

        param_ofstream << "iter" << C_DELIM << "loglike";
        for (std::string augmID : augm_to_be_inferred)
        {
            param_ofstream << C_DELIM << "augm_" << augmID << "__a%";
        }
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << C_DELIM << parID << C_DELIM << parID << "__a%";
        }
        param_ofstream << std::endl;
    }
    else
    {
        raise_error_message("Unable to open file: " + output_filename);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the model parameters at the given iteration time */
void MCMC::print_parameters(int step)
{
    if (param_ofstream.is_open())
    {
        param_ofstream << step << C_DELIM << std::setprecision(2) << std::fixed << get_global_loglike();
        for (std::string augmID : augm_to_be_inferred)
        {
            param_ofstream << C_DELIM << std::setprecision(2) << std::fixed << 100.0 * get_augm_accept_rate(augmID);
        }
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << C_DELIM << std::setprecision(5) << std::fixed << get_value(parID);
            param_ofstream << C_DELIM << std::setprecision(2) << std::fixed << 100.0 * get_accept_rate(parID);
        }
        param_ofstream << std::endl;
    }
    else
    {
        raise_error_message("Output file closed unexpectedly.");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : print augm variables  ///////////////////////////////////
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
            display_decorator("Output::augmented database", 1, '-');
            display_stdout("save augmented database to:\n" + save_augm_pathtofile);

            // change format of augmented data
            format_augmented_data();

            //=====Clean up list of headers for output files=====//
            data_headers.erase(find(data_headers.begin(), data_headers.end(), "index"));
            data_headers.erase(find(data_headers.begin(), data_headers.end(), "status"));
            data_headers.erase(find(data_headers.begin(), data_headers.end(), "house_index"));
            for (std::string augmID : augm_full_set)
            {
                auto it = find(data_headers.begin(), data_headers.end(), augmID);
                if (it != data_headers.end())
                {
                    data_headers.erase(it);
                }
            }

            //=====Print column headers=====//
            augm_ofstream << "index" << C_DELIM << "status" << C_DELIM << "house_index";
            for (std::string dataID : data_headers)
            {
                augm_ofstream << C_DELIM << dataID;
            }
            for (std::string augmID : augm_full_set)
            {
                augm_ofstream << C_DELIM << augmID;
            }
            augm_ofstream << std::endl;

            //=====Print augmented variables for each individual=====//
            for (int h = 0; h < num_househ; h++)
            {
                for (int i = 0; i < househ[h].house_size; i++)
                {
                    augm_ofstream << househ[h].indiv[i].get_index()
                                  << C_DELIM << househ[h].indiv[i].get_status()
                                  << C_DELIM << househ[h].get_index();
                    for (std::string dataID : data_headers)
                    {
                        augm_ofstream << C_DELIM << househ[h].indiv[i].get_data(dataID);
                    }
                    for (std::string augmID : augm_full_set)
                    {
                        augm_ofstream << C_DELIM << househ[h].indiv[i].get_augm(augmID);
                    }
                    augm_ofstream << std::endl;
                }
            }

            augm_ofstream.close();
        }
        else
        {
            raise_error_message("Unable to open file: " + save_augm_pathtofile);
        }
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//