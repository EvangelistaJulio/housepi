// build.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Unknown_Status_Warning = false;       // unknown status
bool Fix_Init_Mult_Chains_Warning = false; // fixed paramter values when multiple chains
bool Unkown_Setup_Entry_Warning = false;   // wrong entries in setup file

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MCMC : build  ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::build(std::string setup_file, std::string params_file,
                 std::string data_file, std::string augm_file,
                 std::string extra_file, bool verbose)
{
    //===== reset internal status =====//
    intern_is_infer_complete = false;

    // 1) Load Setup (optional)
    build_setup(setup_file, verbose);

    // 2) Load Parameters and Check (required)
    build_parameters(params_file, verbose);

    // 3) Load Data and Check (required)
    build_data(data_file, verbose);

    // 4) Load Augmented Variables (optional)
    build_augmented(augm_file, verbose);

    // 5) Load Extra File (optional, user defined)
    // UPDATE: there could be plenty extra files to be store somewhere
    build_extra(extra_file, verbose);
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
void MCMC::build_parameters(std::string pathtofile, bool verbose)
{
    if (pathtofile == UNDEFINED)
    {
        // error if nothing was previously defined
        if (param_to_be_static.empty() && param_to_be_inferred.empty())
        {
            screen::raise_error("Undefined parameter path/to/file.");
        }
    }
    else
    {
        // open input file stream
        std::ifstream input(pathtofile.c_str());

        if (!input.is_open())
        {
            screen::raise_error("Impossible to read file: " + pathtofile);
        }
        else
        {
            screen::display_decorator("Parameters", screen::C_DECOR, "import: " + pathtofile, verbose);

            //===== Load column headers =====//
            std::vector<std::string> headers = read_column_headers(input, param_file_headers, pathtofile);

            //===== Clean original list of parameters =====//
            drop_all_params();

            //===== Load parameters from the row with options =====//
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row); this_row.clear())
            {
                // add entry to shared parameter dictionary
                if (stobool(this_row.at("infer")))
                {
                    assign_parameter(
                        this_row["name"],
                        std::stod(this_row["value"]),
                        this_row["prior"],
                        parse_dict_args(this_row["pargs"]),
                        this_row["strategy"],
                        parse_dict_args(this_row["sargs"]),
                        parse_dict_args(this_row["range"]));
                }
                else
                {
                    assign_parameter(this_row.at("name"), stod(this_row.at("value")));
                }
            }

            //===== Display a summary of the used parameters =====//
            if (verbose)
            {
                screen::display_decorator("Parameters::definitions", '-');
                show_parameters();
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
////////////////////////////////////  MCMC : build data  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Use the loaded data to fill-in the respective report. */
void MCMC::update_report_data()
{
    report_data = {{"n.cases", 0},
                   {"n.unk.stat", 0},
                   {"no.HH.case", 0},
                   {"lst.1.case", 0},
                   {"lst.2.case", 0},
                   {"max.n.HH.cases", 0},
                   {"n.HH.with.max.cases", 0},
                   {"min.HH.size", 1000},
                   {"n.HH.with.min.size", 0},
                   {"max.HH.size", 0},
                   {"n.HH.with.max.size", 0}};

    //===== Set Contact Rates and Perform Database Analysis =====//
    for (int h = 0; h < num_househ; h++)
    {
        int num_infec = househ[h].infectious_list.size();

        report_data["n.cases"] += num_infec;
        report_data["n.unk.stat"] += househ[h].unknown_status_list.size();

        report_data["no.HH.case"] += (num_infec == 0);
        report_data["lst.1.case"] += (num_infec > 0);
        report_data["lst.2.case"] += (num_infec > 1);

        // max number of house cases
        if (report_data["max.n.HH.cases"] < num_infec)
        {
            report_data["n.HH.with.max.cases"] = 1;
            report_data["max.n.HH.cases"] = num_infec;
        }
        else if (report_data["max.n.HH.cases"] == num_infec)
        {
            report_data["n.HH.with.max.cases"] += 1;
        }

        // minimum house number of members
        if (report_data["min.HH.size"] > househ[h].house_size)
        {
            report_data["n.HH.with.min.size"] = 1;
            report_data["min.HH.size"] = househ[h].house_size;
        }
        else if (report_data["min.HH.size"] == househ[h].house_size)
        {
            report_data["n.HH.with.min.size"] += 1;
        }

        // maximum house number of members
        if (report_data["max.HH.size"] < househ[h].house_size)
        {
            report_data["n.HH.with.max.size"] = 1;
            report_data["max.HH.size"] = househ[h].house_size;
        }
        else if (report_data["max.HH.size"] == househ[h].house_size)
        {
            report_data["n.HH.with.max.size"] += 1;
        }
    }
}

/* Read data file and pre-process the data. */
void MCMC::build_data(std::string pathtofile, bool verbose)
{
    if (pathtofile == UNDEFINED)
    {
        // error if nothing was previously defined
        if (househ.empty())
        {
            screen::raise_error("undefined database path/to/file.");
        }
    }
    else
    {
        // open input file stream
        std::ifstream input(pathtofile.c_str());
        std::vector<std::string> required_headers = {"index", "house_index"};

        if (!input.is_open())
        {
            screen::raise_error("Impossible to open file: " + pathtofile);
        }
        else
        {
            screen::display_decorator("Database", screen::C_DECOR, "import: " + pathtofile, verbose);

            //=====Load column headers=====//
            intern_data_headers = read_column_headers(input, required_headers, pathtofile);

            //===== Drop all previous data =====//
            househ.clear();
            num_indiv = 0;
            num_househ = 0;
            report_group_classes.clear();

            //=====Load individuals from the row with values=====//
            std::vector<Individual> aux_individuals = {};
            for (ssargs this_row = {}; read_row_by_headers(input, intern_data_headers, this_row); this_row.clear())
            {

                Individual this_indiv(row_to_sdargs(this_row));
                assign_rate_per_group_to_individual(this_indiv);
                set_rate_per_group(this_indiv);
                aux_individuals.push_back(this_indiv);
                num_indiv++;

                //=====check for the validity of the status=====//
                std::set<int> valid_statuses = {INFECTIOUS, MISSING_STATUS, SUSCEPTIBLE};
                int this_status = this_indiv.get_status();
                if(valid_statuses.find(this_status) == valid_statuses.end())
                {
                    this_indiv.set_status(MISSING_STATUS);

                    std::ostringstream oss;
                    oss << "Status `" << this_indiv.get_status() << "` interpreted as missing."
                        << "\nUse " << SUSCEPTIBLE << "/" << INFECTIOUS << " for susceptible/infectious states."
                        << "\nOtherwise use " << MISSING_STATUS << " for missing.";
                    screen::raise_warning(oss.str(), Unknown_Status_Warning);
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

            //===== Set Contact Rates =====//
            for (int h = 0; h < num_househ; h++)
            {
                set_rate_per_contact(h);
            }

            //===== Update Report for Data =====//
            update_report_data();

            //=====Display Database and Model Reports=====//
            if (verbose)
            {
                screen::display_decorator("Database::summary", '-');
                show_summary();
                show_group_classes();
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
            screen::raise_error(oss.str());
        }

        std::string parID = which.second;

        //=====initialize model report=====//
        if (num_indiv == 0)
        {
            // if first individual, initialize rateID-report
            report_group_classes.emplace(rateID, std::map<std::string, int>());
            report_group_classes.at(rateID).emplace("[reference]", 0);
        }
        else
        {
            // all individuals must have the same set of rates of first individual
            if (report_group_classes.find(rateID) == report_group_classes.end())
            {
                screen::raise_error("Group rate `" + rateID + "` has not been provided for all individuals.");
            }
        }

        //=====set rates and update report=====//
        if (parID == UNDEFINED)
        {
            // error : undefined rate
            std::ostringstream oss;
            oss << "Rate `" << rateID << "` is not provided for individual " << this_indiv.get_index() << ".";
            screen::raise_error(oss.str());
        }
        else if (parID == REFERENCE)
        {
            // rate set to reference
            this_indiv.set_rate(rateID, get_reference());
            report_group_classes.at(rateID).at("[reference]") += 1;
        }
        else if (parID != REFERENCE)
        {
            // rate associated to parameter
            assert_param(parID);
            this_indiv.set_rate(rateID, get_address(parID));
            if (report_group_classes.at(rateID).find(parID) == report_group_classes.at(rateID).end())
            {
                report_group_classes.at(rateID).emplace(parID, 1);
            }
            else
            {
                report_group_classes.at(rateID).at(parID) += 1;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : set rate per contact  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set contact rate according to user assignement.
FIXME : is it necessary occupy this space even if all REFERENCE?
*/
void MCMC::set_rate_per_contact(int h)
{
    // allocate contact patterns
    int s = househ[h].house_size;
    househ[h].which_rate.resize(s, std::vector<std::string>(s, UNDEFINED));
    househ[h].contact_rates.resize(s, std::vector<const double *>(s, &REFERENCE_VALUE));

    if (h == 0)
    {
        // if first household, initialize contact rate report
        report_group_classes.emplace(CONTACT_RATE, std::map<std::string, int>());
        report_group_classes.at(CONTACT_RATE).emplace("[reference]", 0);
    }
    else
    {
        // all household must have a contact rate
        if (report_group_classes.find(CONTACT_RATE) == report_group_classes.end())
        {
            screen::raise_error("Contact rate has not been provided for all households.");
        }
    }

    for (int k = 0; k < househ[h].house_size; k++)
    {
        for (int i = 0; i < househ[h].house_size; i++)
        {
            // note : from `k` --> to `i`
            if (k != i)
            {
                assign_rate_per_contact_to_household(h, k, i);
                std::string parID = househ[h].which_rate[k][i];

                if (parID == UNDEFINED)
                {
                    std::ostringstream oss;
                    oss << "Contact rate not provided for individuals " << househ[h].indiv[k].get_index()
                        << " and " << househ[h].indiv[i].get_index() << " within household " << househ[h].get_index() << ".";
                    screen::raise_error(oss.str());
                }
                else if (parID == REFERENCE)
                {
                    househ[h].set_rate(k, i, get_reference()); // reference rate
                    report_group_classes.at(CONTACT_RATE).at("[reference]") += 1;
                }
                else if (parID != EMPTY_RATE)
                {
                    assert_param(parID);
                    househ[h].set_rate(k, i, get_address(parID)); // parameter rate
                    if (report_group_classes.at(CONTACT_RATE).find(parID) == report_group_classes.at(CONTACT_RATE).end())
                    {
                        report_group_classes.at(CONTACT_RATE).emplace(parID, 1);
                    }
                    else
                    {
                        report_group_classes.at(CONTACT_RATE).at(parID) += 1;
                    }
                }

                // count the total number of pairs
                intern_n_contacts++;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : display summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the report of the database and of the model group rates. */
void MCMC::show_summary()
{
    using namespace screen;

    //===== PRINT =====//

    std::ostringstream oss;

    oss << "n. of individuals: " << num_indiv
        << "\n    * n. of cases: " << report_data["n.cases"]
        << " " << percent(report_data["n.cases"], num_indiv);
    if (report_data["n.unk.stat"] > 0)
    {
        oss << "\n    * n. of missing status : " << report_data["n.unk.stat"]
        << " " << percent(report_data["n.unk.stat"], num_indiv);
    }

    oss << "\nn. of households: " << num_househ;
    if (report_data["no.HH.case"] > 0)
    {
        oss << "\n    * with no cases: " << report_data["no.HH.case"]
            << " " << percent(report_data["no.HH.case"], num_househ);
    }

    oss << "\n    * with at least 1 case: " << report_data["lst.1.case"]
        << " " << percent(report_data["lst.1.case"], num_househ);
    if (report_data["lst.1.case"] > 0)
    {
        oss << "\n    * with at least 2 cases: " << report_data["lst.2.case"]
            << " " << percent(report_data["lst.2.case"], num_househ);
    }

    oss << "\nmax n. house cases: " << report_data["max.n.HH.cases"]
        << " " << percent(report_data["n.HH.with.max.cases"], num_househ)
        << "\nmin n. house members: " << report_data["min.HH.size"]
        << " " << percent(report_data["n.HH.with.min.size"], num_househ)
        << "\nmax n. house members: " << report_data["max.HH.size"]
        << " " << percent(report_data["n.HH.with.max.size"], num_househ);

    print_stdout(oss.str() + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : display report  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Display model report. */
void MCMC::show_group_classes()
{
    using namespace screen;

    std::ostringstream oss;
    for (const auto &entry : report_group_classes)
    {
        if ((entry.first != CONTACT_RATE) && (entry.first != EMPTY_RATE))
        {
            oss << entry.first << ":";
            for (auto const &x : entry.second)
            {
                oss << "\n    * " << x.first << ": " << x.second
                    << " " << percent(x.second, num_indiv);
            }
            oss << "\n";
        }
    }

    if (report_group_classes.at(CONTACT_RATE).find(EMPTY_RATE) != report_group_classes.at(CONTACT_RATE).end())
    {
        oss << CONTACT_RATE << ":";
        for (auto const &x : report_group_classes.at(CONTACT_RATE))
        {
            oss << "\n    * " << x.first << ": " << x.second
                << " " << percent(x.second, intern_n_contacts);
        }
        oss << "\n";
    }

    if (oss.str().size() > 0)
    {
        print_stdout(oss.str());
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
void MCMC::build_augmented(std::string pathtofile, bool verbose)
{
    if (pathtofile != UNDEFINED)
    {
        std::ifstream input(pathtofile.c_str());
        std::vector<std::string> required_headers = {
            "name", "value", "model", "model_args", "strategy", "strategy_args", "move_from", "range", "init_rand"};

        if (input.is_open())
        {
            screen::display_decorator("Augmented", screen::C_DECOR, "import: " + pathtofile, verbose);

            //=====Load column headers=====//
            std::vector<std::string> headers = read_column_headers(input, required_headers, pathtofile);

            //=====Load augmented from the row with options=====//
            load_augm_definition.clear();
            augm_full_set.clear();
            augm_to_be_inferred.clear();
            augm_with_model_set.clear();
            augm_rand_init_set.clear();
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row); this_row.clear())
            {
                load_augm_definition.push_back(this_row);
            }
        }
        else
        {
            screen::raise_error("Impossible to read file: " + pathtofile);
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
                screen::raise_error("Unknown random initialization without a model for `" + augmID + "`.");
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
    if (!load_augm_definition.empty())
    {
        //=====assign augmented variables to all individuals=====//
        assign_augmented_variables();

        //=====random initialization on demand=====//
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
void MCMC::build_setup(std::string pathtofile, bool verbose)
{
    // default number of threads
    unsigned num_threads = 1;

    if (pathtofile == UNDEFINED)
    {
        screen::display_decorator("Setup", screen::C_DECOR, UNDEFINED, verbose);
    }
    else
    {
        std::ifstream input(pathtofile.c_str());

        if (!input.is_open())
        {
            screen::raise_error("Impossible to read setup file: " + pathtofile);
        }
        else
        {
            screen::display_decorator("Setup", screen::C_DECOR, "import: " + pathtofile, verbose);

            // Create a std::stringstream of the current line
            ssargs safe_load = {};
            for (std::string line; getline(input, line);)
            {
                std::string name, choice;
                std::istringstream iss(line);

                //===== NAME =====//
                if (!getline(iss, name, screen::C_DELIM))
                {
                    screen::raise_error("Unable to extract setup entry.\nCheck file line: " + line);
                }

                //===== CHOICE =====//
                if (!getline(iss, choice))
                {
                    screen::raise_error("Unable to extract setup choice.\nCheck file line: " + line);
                }

                // add the new setup
                safe_load.emplace(string_cleaner(name), string_cleaner(choice));
            }

            //---- Check that all provided setups are correct -----//
            std::vector<std::string> file_entries = {
                "num_threads", "seed", "num_chains",
                "num_iter", "num_adapt", "thinning",
                "burn_in", "rnd_init_params", "inv_temper"};

            std::vector<std::string> wrong_entries = {};
            for (const auto &pair : safe_load)
            {
                if (find(file_entries.begin(), file_entries.end(), pair.first) == file_entries.end())
                {
                    wrong_entries.push_back(pair.first);
                }
            }
            if (wrong_entries.size() > 0)
            {
                std::ostringstream oss;
                oss << "Unknown setup option(s): `" << wrong_entries[0] << "`";
                for (size_t k = 1; k < wrong_entries.size(); k++)
                {
                    oss << ", `" << wrong_entries[k] << "`";
                }
                oss << ". Available setup options are: `" << file_entries[0] << "`";
                for (size_t k = 1; k < file_entries.size(); k++)
                {
                    oss << ", `" << file_entries[k] << "`";
                }
                screen::raise_warning(oss.str() + ".", Unkown_Setup_Entry_Warning);
            }

            //===== NUMBER OF THREADS =====//
            if (safe_load.find("num_threads") != safe_load.end())
            {
                num_threads = std::stoi(safe_load.at("num_threads"));
            }

            //===== PRNG SEED =====//
            if ((safe_load.find("seed") != safe_load.end()) && (safe_load.at("seed") != UNDEFINED))
            {
                setup_seed = std::stoull(safe_load.at("seed"));
            }

            //===== NUMBER OF CHAINS =====//
            if (safe_load.find("num_chains") != safe_load.end())
            {
                setup_num_chains = std::stoul(safe_load.at("num_chains"));
                if (setup_num_iter < 1)
                {
                    screen::raise_error(
                        "number of chains must be greater than 0. Check your setup.");
                }
            }

            //===== NUMBER OF ITERATIONS =====//
            if (safe_load.find("num_iter") != safe_load.end())
            {
                setup_num_iter = std::stoul(safe_load.at("num_iter"));
                if (setup_num_iter < 2)
                {
                    screen::raise_error(
                        "number of iterations must be greater than 1. Check your setup.");
                }
            }

            //===== NUMBER OF (FINITE) ADAPTIVE ITERATIONS =====//
            if (safe_load.find("num_adapt") != safe_load.end())
            {
                setup_num_adapt = std::stoul(safe_load.at("num_adapt"));
            }

            //===== THINNING =====//
            if (safe_load.find("thinning") != safe_load.end())
            {
                setup_thinning = std::stoul(safe_load.at("thinning"));
                if (setup_thinning >= setup_num_iter)
                {
                    screen::raise_error(
                        "thinning must be smaller than number of iterations. Check your setup.");
                }
                else if (setup_thinning == 0)
                {
                    screen::raise_error(
                        "thinning must be greater than 0. Check your setup.");
                }
            }

            //===== BURN-IN =====//
            if (safe_load.find("burn_in") != safe_load.end())
            {
                setup_burn_in = std::stoul(safe_load.at("burn_in"));
                if (setup_burn_in >= setup_num_iter)
                {
                    screen::raise_error(
                        "burn-in must be smaller than number of iterations. Check your setup.");
                }
            }

            //===== RANDOM INITIALIZATION OF PARAMETERS =====//
            if (safe_load.find("rnd_init_params") != safe_load.end())
            {
                setup_rnd_init_params = stobool(safe_load.at("rnd_init_params"));
                if ((setup_num_chains > 1) && (setup_rnd_init_params == false))
                {
                    screen::raise_warning("Unable to fix initial parameters when running multiple chains.",
                                          Fix_Init_Mult_Chains_Warning);
                    setup_rnd_init_params = true;
                }
            }

            //===== INVERSE TEMPERATURE =====//
            if (safe_load.find("inv_temper") != safe_load.end())
            {
                setup_inv_temper = std::stod(safe_load.at("inv_temper"));
                if ((setup_inv_temper < 0.) || (setup_inv_temper > 1.))
                {
                    screen::raise_error(
                        "Inverse temperature " + std::to_string(setup_inv_temper) + " out of range [0, 1]. Check your setup.");
                }
            }
        }
    }

    //===== Define Number of Threads (Monte_Carlo_Core) =====//
    set_num_threads(num_threads);

    //===== Pseudo-Random Number Generator Seeding =====//
    set_prng_seed(setup_seed);

    //===== Display =====//
    if (verbose)
    {
        screen::display_decorator("Setup::definitions", '-');
        show_setup();
    }
}

/* Display the complete choices for the setup. */
void MCMC::show_setup()
{
    using namespace screen;

    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nnum_threads" << C_DELIM << get_num_threads()
        << "\nnum_chains" << C_DELIM << setup_num_chains
        << "\nnum_adapt" << C_DELIM << setup_num_adapt
        << "\nburn_in" << C_DELIM << setup_burn_in
        << "\nnum_iter" << C_DELIM << setup_num_iter
        << "\nthinning" << C_DELIM << setup_thinning
        << "\nrnd_init_params" << C_DELIM << truth(setup_rnd_init_params);
    if (setup_inv_temper != 1.)
    {
        oss << "\ninv_temper" << C_DELIM << setup_inv_temper;
    }
    oss << "\n";

    print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : build extra  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Extra File (default unused) */
void MCMC::build_extra(std::string pathtofile, bool verbose)
{
    UNUSED(pathtofile);
    UNUSED(verbose);
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
void MCMC::init_inferstate_file(std::string output_filename)
{
    using namespace screen;

    if (output_filename == UNDEFINED)
    {
        raise_error("Undefined output file.");
    }
    else
    {
        inferstate_ofstream.open(output_filename.c_str());

        if (inferstate_ofstream.is_open())
        {
#ifndef USE_RCPP
            print_stdout("save: " + output_filename + "\n");
#endif
            // print loglikelihood
            inferstate_ofstream << "loglike";

            /*
            // augmented variable acceptance rate
            for (std::string augmID : augm_to_be_inferred)
            {
                inferstate_ofstream << C_DELIM << "augm_" << augmID << "__a%";
            }
            */

            for (std::string parID : param_to_be_inferred)
            {
                inferstate_ofstream << C_DELIM << parID;
                // variable acceptance rate
                // << C_DELIM << parID << "__a%";
            }

            // which chain
            if (setup_num_chains > 1)
            {
                inferstate_ofstream << C_DELIM << "CHAIN__";
            }

            inferstate_ofstream << std::endl;
        }
        else
        {
            raise_error("Unable to open file: " + output_filename);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print the model parameters at the given iteration time */
void MCMC::print_inferstate()
{
    if (inferstate_ofstream.is_open())
    {
        // print loglikelihood
        inferstate_ofstream << std::setprecision(2) << std::fixed << get_global_loglike();

        /*
        for (std::string augmID : augm_to_be_inferred)
        {
            inferstate_ofstream << C_DELIM << std::setprecision(2) << std::fixed << 100.0 * get_augm_accept_rate(augmID);
        }
        */

        // print parameter values
        for (std::string parID : param_to_be_inferred)
        {
            inferstate_ofstream << screen::C_DELIM << std::setprecision(5) << std::fixed << get_value(parID);
            // inferstate_ofstream << C_DELIM << std::setprecision(2) << std::fixed << 100.0 * get_accept_rate(parID);
        }

        // print which chain
        if (setup_num_chains > 1)
        {
            inferstate_ofstream << screen::C_DELIM << intern_chain;
        }

        inferstate_ofstream << std::endl;
    }
    else
    {
        screen::raise_error("Output file closed unexpectedly.");
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
            // change format of augmented data
            format_augmented_data();

            //=====Clean up list of headers for output files=====//
            std::vector<std::string> headers = intern_data_headers;

            // remove default headers
            headers.erase(find(headers.begin(), headers.end(), "index"));
            headers.erase(find(headers.begin(), headers.end(), "status"));
            headers.erase(find(headers.begin(), headers.end(), "house_index"));

            // remove eventual variable with the same name of augmented variables
            for (std::string augmID : augm_full_set)
            {
                auto it = find(headers.begin(), headers.end(), augmID);
                if (it != headers.end())
                {
                    headers.erase(it);
                }
            }

            //=====Print column headers=====//
            // defaults
            augm_ofstream << "index" << screen::C_DELIM << "status" << screen::C_DELIM << "house_index";
            // data
            for (std::string dataID : headers)
            {
                augm_ofstream << screen::C_DELIM << dataID;
            }
            // augmented variables
            for (std::string augmID : augm_full_set)
            {
                augm_ofstream << screen::C_DELIM << augmID;
            }
            augm_ofstream << std::endl;

            //=====Print augmented variables for each individual=====//
            for (int h = 0; h < num_househ; h++)
            {
                for (int i = 0; i < househ[h].house_size; i++)
                {
                    augm_ofstream << househ[h].indiv[i].get_index()
                                  << screen::C_DELIM << househ[h].indiv[i].get_status()
                                  << screen::C_DELIM << househ[h].get_index();
                    for (std::string dataID : headers)
                    {
                        augm_ofstream << screen::C_DELIM << househ[h].indiv[i].get_data(dataID);
                    }
                    for (std::string augmID : augm_full_set)
                    {
                        augm_ofstream << screen::C_DELIM << househ[h].indiv[i].get_augm(augmID);
                    }
                    augm_ofstream << std::endl;
                }
            }

            augm_ofstream.close();
        }
        else
        {
            screen::raise_error("Unable to open file: " + save_augm_pathtofile);
        }
    }
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//