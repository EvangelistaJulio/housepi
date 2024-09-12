/*
 * build.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Param_Infer_Fixed_Move_Warning = false;
bool Non_Specified_Type_Warning = false;
bool Non_Specified_Value_Warning = false;
bool Unknown_Status_Warning = false;       // unknown status
bool Fix_Init_Mult_Chains_Warning = false; // fixed paramter values when multiple chains
bool Unkown_Setup_Entry_Warning = false;   // wrong entries in setup file
bool Unknown_Household_Data_Structure_Warning = false;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// declaration for type in data file
const std::string TYPE_PARAM = "P";   // Parameter type declaration
const std::string TYPE_DATA = "D";    // Data type declaration
const std::string TYPE_SCALAR = "S";  // Scalar type declaration
const std::string TYPE_COMMAND = "C"; // Command type declaration

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : methods (auxiliary functions)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Clean std::string from forbidden characters and convert to UNDEFINED. */
std::string MCMC::string_cleaner(std::string tkn)
{
    // forbidden characters
    for (char r : {'"', '\r', '\n'})
    {
        tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
    }

    // check for undefined
    if ((tkn == UNDEFINED) || (tkn.size() == 0))
    {
        tkn = UNDEFINED;
    }

    return tkn;
}

/* Convert string to boolean variable. */
bool MCMC::stobool(const std::string tkn)
{
    bool output = false;

    if (tkn == "TRUE" || tkn == "true" || tkn == "True" || tkn == "1")
    {
        output = true;
    }
    else if (tkn != "FALSE" && tkn != "false" && tkn != "False" && tkn != "0")
    {
        screen::raise_error("Unrecognized boolean value for `" + tkn + "`.");
    }

    return output;
}

/* Parse the list from the string. */
std::vector<std::string> parse_list(const std::string tkn, const char delim)
{
    std::istringstream args_ss(tkn);
    std::vector<std::string> output = {};
    std::string tmp;

    // rea the entries of the dictionary
    while (getline(args_ss, tmp, delim))
    {
        output.push_back(tmp);
    }

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Read column headers. */
std::vector<std::string> MCMC::read_column_headers(std::ifstream &input,
                                                   std::string pathtofile,
                                                   std::vector<std::string> required_headers)
{
    std::string line, tkn;
    getline(input, line);
    std::istringstream iss(line);
    std::vector<std::string> headers = {};

    while (getline(iss, tkn, screen::C_DELIM))
    {
        headers.push_back(string_cleaner(tkn));
    }

    // check that all required headers are present
    if (!required_headers.empty())
    {
        bool is_missing = false;
        std::ostringstream oss;

        // check the validity of the headers
        for (const std::string &req : required_headers)
        {
            if (std::find(headers.begin(), headers.end(), req) == headers.end())
            {
                if (!is_missing)
                {
                    // point to the faulty file
                    oss << "Check file: " << pathtofile << "\nMissing column(s): `" + req + "`";
                    is_missing = true;
                }
                else
                {
                    oss << ", `" << req << "`";
                }
            }
        }

        // raise the error if something was missing
        if (is_missing)
        {
            screen::raise_error(oss.str() + ".");
        }
    }

    return headers;
}

/* Load files from the row with headers */
bool MCMC::read_row_by_headers(std::ifstream &input, std::vector<std::string> headers, ssargs &this_row)
{
    std::string line, tkn;
    if (getline(input, line))
    {
        std::istringstream iss(line);
        int num_cols = 0;

        while (getline(iss, tkn, screen::C_DELIM))
        {
            this_row.emplace(headers[num_cols], string_cleaner(tkn));
            num_cols++;
        }
        return true;
    }
    else
    {
        return false;
    }
}

/* convert the output of a row to a string-to-double dictionary */
sdargs MCMC::row_to_sdargs(ssargs &this_row)
{
    sdargs row_with_double = {};
    for (auto const &x : this_row)
    {
        row_with_double.emplace(x.first, stod(x.second));
    }
    return row_with_double;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : build setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Setup file */
void MCMC::build_setup(const std::string pathtofile, bool verbose)
{
    // default number of threads
    unsigned num_threads = 1;

    if (pathtofile == UNDEFINED)
    {
        screen::print_decorator("Setup", screen::C_DECOR, UNDEFINED, verbose);
    }
    else
    {
        std::ifstream input(pathtofile.c_str());

        if (input.is_open())
        {
            screen::print_decorator("Setup", screen::C_DECOR, "import: " + pathtofile, verbose);

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
                "burn_in", "rnd_init_params",
                "num_rungs", "inv_temper"};

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
                setup_inv_T = std::stold(safe_load.at("inv_temper"));
                if ((setup_inv_T < 0.) || (setup_inv_T > 1.))
                {
                    screen::raise_error(
                        "Inverse temperature " + std::to_string(setup_inv_T) + " out of range [0, 1]. Check your setup.");
                }
            }
        }
        else
        {
            screen::raise_error("Impossible to open setup file `" + pathtofile + "`.");
        }
    }

    //===== Define Number of Threads (Monte_Carlo_Core) =====//
    set_num_threads(num_threads);

    //===== Pseudo-Random Number Generator Seeding =====//
    set_prng_seed(setup_seed);

    //===== Number of MCMC sweeps per single chain =====//
    // note : Euclidean divisions!
    __mcmc_sweeps = setup_num_adapt / setup_thinning;
    __mcmc_sweeps += setup_burn_in / setup_thinning;
    __mcmc_sweeps += (setup_num_iter - setup_burn_in) / setup_thinning;
    __mcmc_sweeps *= setup_thinning * setup_num_chains;

    //===== MCMC Sample Size =====//
    dgn_ss = (setup_num_iter - setup_burn_in) / setup_thinning;
    dgn_ss *= setup_num_chains;

    //===== Display =====//
    if (verbose)
    {
        screen::print_decorator("Setup::definitions", '-');
        show_setup();
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Parse the list from the string. */
std::vector<std::string> MCMC::parse_list_args(const std::string tkn, const char delim)
{
    std::istringstream args_ss(tkn);
    std::vector<std::string> output = {};
    std::string tmp;

    // rea the entries of the dictionary
    while (getline(args_ss, tmp, delim))
    {
        output.push_back(tmp);
    }

    return output;
}

/* Combine after control two lists into a dictionary */
sdargs MCMC::merge_list_args(const std::string tkn_keys, const std::string tkn_vals)
{
    std::vector<std::string> keys = parse_list_args(tkn_keys);
    std::vector<std::string> vals = parse_list_args(tkn_vals);
    sdargs output = {};

    if (keys.size() != vals.size())
    {
        // inform parameter about the error
        output.emplace(IO_MISMATCH, INF_DOUBLE);
    }
    else
    {
        for (size_t x = 0; x < keys.size(); x++)
        {
            output.emplace(keys[x], std::stod(vals[x]));
        }
    }
    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : build parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read model parameter file and build paramerter dictionary. */
void MCMC::build_parameters(const std::string pathtofile, bool verbose)
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

        if (input.is_open())
        {
            screen::print_decorator("Parameters", screen::C_DECOR, "import: " + pathtofile, verbose);

            //===== Load column headers =====//
            std::vector<std::string> headers = read_column_headers(input, pathtofile, param_file_headers);

            //===== Clean original list of parameters =====//
            drop_all_params();

            //===== Load parameters from the row with options =====//
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row); this_row.clear())
            {
                // add entry to shared parameter dictionary
                if (stobool(this_row.at("infer")))
                {
                    // load range
                    sdargs range = {};
                    if (this_row["min"] != UNDEFINED)
                    {
                        range.emplace("min", std::stod(this_row["min"]));
                    }
                    if (this_row["max"] != UNDEFINED)
                    {
                        range.emplace("max", std::stod(this_row["max"]));
                    }

                    // assign parameter to the set of variables
                    assign_parameter(
                        this_row["name"],
                        std::stod(this_row["value"]),
                        this_row["prior"],
                        merge_list_args(this_row["args"], this_row["argsval"]),
                        this_row["rule"],
                        merge_list_args(this_row["rargs"], this_row["rargsval"]),
                        range);
                }
                else
                {
                    assign_parameter(this_row.at("name"), stod(this_row.at("value")));
                }
            }

            //===== Show a summary of the used parameters =====//
            if (verbose)
            {
                screen::print_decorator("Parameters::definitions", '-');
                show_parameters();
            }
        }
        else
        {
            screen::raise_error("Impossible to open parameter file `" + pathtofile + "`.");
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : build data  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign group rate to an individual  */
void MCMC::assign_rate_per_group_to_individual(Individual &this_indiv)
{
    // this_indiv.assign_rate_per_group(EMPTY_RATE, EMPTY_RATE);
    UNUSED(this_indiv);
}

/* Assign the contact patterns between a couple of individuals. */
void MCMC::assign_rate_per_contact_to_household(int h, int from, int to)
{
    househ[h].assign_rate_per_contact(from, to, EMPTY_RATE);
}

/* Use the loaded data to fill-in the respective report. */
void MCMC::update_report_data()
{
    report_data = {{"n.cases", 0},
                   {"n.miss.stat", 0},
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
        int num_infec = househ[h].num_infectious;

        report_data["n.cases"] += num_infec;
        report_data["n.miss.stat"] += househ[h].num_missing_status;

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
void MCMC::build_data(const std::string pathtofile, bool verbose)
{
    // define valid status set
    std::set<int> valid_statuses = {INFECTIOUS, MISSING_STATUS, SUSCEPTIBLE};
    // define standard headers set
    std::set<std::string> standard_headers = {"house_index", "index", "status"};

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

        if (input.is_open())
        {
            screen::print_decorator("Database", screen::C_DECOR, "import: " + pathtofile, verbose);

            //===== Load Column Headers =====//
            __data_headers = read_column_headers(input, pathtofile);

            // boolean to assess the presence of indexes
            bool is_index = (std::find(__data_headers.begin(), __data_headers.end(), "index") != __data_headers.end());
            bool is_house_index = (std::find(__data_headers.begin(), __data_headers.end(), "house_index") != __data_headers.end());

            if (!is_house_index)
            {
                std::ostringstream oss;
                oss << "Impossible to associate a household structure: missing data column `house_index`."
                    << "\nIf it was intentional, ignore this message.";
                screen::raise_warning(oss.str(), Unknown_Household_Data_Structure_Warning);
            }

            //===== Drop all previous data =====//
            num_indiv = 0;
            num_househ = 0;
            househ.clear();
            report_group_classes.clear();

            //===== Load individuals from the row with values =====//

            std::vector<Individual> aux_individuals = {};
            for (ssargs this_row = {}; read_row_by_headers(input, __data_headers, this_row); this_row.clear())
            {
                sdargs args = row_to_sdargs(this_row);
                Individual this_indiv(args);
                // associate data index or default provided index
                this_indiv.set_index((is_index) ? args.at("index") : num_indiv);

                assign_rate_per_group_to_individual(this_indiv);
                set_rate_per_group(this_indiv);
                aux_individuals.push_back(this_indiv);
                num_indiv++;

                //===== Check for the validity of the status =====//
                if (valid_statuses.find(this_indiv.get_status()) == valid_statuses.end())
                {
                    this_indiv.set_status(MISSING_STATUS);

                    std::ostringstream oss;
                    oss << "Status `" << this_indiv.get_status() << "` of individual `" << this_indiv.get_index() << "` interpreted as missing."
                        << "\nUse " << SUSCEPTIBLE << "/" << INFECTIOUS << " for susceptible/infectious status."
                        << " Use " << MISSING_STATUS << " for missing status.";
                    screen::raise_warning(oss.str(), Unknown_Status_Warning);
                }
            }

            //===== Safe Load of Households =====//
            int hh_idx;
            std::unordered_map<int, int> aux_hh_idx = {};
            for (Individual this_indiv : aux_individuals)
            {
                // user provided structure or no structure
                hh_idx = (is_house_index) ? this_indiv.get_data("house_index") : num_househ;

                // if household index doesn't exist, add new household
                // note : this is redundant if no structure
                if (aux_hh_idx.find(hh_idx) == aux_hh_idx.end())
                {
                    househ.emplace_back();                  // add a new empty household
                    househ[num_househ].set_index(hh_idx);   // associate its index
                    aux_hh_idx.emplace(hh_idx, num_househ); // map index to household number
                    num_househ++;                           // increment household number
                }
                househ[aux_hh_idx.at(hh_idx)].assign_individual(this_indiv);
            }

            //===== Set Contact Rates =====//
            for (int h = 0; h < num_househ; h++)
            {
                set_rate_per_contact(h);
            }

            //===== Update Report for Data =====//
            update_report_data();

            //===== Display Database and Model Reports =====//
            if (verbose)
            {
                screen::print_decorator("Database::summary", '-');
                show_summary();
                show_group_classes();
            }

            //===== Clean up data entries for output files =====//
            // eventually remove standard headers to reorder them at the end
            for (const std::string &col : standard_headers)
            {
                __data_headers.erase(
                    std::remove(__data_headers.begin(), __data_headers.end(), col),
                    __data_headers.end());
            }
        }
        else
        {
            screen::raise_error("Impossible to open data file `" + pathtofile + "`.");
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

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
            this_indiv.set_rate(rateID, get_param_address(parID));
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

/* Set contact rate according to user assignement. */
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
                    househ[h].set_rate(k, i, get_param_address(parID)); // parameter rate
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
                __n_contacts++;
            }
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a typed string to a double */
double MCMC::type_to_double(const std::string val, const std::string typ, int h, int i)
{
    double output = INF_DOUBLE;

    if ((val == UNDEFINED) && (typ == UNDEFINED))
    {
        output = INF_DOUBLE; // default value
    }
    else if ((val != UNDEFINED) && (typ == UNDEFINED))
    {
        try
        {
            output = std::stod(val);
        }
        catch (const std::invalid_argument &e)
        {
            std::ostringstream oss;
            oss << "Value `" + val + "` can't be cast to scalar value. Provide the respective type.";
            screen::raise_warning(oss.str(), Non_Specified_Type_Warning);
        }
    }
    else if ((val == UNDEFINED) && (typ != UNDEFINED))
    {
        screen::raise_warning("Found a type without a paired value. Falling back to default.",
                              Non_Specified_Value_Warning);

        output = INF_DOUBLE; // default value
    }
    else
    {
        if (typ == TYPE_SCALAR)
        {
            output = std::stod(val);
        }
        else if (typ == TYPE_DATA)
        {
            output = get_data(h, i, val);
        }
        else if (typ == TYPE_PARAM)
        {
            output = get_param(val);
            if (get_infer(val))
            {
                std::ostringstream oss;
                // UPDATE : this error is not clear as exposed here
                oss << "Parameter to be inferred `" << val << "` is chosen for a fixed value move."
                    << "Eventual changes of the parameter will not affect the move.";
                screen::raise_warning(oss.str(), Param_Infer_Fixed_Move_Warning);
            }
        }
        else if (typ == TYPE_COMMAND)
        {
            if (val == "default")
            {
                output = INF_DOUBLE;
            }
            else
            {
                screen::raise_error("Unknown command type `" + val + "`.");
            }
        }
        else
        {
            screen::raise_error("Unknown variable type specifier `" + typ + "`.");
        }
    }

    return output;
}

/* Convert a type-string to a double-pointer */
const double *MCMC::type_to_doubleptr(const std::string val, const std::string typ, int h, int i)
{
    const double *output = nullptr;

    if ((val == UNDEFINED) && (typ == UNDEFINED))
    {
        output = nullptr; // default value
    }
    else if ((val != UNDEFINED) && (typ == UNDEFINED))
    {
        std::ostringstream oss;
        oss << "Impossible to interpret default value of  `" << val << "`. Please provide a type.";
        screen::raise_error(oss.str());
    }
    else if ((val == UNDEFINED) && (typ != UNDEFINED))
    {
        screen::raise_warning("Found a type without a paired value. Falling back to default.",
                              Non_Specified_Value_Warning);

        output = nullptr; // default value
    }
    else
    {
        if (typ == TYPE_DATA)
        {
            output = househ[h].indiv[i].get_data_address(val);
        }
        else if (typ == TYPE_PARAM)
        {
            output = get_param_address(val);
        }
        else
        {
            std::ostringstream oss;
            oss << "Unsopported type specifier `" + typ + "` found while loading `" << val << "`. "
                << "Expected type specifier `" << TYPE_PARAM << "` instead.";
            screen::raise_error(oss.str());
        }
    }

    return output;
}

/* Convert a string with type-dictionary to string-to-double */
sdargs MCMC::typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, int h, int i)
{
    std::vector<std::string> keys = parse_list_args(tkn_key);
    std::vector<std::string> vals = parse_list_args(tkn_val);
    std::vector<std::string> typs = parse_list_args(tkn_typ);

    sdargs output = {};

    if (!((keys.size() == vals.size()) && (typs.size() == vals.size())))
    {
        // inform augmented variable about the error
        output.emplace(IO_MISMATCH, INF_DOUBLE);
    }
    else
    {
        for (size_t x = 0; x < keys.size(); x++)
        {
            output.emplace(keys[x], type_to_double(vals[x], typs[x], h, i));
        }
    }

    return output;
}

/* Convert a string with type-dictionary to string-to-double-pointer */
sdargs_ptr MCMC::typedict_to_sdargs_ext(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, int h, int i)
{
    std::vector<std::string> keys = parse_list_args(tkn_key);
    std::vector<std::string> vals = parse_list_args(tkn_val);
    std::vector<std::string> typs = parse_list_args(tkn_typ);

    sdargs_ptr output = {};

    if (!((keys.size() == vals.size()) && (typs.size() == vals.size())))
    {
        // inform augmented variable about the error
        output.emplace(IO_MISMATCH, nullptr);
    }
    else
    {
        for (size_t x = 0; x < keys.size(); x++)
        {
            output.emplace(keys[x], type_to_doubleptr(vals[x], typs[x], h, i));
        }
    }

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : build augmented  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file */
void MCMC::build_augmented(const std::string pathtofile, bool verbose)
{
    if (pathtofile != UNDEFINED)
    {
        std::ifstream input(pathtofile.c_str());
        std::vector<std::string> required_headers = {"name", "value", "valuetype",
                                                     "model", "args", "argsval", "argstype", "rule", "init_rand",
                                                     "rargs", "rargsval", "rargstype", "rfrom", "rfromtype",
                                                     "min", "mintype", "max", "maxtype"};

        if (input.is_open())
        {
            screen::print_decorator("Augmented", screen::C_DECOR, "import: " + pathtofile, verbose);

            //===== Load column headers =====//
            std::vector<std::string> headers = read_column_headers(input, pathtofile, required_headers);

            //===== Load augmented from the row with options =====//
            load_augm_definition.clear();
            augm_full_set.clear();
            augm_to_be_inferred.clear();
            augm_with_model_set.clear();
            augm_rand_init_set.clear();
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row); this_row.clear())
            {
                load_augm_definition.push_back(this_row);
            }

            assign_augmented_variables();
        }
        else
        {
            screen::raise_error("Impossible to open augmented file `" + pathtofile + "`.");
        }
    }
}

/* Assign augmented variables to all individuals. */
void MCMC::assign_augmented_variables()
{
    std::string augmID;

    // load the definitions of agumetned variablews and define lists of usage
    for (ssargs loader : load_augm_definition)
    {
        augmID = loader.at("name");

        if (augm_full_set.find(augmID) == augm_full_set.end())
        {
            augm_full_set.emplace(augmID);
        }
        else
        {
            std::ostringstream oss;
            oss << "Ambigous definitions for augmented variable `" << augmID << "`."
                << " Correct your augmented variable file.";
            screen::raise_error(oss.str());
        }

        // list of augmented variables to be inferred
        if (loader.at("rule") != UNDEFINED)
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
    __is_data_augm = !augm_to_be_inferred.empty();

    //===== Eventually assign augmented variables to each individual =====//
    if (!load_augm_definition.empty())
    {
        for (int h = 0; h < num_househ; h++)
        {
            for (int i = 0; i < househ[h].house_size; i++)
            {
                for (ssargs loader : load_augm_definition)
                {
                    // load range
                    sdargs range = {};
                    if (type_to_double(loader.at("min"), loader.at("mintype"), h, i) != INF_DOUBLE)
                    {
                        range.emplace("min", type_to_double(loader.at("min"), loader.at("mintype"), h, i));
                    }
                    if (type_to_double(loader.at("max"), loader.at("maxtype"), h, i) != INF_DOUBLE)
                    {
                        range.emplace("max", type_to_double(loader.at("max"), loader.at("maxtype"), h, i));
                    }

                    // assign augmented variable
                    househ[h].indiv[i].assign_augm(
                        loader.at("name"),
                        type_to_double(loader.at("value"), loader.at("valuetype"), h, i),
                        loader.at("model"),
                        typedict_to_sdargs_ext(loader.at("args"), loader.at("argsval"), loader.at("argstype"), h, i),
                        loader.at("rule"),
                        typedict_to_sdargs(loader.at("rargs"), loader.at("rargsval"), loader.at("rargstype"), h, i),
                        type_to_double(loader.at("rfrom"), loader.at("rfromtype"), h, i),
                        range);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : build extra  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Extra File (default unused) */
void MCMC::build_extra(const std::string pathtofile, bool verbose)
{
    UNUSED(pathtofile);
    UNUSED(verbose);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MCMC : BUILD  ///////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::build(const std::string setup_file, const std::string params_file,
                 const std::string data_file, const std::string augm_file,
                 const std::string extra_file, bool verbose)
{
    //===== reset internal status =====//
    __is_infer_complete = false;

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