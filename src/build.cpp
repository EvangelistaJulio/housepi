// build.cpp //

#include "mcmc.h"

// declaration for type in data file
std::string TYPE_PARAM = "P";
std::string TYPE_DATA = "D";
std::string TYPE_SCALAR = "S";

///////////////////////////////////////////////////////////////////////
// ///////////////////// Auxiliary Functions /////////////////////// //
///////////////////////////////////////////////////////////////////////

/* Cleaner std::string from forbidden characters*/
std::string string_cleaner(std::string token)
{
    for (char r : {'\r', '\n'})
    {
        token.erase(remove(token.begin(), token.end(), r), token.end());
    }
    return token;
}

/* Convert string to boolean variable */
bool stobool(std::string token)
{
    token = string_cleaner(token);

    if (token == "TRUE" || token == "true" || token == "True" || token == "1")
    {
        return true;
    }
    else if (token == "FALSE" || token == "false" || token == "False" || token == "0")
    {
        return false;
    }
    else
    {
        std::string err_text = "Unrecognized boolean value for: " + token + ".";
        exit_with_failure_error(err_text);
        return false;
    }
}

/* Parse the dictionary from the string */
ssargs parse_dict_string(std::string token, char dict_delim)
{
    std::istringstream args_ss(token);
    std::vector<std::string> args_stock;
    ssargs output = {};
    std::string tmp;

    // rea the entryes of the dictionary
    while (getline(args_ss, tmp, dict_delim))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output.clear(); // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
    {
        std::string err_text = "Incomplete dictionary: " + token + ".";
        exit_with_failure_error(err_text);
    }
    else
    {
        for (std::size_t i = 0; i < args_stock.size(); i += 2)
        {
            try
            {
                output.emplace(args_stock[i], args_stock[i + 1]);
            }
            catch (const std::exception &)
            {
                std::string err_text = "Invalid entry of dictionary: " + token + ".";
                exit_with_failure_error(err_text);
            }
        }
    }

    return output;
}

/* Parse the scalar arguments dictionary from the std::string */
sdargs parse_dict_args(std::string token, char dict_delim)
{
    ssargs input = parse_dict_string(token, dict_delim);
    sdargs output = {};
    for (auto const &x : input)
    {
        output.emplace(x.first, std::stod(x.second));
    }
    return output;
}

/* Read column headers */
std::vector<std::string> read_column_headers(std::ifstream &input, char delim = DELIM)
{
    std::string line, token;
    getline(input, line);
    std::istringstream iss(line);
    std::vector<std::string> headers = {};

    while (getline(iss, token, delim))
    {
        token = string_cleaner(token);
        headers.push_back(token);
    }

    return headers;
}

/* Load parameters from the row with options */
bool read_row_by_headers(std::ifstream &input, std::vector<std::string> headers, ssargs &this_row, char delim = DELIM)
{
    std::string line, token;
    if (getline(input, line))
    {
        std::istringstream iss(line);
        std::size_t num_cols = 0;

        while (getline(iss, token, delim))
        {
            this_row.emplace(headers[num_cols], token);
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
sdargs row_to_sdargs(ssargs &this_row)
{
    sdargs row_with_double = {};
    for (auto const &x : this_row)
    {
        row_with_double.emplace(x.first, stod(x.second));
    }
    return row_with_double;
}

/* Convert a type-string to a double */
double MCMC::type_to_double(std::string token, int i, char option_delim)
{
    double output = INFINITY;
    if (token != UNDEFINED)
    {
        std::size_t pos = token.find(option_delim);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_SCALAR)
        {
            output = std::stod(v);
        }
        else if (k == TYPE_DATA)
        {
            output = indiv[i].get_data(v);
        }
        else if (k == TYPE_PARAM)
        {
            output = get_value(v);
        }
        else
        {
            std::string err_text = "Unknown variable type specifier `" + k + "`.";
            exit_with_failure_error(err_text);
        }
    }
    return output;
}

/* Convert a type-string to a double-pointer */
const double *MCMC::type_to_doubleptr(std::string token, char option_delim)
{
    const double *output = nullptr;
    if (token != UNDEFINED)
    {
        std::size_t pos = token.find(option_delim);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_PARAM)
        {
            output = get_address(v);
        }
        else
        {
            std::string err_text = "Expected variable type `" + TYPE_PARAM + "` while loading `" + token + "`. Specifier `" + k + "` not supported here.";
            exit_with_failure_error(err_text);
        }
    }
    return output;
}

/* Convert a string with type-dictionary to string-to-double */
sdargs MCMC::typedict_to_sdargs(std::string token, int, char dict_delim, char option_delim)
{
    ssargs input = parse_dict_string(token, dict_delim);
    sdargs output = {};
    for (auto const &x : input)
    {
        output.emplace(x.first, type_to_double(x.second, 0, option_delim));
    }
    return output;
}

/* Convert a string with type-dictionary to string-to-double-pointer */
sdargs_ext MCMC::typedict_to_sdargs_ext(std::string token, char dict_delim, char option_delim)
{
    ssargs input = parse_dict_string(token, dict_delim);
    sdargs_ext output = {};
    for (auto const &x : input)
    {
        output.emplace(x.first, type_to_doubleptr(x.second, option_delim));
    }
    return output;
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////// Load Setup File /////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Setup file */
Setup_MCMC_Struct load_setup(std::string data_pathtofile)
{
    char delim = ',';
    Setup_MCMC_Struct loader;
    std::ifstream input(data_pathtofile.c_str());
    ssargs safe_load = {};

    if (input.is_open())
    {
        if (Verbose)
        {
            std::cout << "\n================SETUP=================" << std::endl;
            std::cout << "Setup file: " << data_pathtofile << "\n"
                      << std::endl;
        }

        // Create a std::stringstream of the current line
        for (std::string line; getline(input, line);)
        {
            std::string name, choice;
            std::istringstream iss(line);

            //=====NAME=====//
            if (!getline(iss, name, delim))
            {
                std::string err_text = "Unable to extract setup name.\nCheck file line: " + line;
                exit_with_failure_error(err_text);
            }

            //=====CHOICE=====//
            if (!getline(iss, choice))
            {
                std::string err_text = "Unable to extract setup choice.\nCheck file line: " + line;
                exit_with_failure_error(err_text);
            }

            // add the new setup
            safe_load.emplace(string_cleaner(name), string_cleaner(choice));
        }

        std::string seed_str = safe_load.at("seed");
        loader.seed = (seed_str == UNDEFINED) ? time(NULL) : stoi(seed_str);
        loader.num_iter = stoi(safe_load.at("num_iter"));
        loader.num_adapt_iter = stoi(safe_load.at("num_adapt_iter"));
        loader.thinning = stoi(safe_load.at("thinning"));
        loader.delta_t = stod(safe_load.at("delta_t"));

        //=====Display=====//
        if (Verbose)
        {
            std::cout << "seed : " << loader.seed << std::endl;
            std::cout << "num_iter : " << loader.num_iter << std::endl;
            std::cout << "num_adapt_iter : " << loader.num_adapt_iter << std::endl;
            std::cout << "thinning : " << loader.thinning << std::endl;
            std::cout << "delta_t : " << loader.delta_t << std::endl;
        }
    }
    else
    {
        std::string err_text = "Impossible to read file: " + data_pathtofile + ".";
        exit_with_failure_error(err_text);
    }

    return loader;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// MCMC : load parameters ///////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Parameters file */
void MCMC::build_parameters(std::string data_pathtofile)
{
    std::ifstream input(data_pathtofile.c_str());

    if (input.is_open())
    {
        if (Verbose)
        {
            std::cout << "\n==============PARAMETERS===============" << std::endl;
            std::cout << "Parameters file: " << data_pathtofile << "\n"
                      << std::endl;
        }

        //=====Load column headers=====//
        std::vector<std::string> headers = read_column_headers(input, DELIM);
        // check the validity of the headers
        for (std::string req : param_file_required_headers)
        {
            if (std::find(headers.begin(), headers.end(), req) == headers.end())
            {
                std::string err_text = "Missing column `" + req + "`. Please check file : " + data_pathtofile + ".";
                exit_with_failure_error(err_text);
            }
        }

        //=====Load parameters from the row with options=====//
        for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row, DELIM); this_row.clear())
        {
            // add entry to shared parameter dictionary
            if (stobool(this_row.at("infer")))
            {
                add_param(
                    this_row["name"],
                    std::stod(this_row["init_value"]),
                    this_row["prior"],
                    parse_dict_args(this_row["prior_args"], DICT_DELIM),
                    this_row["strategy"],
                    parse_dict_args(this_row["strategy_args"], DICT_DELIM),
                    parse_dict_args(this_row["range"], DICT_DELIM));
            }
            else
            {
                add_param(this_row.at("name"), stod(this_row.at("init_value")));
            }
        }

        //=====Initialize the uniform parameter picker=====//
        param_picker = std::uniform_int_distribution<std::size_t>(0, num_param_to_infer - 1);

        //=====display a summary of the used parameters=====//
        if (Verbose)
        {
            display_parameters();
        }
    }
    else
    {
        std::string err_text = "Impossible to read file:" + data_pathtofile + ".";
        exit_with_failure_error(err_text);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////  MCMC : load augmented //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Model Augmented Variable file */
void MCMC::build_augmented(std::string data_pathtofile)
{
    std::ifstream input(data_pathtofile.c_str());
    std::vector<std::string> required_headers = {
        "name", "init_value", "model", "model_args", "strategy", "strategy_args", "move_from", "range"};

    if (data_pathtofile != UNDEFINED)
    {
        if (input.is_open())
        {
            if (Verbose)
            {
                std::cout << "\n==============AUGMENTED===============" << std::endl;
                std::cout << "Augmented Variable file: " << data_pathtofile << "\n"
                          << std::endl;
            }

            //=====Load column headers=====//
            std::vector<std::string> headers = read_column_headers(input, DELIM);
            // check the validity of the headers
            for (std::string req : required_headers)
            {
                if (std::find(headers.begin(), headers.end(), req) == headers.end())
                {
                    std::string err_text = "Missing column `" + req + "`. Please check file : " + data_pathtofile + ".";
                    exit_with_failure_error(err_text);
                }
            }

            //=====Load augmented from the row with options=====//
            for (ssargs this_row = {}; read_row_by_headers(input, headers, this_row, DELIM); this_row.clear())
            {
                load_augm_definition.push_back(this_row);
            }
        }
        else
        {
            std::string err_text = "Impossible to read file:" + data_pathtofile + ".";
            exit_with_failure_error(err_text);
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
    for (std::size_t i = 0; i < num_indiv; i++)
    {
        int status = indiv[i].get_status();

        for (ssargs loader : load_augm_definition)
        {
            augmID = loader.at("name");
            Augmented augm_var(
                type_to_double(loader.at("init_value"), i),
                loader.at("model_name"),
                typedict_to_sdargs_ext(loader.at("model_args")),
                loader.at("strategy"),
                typedict_to_sdargs(loader.at("strategy_args"), i),
                type_to_double(loader.at("move_from"), i),
                typedict_to_sdargs(loader.at("range"), i));

            // if the individual is not healthy, randomly initializes the augmented variable
            (status == HEALTHY) ? indiv[i].add_augm(augmID, augm_var) : indiv[i].add_augm(augmID, augm_var, rng);
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
// ////////////////////////////////// MCMC : build_data //////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Extra File (default unused) */
void MCMC::build_extra(std::string data_pathtofile)
{
    UNUSED(data_pathtofile);
}

/* Read data file and pre-process the data*/
void MCMC::build_data(std::string data_pathtofile)
{
    std::ifstream input(data_pathtofile.c_str());
    std::vector<std::string> required_headers = {"index", "status", "HH_idx"};

    if (input.is_open())
    {
        if (Verbose)
        {
            std::cout << "\n================DATA=================" << std::endl;
            std::cout << "Data file: " << data_pathtofile << std::endl;
        }

        //=====Load column headers=====//
        data_headers = read_column_headers(input, DELIM);
        // check the validity of the headers
        for (std::string req : required_headers)
        {
            if (std::find(data_headers.begin(), data_headers.end(), req) == data_headers.end())
            {
                std::string err_text = "Missing column `" + req + "`. Please check file : " + data_pathtofile + ".";
                exit_with_failure_error(err_text);
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
            if (this_indiv.get_status() == SYMPTOMATIC)
            {
                num_cases++;
            }
            else if (this_indiv.get_status() != HEALTHY)
            {
                num_unknown_status++;
            }
        }

        //=====assign augmented variables to all individuals=====//
        assign_augmented_variables();

        //=====assign structure of the community=====//
        assign_community_structure("HH_idx"); // FIXME : see above

        //=====initialize contact patterns=====//
        // FIXME : this should be optional
        // FIXME : missing report
        contact_patterns = new const double *[num_indiv * num_indiv];
        std::fill_n(contact_patterns, num_indiv * num_indiv, nullptr);

        for (std::size_t i = 0; i < num_indiv; i++)
        {
            //-----contact patters-----//
            // FIXME : there should be a control also here
            for (std::size_t k = 0; k < i; k++)
            {
                assign_rate_per_contact(i, k);
                assign_rate_per_contact(k, i);
            }
        }

        //=====clean up list of headers for output files=====//
        data_headers.erase(find(data_headers.begin(), data_headers.end(), "index"));
        data_headers.erase(find(data_headers.begin(), data_headers.end(), "status"));

        for (std::string augmID : augm_set)
        {
            auto it = find(data_headers.begin(), data_headers.end(), augmID);
            if (it != data_headers.end())
            {
                data_headers.erase(it);
            }
        }

        //=====Display=====//
        if (Verbose)
        {
            std::cout << "\n============MODEL SUMMARY=============" << std::endl;
            std::cout << "N. of individuals : " << num_indiv << std::endl;
            std::cout << "N. of households : " << num_househ << std::endl;
            std::cout << "N. of cases : " << num_cases << std::endl;
            if (num_unknown_status > 0)
            {
                std::cout << "N. of missing status : " << num_unknown_status << std::endl;
            }

            for (const auto &which : model_report)
            {
                std::cout << which.first << " : " << std::endl;
                for (auto const &x : which.second)
                {
                    std::cout << "    * " << x.first << " : " << x.second << " (" << (100.0 * x.second / num_indiv) << " %)" << std::endl;
                }
            }
        }
        model_report.clear();
    }
    else
    {
        std::string err_text = "Impossible to open file: " + data_pathtofile + ".";
        exit_with_failure_error(err_text);
    }
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////  MCMC : print parameters //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void MCMC::init_param_file(std::string output_filename)
{
    param_ofstream.open(output_filename.c_str());

    if (param_ofstream.is_open())
    {
        if (Verbose)
        {
            std::cout << "Inferred parameter file: " << output_filename << std::endl;
        }

        param_ofstream << "iter" << DELIM << "loglike";
        for (std::string augmID : augm_to_infer_set)
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
        std::string err_text = "Unable to open file `" + output_filename + "`.";
        exit_with_failure_error(err_text);
    }
}

/* Print the model parameters at the given iteration time */
void MCMC::print_parameters(unsigned step)
{
    if (param_ofstream.is_open())
    {
        param_ofstream << step << DELIM << get_global_loglike();
        for (std::string augmID : augm_to_infer_set)
        {
            param_ofstream << DELIM << 100.0 * ((double)m_num_accept_data[augmID] / m_num_propos_data[augmID]);
        }
        for (std::string parID : param_to_be_inferred)
        {
            param_ofstream << DELIM << get_value(parID) << DELIM << 100.0 * get_accept_rate(parID);
        }
        param_ofstream << std::endl;
    }
    else
    {
        std::string err_text = "Output file closed unexpectedly.";
        exit_with_failure_error(err_text);
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
        augm_ofstream << "index" << DELIM << "status";
        for (std::string dataID : data_headers)
        {
            augm_ofstream << DELIM << dataID;
        }
        for (std::string augmID : augm_set)
        {
            augm_ofstream << DELIM << augmID;
        }
        augm_ofstream << std::endl;

        for (std::size_t i = 0; i < num_indiv; i++)
        {
            augm_ofstream << indiv[i].get_index() << DELIM << indiv[i].get_status();
            for (std::string dataID : data_headers)
            {
                augm_ofstream << DELIM << indiv[i].get_data(dataID);
            }
            for (std::string augmID : augm_set)
            {
                augm_ofstream << DELIM << indiv[i].get_augm(augmID);
            }
            augm_ofstream << std::endl;
        }

        augm_ofstream.close();
    }
    else
    {
        std::string err_text = "Unable to open file `" + output_filename + "`.";
        exit_with_failure_error(err_text);
    }
}