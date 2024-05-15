// build-read.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

// declaration for type in data file
std::string TYPE_PARAM = "P";   // Parameter type declaration
std::string TYPE_DATA = "D";    // Data type declaration
std::string TYPE_SCALAR = "S";  // Scalar type declaration
std::string TYPE_COMMAND = "C"; // Command type declaration

bool Param_Infer_Fixed_Move_Warning = false;

///////////////////////////////////////////////////////////////////////
// ///////////////////// Auxiliary Functions /////////////////////// //
///////////////////////////////////////////////////////////////////////

/* Cleaner std::string from forbidden characters*/
std::string MCMC::string_cleaner(std::string token)
{
    for (char r : {'\r', '\n'})
    {
        token.erase(remove(token.begin(), token.end(), r), token.end());
    }
    return token;
}

/* Convert string to boolean variable */
bool MCMC::stobool(std::string token)
{
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
        display_error_message("Unrecognized boolean value for: " + token + ".");
        std::exit(EXIT_FAILURE);
    }
}

/* Parse the dictionary from the string */
ssargs MCMC::parse_dict_string(std::string token)
{
    std::istringstream args_ss(token);
    std::vector<std::string> args_stock;
    ssargs output = {};
    std::string tmp;

    // rea the entryes of the dictionary
    while (getline(args_ss, tmp, DICT_DELIM))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output.clear(); // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
    {
        display_error_message("Incomplete dictionary: " + token + ".");
        std::exit(EXIT_FAILURE);
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
                display_error_message("Invalid entry of dictionary: " + token + ".");
                std::exit(EXIT_FAILURE);
            }
        }
    }

    return output;
}

/* Parse the scalar arguments dictionary from the std::string */
sdargs MCMC::parse_dict_args(std::string token)
{
    ssargs input = parse_dict_string(token);
    sdargs output = {};
    for (auto const &x : input)
    {
        output.emplace(x.first, std::stod(x.second));
    }
    return output;
}

/* Read column headers */
std::vector<std::string> MCMC::read_column_headers(std::ifstream &input)
{
    std::string line, token;
    getline(input, line);
    std::istringstream iss(line);
    std::vector<std::string> headers = {};

    while (getline(iss, token, DELIM))
    {
        headers.push_back(string_cleaner(token));
    }

    return headers;
}

/* Load files from the row with headers */
bool MCMC::read_row_by_headers(std::ifstream &input, std::vector<std::string> headers, ssargs &this_row)
{
    std::string line, token;
    if (getline(input, line))
    {
        std::istringstream iss(line);
        int num_cols = 0;

        while (getline(iss, token, DELIM))
        {
            this_row.emplace(headers[num_cols], string_cleaner(token));
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

/* Convert a type-string to a double */
double MCMC::type_to_double(std::string token, int h, int i)
{
    if (token == UNDEFINED)
    {
        // UPDATE : when does it make sense to have token == UNDEFINED ?
        return INF_DOUBLE;
    }
    else
    {
        int pos = token.find(OPTION_DELIM);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_SCALAR)
        {
            return std::stod(v);
        }
        else if (k == TYPE_DATA)
        {
            househ[h].indiv[i].assert_data(v);
            return househ[h].indiv[i].get_data(v);
        }
        else if (k == TYPE_PARAM)
        {
            assert_param(v);
            if (get_infer(v))
            {
                std::ostringstream oss;
                oss << "The parameter `" << v << "` is chosen for a fixed value move."
                    << "Eventual changes of the parameter will not affect the move.";
                raise_warning_message(oss, &Param_Infer_Fixed_Move_Warning);
            }
            return get_value(v);
        }
        else if (k == TYPE_COMMAND)
        {
            if (v == "default")
            {
                return INF_DOUBLE;
            }
            else
            {
                display_error_message("Unknown command `" + v + "`.");
                std::exit(EXIT_FAILURE);
            }
        }
        else
        {
            display_error_message("Unknown variable type specifier `" + k + "`.");
            std::exit(EXIT_FAILURE);
        }
    }
}

/* Convert a type-string to a double-pointer */
const double *MCMC::type_to_doubleptr(std::string token)
{
    if (token == UNDEFINED)
    {
        return nullptr;
    }
    else
    {
        int pos = token.find(OPTION_DELIM);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_PARAM)
        {
            assert_param(v);
            return get_address(v);
        }
        else
        {
            std::ostringstream oss;
            oss << "Unsopported type specifier `" + k + "` found while loading `" << token << "`. "
                << "Expected type specifier `" << TYPE_PARAM << "` instead.";
            display_error_message(oss);
            std::exit(EXIT_FAILURE);
        }
    }
}

/* Convert a string with type-dictionary to string-to-double */
sdargs MCMC::typedict_to_sdargs(std::string token, int h, int i)
{
    sdargs output = {};
    if (token != UNDEFINED)
    {
        ssargs input = parse_dict_string(token);
        for (auto const &x : input)
        {
            output.emplace(x.first, type_to_double(x.second, h, i));
        }
    }
    return output;
}

/* Convert a string with type-dictionary to string-to-double-pointer */
sdargs_ext MCMC::typedict_to_sdargs_ext(std::string token)
{
    sdargs_ext output = {};
    if (token != UNDEFINED)
    {
        ssargs input = parse_dict_string(token);
        for (auto const &x : input)
        {
            output.emplace(x.first, type_to_doubleptr(x.second));
        }
    }
    return output;
}