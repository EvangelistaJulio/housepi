// build-read.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "mcmc.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

bool Param_Infer_Fixed_Move_Warning = false;

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  MCMC : methods (auxiliary read functions)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// declaration for type in data file
const std::string TYPE_PARAM = "P";   // Parameter type declaration
const std::string TYPE_DATA = "D";    // Data type declaration
const std::string TYPE_SCALAR = "S";  // Scalar type declaration
const std::string TYPE_COMMAND = "C"; // Command type declaration

/* Cleaner std::string from forbidden characters. */
std::string MCMC::string_cleaner(std::string token)
{
    for (char r : {'\r', '\n'})
    {
        token.erase(remove(token.begin(), token.end(), r), token.end());
    }
    return token;
}

/* Convert string to boolean variable. */
bool MCMC::stobool(std::string token)
{
    bool output = false;

    if (token == "TRUE" || token == "true" || token == "True" || token == "1")
    {
        output = true;
    }
    else if (token != "FALSE" || token != "false" || token != "False" || token != "0")
    {
        raise_error_message("Unrecognized boolean value for: " + token + ".");
    }

    return output;
}

/* Parse the dictionary from the string. */
ssargs MCMC::parse_dict_string(std::string token)
{
    std::istringstream args_ss(token);
    std::vector<std::string> args_stock;
    ssargs output = {};
    std::string tmp;

    // rea the entryes of the dictionary
    while (getline(args_ss, tmp, C_DICT_DELIM))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output.clear(); // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
    {
        raise_error_message("Incomplete dictionary: " + token + ".");
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
                raise_error_message("Invalid entry of dictionary: " + token + ".");
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

    while (getline(iss, token, C_DELIM))
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

        while (getline(iss, token, C_DELIM))
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

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : methods (type-defined variables)  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a type-string to a double */
double MCMC::type_to_double(std::string token, int h, int i)
{
    double output = INF_DOUBLE; // default value

    // Note : when does it make sense to have token == UNDEFINED ?
    if (token != UNDEFINED)
    {
        int pos = token.find(C_OPTION_DELIM);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);

        if (k == TYPE_SCALAR)
        {
            output = std::stod(v);
        }
        else if (k == TYPE_DATA)
        {
            househ[h].indiv[i].assert_data(v);
            output = househ[h].indiv[i].get_data(v);
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
            output = get_value(v);
        }
        else if (k == TYPE_COMMAND)
        {
            if (v != "default")
            {
                raise_error_message("Unknown command `" + v + "`.");
            }
        }
        else
        {
            raise_error_message("Unknown variable type specifier `" + k + "`.");
        }
    }

    return output;
}

/* Convert a type-string to a double-pointer */
const double *MCMC::type_to_doubleptr(std::string token)
{
    const double * output = nullptr;
    
    if (token != UNDEFINED)
    {
        int pos = token.find(C_OPTION_DELIM);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_PARAM)
        {
            assert_param(v);
            output = get_address(v);
        }
        else
        {
            std::ostringstream oss;
            oss << "Unsopported type specifier `" + k + "` found while loading `" << token << "`. "
                << "Expected type specifier `" << TYPE_PARAM << "` instead.";
            raise_error_message(oss);
        }
    }

    return output;
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

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//