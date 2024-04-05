// build_aux_read.cpp //

#include "mcmc.h"

// declaration for type in data file
std::string TYPE_PARAM = "P";
std::string TYPE_DATA = "D";
std::string TYPE_SCALAR = "S";

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
        raise_error_message(err_text);
        std::exit(EXIT_FAILURE);
    }
}

/* Parse the dictionary from the string */
ssargs MCMC::parse_dict_string(std::string token, char dict_delim)
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
        raise_error_message(err_text);
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
                std::string err_text = "Invalid entry of dictionary: " + token + ".";
                raise_error_message(err_text);
                std::exit(EXIT_FAILURE);
            }
        }
    }

    return output;
}

/* Parse the scalar arguments dictionary from the std::string */
sdargs MCMC::parse_dict_args(std::string token, char dict_delim)
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
std::vector<std::string> MCMC::read_column_headers(std::ifstream &input, char delim)
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
bool MCMC::read_row_by_headers(std::ifstream &input, std::vector<std::string> headers, ssargs &this_row, char delim)
{
    std::string line, token;
    if (getline(input, line))
    {
        std::istringstream iss(line);
        int num_cols = 0;

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
double MCMC::type_to_double(std::string token, int i, char option_delim)
{
    if (token == UNDEFINED)
    {
        // UPDATE : when does it make sense to have token == UNDEFINED ?
        return INFINITY;
    }
    else
    {
        int pos = token.find(option_delim);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_SCALAR)
        {
            return std::stod(v);
        }
        else if (k == TYPE_DATA)
        {
            return indiv[i].get_data(v);
        }
        else if (k == TYPE_PARAM)
        {
            return get_value(v);
        }
        else
        {
            std::string err_text = "Unknown variable type specifier `" + k + "`.";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
    }
}

/* Convert a type-string to a double-pointer */
const double *MCMC::type_to_doubleptr(std::string token, char option_delim)
{
    if (token == UNDEFINED)
    {
        return nullptr;
    }
    else
    {
        int pos = token.find(option_delim);
        std::string k = token.substr(0, pos);
        std::string v = token.substr(pos + 1);
        if (k == TYPE_PARAM)
        {
            return get_address(v);
        }
        else
        {
            std::string err_text = "Expected variable type `" + TYPE_PARAM + "` while loading `" + token + "`. Specifier `" + k + "` not supported here.";
            raise_error_message(err_text);
            std::exit(EXIT_FAILURE);
        }
    }
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