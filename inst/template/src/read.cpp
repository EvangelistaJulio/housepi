/*
 * read.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "read.h"

using namespace basics;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Convert std::string to boolean variable. */
bool stobool(const std::string &tkn)
{
    bool output = false;

    if (tkn == "TRUE" || tkn == "true" || tkn == "True" || tkn == "1")
    {
        output = true;
    }
    else if (tkn == "FALSE" || tkn == "false" || tkn == "False" || tkn == "0")
    {
        output = false;
    }
    else
    {
        raise_error("Unrecognized boolean value for '" + tkn + "'.");
    }

    return output;
}

/* Clean std::string from forbidden characters and convert to UNDEFINED. */
std::string strcleaner(std::string tkn)
{
    // forbidden characters
    for (char r : {'"', '\r', '\n'})
    {
        tkn.erase(remove(tkn.begin(), tkn.end(), r), tkn.end());
    }

    // check for undefined
    if (is_na(tkn))
    {
        tkn = UNDEFINED;
    }

    return tkn;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef USE_RCPP
void rlist_to_map(Rcpp::List list, std::map<std::string, std::string> &output)
{
    using namespace Rcpp;

    CharacterVector keys = list.names();

    output.clear();
    for (int i = 0; i < list.size(); ++i)
    {
        RObject element = list[i];
        output[as<std::string>(keys[i])] = as<std::string>(wrap(element));
    }
}

void rdataframe_to_map(Rcpp::DataFrame df, std::map<std::string, std::vector<std::string>> &output)
{
    using namespace Rcpp;

    CharacterVector keys = df.names();

    output.clear();
    for (int i = 0; i < df.size(); ++i)
    {
        CharacterVector col_str = as<CharacterVector>(df[i]);
        std::vector<std::string> col_data(col_str.begin(), col_str.end());
        output[as<std::string>(keys[i])] = std::move(col_data);
    }
}
#endif

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

std::vector<std::string> read_column_headers(std::ifstream &input, const std::vector<std::string> &required_headers)
{
    // clean headers
    std::vector<std::string> headers = {};

    // read first line
    std::string line, tkn;
    getline(input, line);
    std::istringstream iss(line);

    // get headers from first line
    while (getline(iss, tkn, C_DELIM))
    {
        headers.push_back(strcleaner(tkn));
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
                    oss << "Missing column(s): '" + req + "'";
                    is_missing = true;
                }
                else
                {
                    oss << ", '" << req << "'";
                }
            }
        }

        // raise the error if something was missing
        if (is_missing)
        {
            raise_error(oss.str());
        }
    }

    return headers;
}

/* Load files from the row with headers */
bool read_row_by_headers(std::ifstream &input, const std::vector<std::string> &headers, std::map<std::string, std::string> &this_row)
{
    std::string line, tkn;

    if (getline(input, line))
    {
        std::istringstream iss(line);
        size_t num_cols = 0;

        while (getline(iss, tkn, C_DELIM))
        {
            this_row.emplace(headers[num_cols++], strcleaner(tkn));
        }

        // last column is empty
        for (size_t i = num_cols; i < headers.size(); i++)
        {
            this_row.emplace(headers[i], UNDEFINED);
        }

        return true;
    }
    else
    {
        return false;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Parse the list from the std::string. */
std::vector<std::string> parse_list_args(const std::string &tkn, const char delim)
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
dictargs merge_list_args(const std::string &tkn_keys, const std::string &tkn_vals)
{
    std::vector<std::string> keys = parse_list_args(tkn_keys);
    std::vector<std::string> vals = parse_list_args(tkn_vals);
    dictargs output = {};

    if (keys.size() != vals.size())
    {
        // inform parameter about the error
        output.emplace(IO_MISMATCH, INF_DOUBLE);
    }
    else
    {
        for (size_t i = 0; i < keys.size(); ++i)
        {
            output.emplace(keys[i], std::stod(vals[i]));
        }
    }
    return output;
}

/* Read csv list store in pathtofile and store it parsed in csv_parsed_list. */
void read_csv_list(const std::string &pathtofile, std::map<std::string, std::string> &csv_parsed_list)
{
    if (!is_na(pathtofile))
    {
        std::ifstream input(pathtofile.c_str());

        if (input.is_open())
        {
            csv_parsed_list.clear();

            // Create a std::stringstream of the current line
            for (std::string line; getline(input, line);)
            {
                std::string name, choice;
                std::istringstream iss(line);

                //===== NAME =====//
                if (!getline(iss, name, C_DELIM))
                {
                    raise_error("Unable to extract list name. Check file line: " + line);
                }

                //===== CHOICE =====//
                if (!getline(iss, choice))
                {
                    raise_error("Unable to extract list value. Check file line: " + line);
                }

                // add the new setup
                csv_parsed_list.emplace(strcleaner(name), strcleaner(choice));
            }
        }
        else
        {
            raise_error("Impossible to open file '" + pathtofile + "'.");
        }
    }
}

/* Read csv table store in pathtofile and store it parsed in csv_parsed_table. */
void read_csv_table(const std::string &pathtofile,
                    std::map<std::string, std::vector<std::string>> &csv_parsed_table,
                    const std::vector<std::string> &required_headers)
{
    if (!is_na(pathtofile))
    {
        std::ifstream input(pathtofile.c_str());

        if (input.is_open())
        {
            std::vector<std::string> headers = read_column_headers(input, required_headers);
            std::map<std::string, std::string> this_row = {};

            csv_parsed_table.clear(); // FIXME: is this correct?
            while (read_row_by_headers(input, headers, this_row))
            {
                for (const auto &x : this_row)
                {
                    csv_parsed_table[x.first].push_back(x.second);
                }
                this_row.clear();
            }
        }
        else
        {
            raise_error("Impossible to open file '" + pathtofile + "'.");
        }
    }
}
