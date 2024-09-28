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

/* Convert string to boolean variable. */
bool stobool(const string &tkn)
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
        raise_error("Unrecognized boolean value for `" + tkn + "`.");
    }

    return output;
}

/* Clean string from forbidden characters and convert to UNDEFINED. */
string string_cleaner(string tkn)
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

/* Read column headers. */
vector<string> read_column_headers(std::ifstream &input,
                                   const string &pathtofile,
                                   const vector<string> &required_headers)
{
    string line, tkn;
    getline(input, line);
    std::istringstream iss(line);
    vector<string> headers = {};

    while (getline(iss, tkn, C_DELIM))
    {
        headers.push_back(string_cleaner(tkn));
    }

    // check that all required headers are present
    if (!required_headers.empty())
    {
        bool is_missing = false;
        std::ostringstream oss;

        // check the validity of the headers
        for (const string &req : required_headers)
        {
            if (std::find(headers.begin(), headers.end(), req) == headers.end())
            {
                if (!is_missing)
                {
                    // point to the faulty file
                    oss << "Missing column(s): `" + req + "`";
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
            raise_error(oss.str() + ". Check file " + pathtofile);
        }
    }

    return headers;
}

/* Load files from the row with headers */
bool read_row_by_headers(std::ifstream &input, const vector<string> &headers, ssargs &this_row)
{
    string line, tkn;

    if (getline(input, line))
    {
        std::istringstream iss(line);
        size_t num_cols = 0;

        while (getline(iss, tkn, C_DELIM))
        {
            this_row.emplace(headers[num_cols++], string_cleaner(tkn));
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

/* convert the output of a row to a string-to-double dictionary */
dictargs row_to_dictargs(const ssargs &this_row)
{
    dictargs row_with_double = {};
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

/* Parse the list from the string. */
vector<string> parse_list_args(const string &tkn, const char delim)
{
    std::istringstream args_ss(tkn);
    vector<string> output = {};
    string tmp;

    // rea the entries of the dictionary
    while (getline(args_ss, tmp, delim))
    {
        output.push_back(tmp);
    }

    return output;
}

/* Combine after control two lists into a dictionary */
dictargs merge_list_args(const string &tkn_keys, const string &tkn_vals)
{
    vector<string> keys = parse_list_args(tkn_keys);
    vector<string> vals = parse_list_args(tkn_vals);
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
