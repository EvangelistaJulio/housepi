/*
 * read.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef READ__H
#define READ__H

#include <map>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "utils.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::map<string, string> ssargs;

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

bool stobool(const string &tkn);
string string_cleaner(string tkn);

vector<string> read_column_headers(std::ifstream &input, const string &pathtofile, const vector<string> &required_headers = {});
bool read_row_by_headers(std::ifstream &input, const vector<string> &headers, ssargs &this_row);
dictargs row_to_dictargs(const ssargs &this_row);

vector<string> parse_list_args(const string &tkn, const char delim = basics::C_LIST_DELIM);
dictargs merge_list_args(const string &tkn_keys, const string &tkn_vals);

#endif // READ__H