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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

bool stobool(const string &tkn);
string string_cleaner(string tkn);

#ifdef USE_RCPP
void rlist_to_map(Rcpp::List list, std::map<string, string> &output);
void rdataframe_to_map(Rcpp::DataFrame df, std::map<string, vector<string>> &output);
#endif

vector<string> read_column_headers(std::ifstream &input, const vector<string> &required_headers);
bool read_row_by_headers(std::ifstream &input, const vector<string> &headers, std::map<string, string> &this_row);
void read_csv_list(const string &pathtofile, std::map<string, string> &csv_parsed_list);
void read_csv_table(const string &pathtofile,
                    std::map<string, vector<string>> &csv_parsed_table,
                    const vector<string> &required_headers = {});

vector<string> parse_list_args(const string &tkn, const char delim = basics::C_LIST_DELIM);
dictargs merge_list_args(const string &tkn_keys, const string &tkn_vals);

#endif // READ__H