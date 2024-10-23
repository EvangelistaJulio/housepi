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

bool stobool(const std::string &tkn);
std::string strcleaner(std::string tkn);

#ifdef USE_RCPP
void rlist_to_map(Rcpp::List list, std::map<std::string, std::string> &output);
void rdataframe_to_map(Rcpp::DataFrame df, std::map<std::string, std::vector<std::string>> &output);
#endif

std::vector<std::string> read_column_headers(std::ifstream &input, const std::vector<std::string> &required_headers);
bool read_row_by_headers(std::ifstream &input, const std::vector<std::string> &headers, std::map<std::string, std::string> &this_row);
void read_csv_list(const std::string &pathtofile, std::map<std::string, std::string> &csv_parsed_list);
void read_csv_table(const std::string &pathtofile,
                    std::map<std::string, std::vector<std::string>> &csv_parsed_table,
                    const std::vector<std::string> &required_headers = {});

std::vector<std::string> parse_list_args(const std::string &tkn, const char delim = tools::C_LIST_DELIM);
dictargs merge_list_args(const std::string &tkn_keys, const std::string &tkn_vals);

#endif // READ__H