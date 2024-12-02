/*
 * build.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__MCMC__BUILD__H
#define HOUSEPI__MCMC__BUILD__H

#include <string>
#include <map>
#include <vector>
#include <memory>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Forward declaration of core implementation. */
class Monte_Carlo_Core;

/* Class to operate on data and augmented variables. */
class Monte_Carlo_Build
{
public:
    //===== Household Data =====//
    size_t n_individuals = 0;
    size_t n_households = 0;
    std::vector<Household> house_vector = {};

    //===== Data Reports =====//
    std::vector<std::string> data_headers = {};
    ulong _n_contacts_ = 0;
    std::map<std::string, size_t> _report_data_;
    std::map<std::string, std::map<std::string, size_t>> _report_group_classes_ = {};
    std::map<std::string, std::string> _report_pathtofile_;

    //===== Augmented Variables =====//            
    std::unordered_set<std::string> augm_list_full = {};
    std::unordered_set<std::string> augm_list_infer = {};
    std::unordered_set<std::string> augm_list_model = {};
    std::map<int, std::unordered_set<std::string>> augm_status_map_infer = {};
    std::map<int, std::unordered_set<std::string>> augm_status_map_model = {};
    bool is_augm_infer = false;
    bool is_augm_model = false;

    //===== CONSTRUCTOR =====//
    Monte_Carlo_Build();
    ~Monte_Carlo_Build();

    //===== IMPORT\EXPORT =====//

#ifdef _HOUSEPI_RCPP
	//===== Rcpp IMPORTS =====//
	void import_setup(const Rcpp::List &x);
	void import_params(const Rcpp::DataFrame &x);
    void import_augm(const Rcpp::DataFrame &x);

    //===== Rcpp EXPORTS =====//
    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_params() const;
    Rcpp::DataFrame export_data() const;
    Rcpp::DataFrame export_traces() const;
    Rcpp::NumericMatrix export_loglike_matrix() const;
#else
    //===== Cpp PRINTS =====//
    void print_infer_state(const std::string &pathtofile) const;
    void print_augm_variables(const std::string &pathtofile) const;
#endif

    //===== READ =====//
    void read_setup(const std::string &pathtofile);
    void read_params(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);

    //===== SHOW =====//
    void show_results() const;
    void show_group_classes() const;
    void show_summary() const;
    void show_full_build() const;

    //===== SETUP =====//

    //----- Setup : Getter -----//
    bool get_rnd_init_params() const;
    ulong get_n_chains() const;
    ulong get_n_steps() const;
    ulong get_n_adapt() const;
    ulong get_thin() const;
    ulong get_n_burnin() const;
    double get_inv_T() const;
    unsigned get_n_threads() const;
    ulong get_c_n_sweeps() const;

    //===== SETUP : setter =====//
    void set_rnd_init_params(bool x);
    void set_inv_T(double x);

    void assign_house(const std::string &ID);
    void update_report_data();
    size_t report_data(const std::string &key) const;

    //===== AUGMENTED =====//
    void clear_all_augm(Household &house, indiv i);
    void assign_all_augm(Household &house, indiv i);

    sys_rnd_gen &get_prng();
    double uniform();

    void init_rand_params();
    double get_param(const std::string &parID);
    bool get_infer(const std::string &parID);
    const double *get_param_pointer(const std::string &parID) const;
    void assert_param(const std::string &parID) const;
    std::vector<std::string> get_param_to_be_inferred() const;

    bool propose_valid_param_move(const std::string &parID);
    void restore_param_value(const std::string &parID);
    void update_param_accept_rate(const std::string &parID, bool is_accepted);
    double param_move_logratio(const std::string &parID) const;
    void adapt_param_rule(const std::string &parID, double ratio);

    void disable_assert();
    void init_move_counts();
    void init_diagnostics(bool is_power_post);
    void update_diagnostics(double glob_loglike, bool is_power_post);
    void dgn_glob_ll_stats(double &avg, double &var) const;

	//----- Metropolis : augmented moves -----//
	double compute_augm_model_loglike(const Household &house) const;
	
    double (Monte_Carlo_Build::*_func_augm_model_loglike_)(const Household &) const = &Monte_Carlo_Build::_no_augm_model_loglike_;
	double _no_augm_model_loglike_(const Household &house) const;
	double _augm_model_loglike_(const Household &house) const;

    void (Monte_Carlo_Build::*_func_augm_model_by_definition_)(const Household &, indiv, double &) const = &Monte_Carlo_Build::_augm_model_by_status_;
	void _augm_model_by_status_(const Household &house, indiv i, double &delta) const;
	void _augm_model_generic_(const Household &house, indiv i, double &delta) const;

private:
    //===== Pointer to Core =====//
    std::unique_ptr<Monte_Carlo_Core> pCore;

    //===== BUILD =====//

    //----- build : auxiliary attributes -----//
    std::map<std::string, std::vector<std::string>> _augm_deftable_ = {};
    
    //----- build : type-specified variables -----//
    double type_to_double(const std::string tkn_val, const std::string tkn_typ, const Household &house, indiv i);
    const double *type_to_doubleptr(const std::string tkn_val, const std::string tkn_typ, const Household &house, indiv i);
    dictargs typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Household &house, indiv i);
    dictargs_ptr typedict_to_sdargs_ptr(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Household &house, indiv i);

    //----- build : Augmented Variables -----//
    void build_augm();
};

#endif // HOUSEPI__MCMC__BUILD__H
