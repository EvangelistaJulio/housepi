/*
 * build.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef MCMC__BUILD__H
#define MCMC__BUILD__H

#include <string>
#include <map>
#include <vector>
#include <memory>

#ifdef USE_RCPP
#include <Rcpp.h>
#endif

#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// Forward declaration of implementation class
class Monte_Carlo_Core;

/* Class to operate on data and augmented variables. */
class Monte_Carlo_Build
{
public:
    //===== Household Data =====//
    size_t num_indiv = 0;
    size_t num_house = 0;
    std::vector<Household> house = {};

    //===== Data Reports =====//
    std::vector<std::string> data_headers = {};
    ulong _num_contacts_ = 0;
    std::map<std::string, int> _report_data_;
    std::map<std::string, std::map<std::string, int>> _report_group_classes_ = {};
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

#ifdef USE_RCPP
	//===== Rcpp IMPORTS =====//
	void import_setup(Rcpp::List x);
	void import_params(Rcpp::DataFrame x);
    void import_augm(Rcpp::DataFrame x);

    //===== Rcpp EXPORTS =====//
    Rcpp::List export_setup() const;
    Rcpp::DataFrame export_params() const;
    Rcpp::DataFrame export_data() const;
    Rcpp::DataFrame export_traces() const;
    Rcpp::NumericMatrix export_loglike_matrix() const;
#else
    void show_results() const;
    void print_infer_state(const std::string &pathtofile) const;
    void print_augm_variables(const std::string &pathtofile) const;
#endif

    void read_setup(const std::string &pathtofile);
    void read_params(const std::string &pathtofile);
    void read_augm(const std::string &pathtofile);

    //===== Show =====//
    void show_group_classes() const;
    void show_summary() const;
    void show_full_build() const;

    //===== SETUP : getter =====//
    bool get_rnd_init_params() const;
    ulong get_num_chains() const;
    ulong get_num_iter() const;
    ulong get_num_adapt() const;
    ulong get_thinning() const;
    ulong get_burn_in() const;
    double get_inv_T() const;
    unsigned get_num_threads() const;
    ulong get_c_num_sweeps() const;

    //===== SETUP : setter =====//
    void set_rnd_init_params(bool x);
    void set_inv_T(double x);

    void assign_house(const std::string &ID);
    void update_report_data();
    int report_data(const std::string &key) const;

    void assign_augm(Individual &this_indiv);
    void clear_augm(Individual &this_indiv);

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
	double compute_augm_model_loglike(const Household &this_house) const { return (this->*_func_augm_model_loglike_)(this_house); }
	
    double (Monte_Carlo_Build::*_func_augm_model_loglike_)(const Household &) const = &Monte_Carlo_Build::_no_augm_model_loglike_;
	double _no_augm_model_loglike_(const Household &this_house) const;
	double _augm_model_loglike_(const Household &this_house) const;

    void (Monte_Carlo_Build::*_func_augm_model_by_definition_)(const Individual &, double &) const = &Monte_Carlo_Build::_augm_model_by_status_;
	void _augm_model_by_status_(const Individual &this_indiv, double &delta) const;
	void _augm_model_generic_(const Individual &this_indiv, double &delta) const;

private:
    //===== Pointer to Core =====//
    std::unique_ptr<Monte_Carlo_Core> pCore;

    //===== BUILD =====//

    //----- build : auxiliary attributes -----//
    std::map<std::string, std::vector<std::string>> _augm_deftable_ = {};
    
    //----- build : type-specified variables -----//
    double type_to_double(const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);
    const double *type_to_doubleptr(const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);
    dictargs typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);
    dictargs_ptr typedict_to_sdargs_ptr(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);

    //----- build : Augmented Variables -----//
    void build_augm();
};

#endif // MCMC__BUILD__H