/*
 * mcmc_R_import_export.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef _HOUSEPI_RCPP

/* Convert R List into map<string, string>. */
void rlist_to_map(const Rcpp::List &list,
                  std::map<std::string, std::string> &output)
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

/* Convert R dataframe into map<string, vector<string>>. */
void rdataframe_to_map(const Rcpp::DataFrame &df,
                       dict_str_vec &output)
{
    using namespace Rcpp;

    CharacterVector keys = df.names();

    output.clear();
    for (size_t i = 0; i < df.size(); ++i)
    {
        CharacterVector col_str = as<CharacterVector>(df[i]);
        str_vec col_data(col_str.begin(), col_str.end());
        output[as<std::string>(keys[i])] = std::move(col_data);
    }
}

Rcpp::String stoRString(const std::string &x)
{
    return (tools::is_na(x)) ? NA_STRING : Rcpp::String(x);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import setup from R list and run build. */
void MCMC::Impl::import_setup(const Rcpp::List &input)
{
    try
    {
        rlist_to_map(input, _setup_deflist_);
    }
    catch (const std::exception &e)
    {
        tools::error("While importing setup definitions: " + std::string(e.what()));
    }
    build_setup();
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
        _report_pathtofile_.erase("setup");
}

/* Import infectivity profile from R list and run build. */
void MCMC::Impl::import_infectivity_profile(const Rcpp::List &input)
{
    try
    {
        rlist_to_map(input, _infec_profile_deflist_);
    }
    catch (const std::exception &e)
    {
        tools::error("While importing infectivity profile definitions: " + std::string(e.what()));
    }
    build_infec_profile();
    if (_report_pathtofile_.find("infec_profile") != _report_pathtofile_.end())
        _report_pathtofile_.erase("infec_profile");
}

/* Import params from R DataFrame and run build. */
void MCMC::Impl::import_parameter(const Rcpp::DataFrame &input)
{
    try
    {
        rdataframe_to_map(input, _param_deftable_);
    }
    catch (const std::exception &e)
    {
        tools::error("While reading parameter definitions: " + std::string(e.what()));
    }
    build_params();
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
        _report_pathtofile_.erase("params");
}

/* Import data from R DataFrame then run build. */
void MCMC::Impl::import_augmented(const Rcpp::DataFrame &input)
{
    try
    {
        rdataframe_to_map(input, _augm_deftable_);
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While reading augmented definitions: " + std::string(e.what()));
    }
    build_augm();
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
        _report_pathtofile_.erase("augm");
}

/* Import data from R DataFrame and run build. */
void MCMC::Impl::import_dataset(const Rcpp::DataFrame &data)
{
    dict_str_vec datatable = {};
    rdataframe_to_map(data, datatable);
    build_data(datatable);
    if (_report_pathtofile_.find("data") != _report_pathtofile_.end())
        _report_pathtofile_.erase("data");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current setup configurations to R List. */
Rcpp::List MCMC::Impl::export_setup() const
{
    using namespace Rcpp;

    return List::create(_["n_chains"] = setup_n_chains,
                        _["sample_size"] = setup_c_sample_size,
                        _["n_steps"] = setup_n_steps,
                        _["n_adapt"] = setup_n_adapt,
                        _["n_burnin"] = setup_n_burnin,
                        _["thin"] = setup_thin,
                        _["rand_init"] = setup_rand_init_params,
                        _["seed"] = setup_seed,
                        _["n_threads"] = setup_n_threads,
                        _["id_trace"] = setup_id_trace,
                        _["progress"] = setup_is_progbar,
                        _["inv_temper"] = setup_inv_T);

}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_parameter() const
{
    using namespace Rcpp;

    std::vector<parID> V_names = param_standard_infer_list + param_extra_infer_list + param_const_list;

    size_t n_rows = V_names.size();

    // check if not empty
    if (V_names.empty())
        tools::error("No parameters to be exported.");

    // pre-allocate memory
    NumericVector V_value(n_rows), V_min(n_rows), V_max(n_rows);
    LogicalVector V_constant(n_rows), V_discrete(n_rows);
    CharacterVector V_prior(n_rows), V_arg_name(n_rows), V_arg_value(n_rows);
    CharacterVector V_rule(n_rows);
    NumericVector V_rule_scale(n_rows);

    // fill datafrane columns
    for (size_t c = 0; c < n_rows; ++c)
    {
        parID p = V_names[c];
        V_value[c] = get_param(p);
        V_constant[c] = param_map.at(p).is_constant();
        V_discrete[c] = param_map.at(p).is_discrete();

        V_prior[c] = stoRString(param_map.at(p).show_distr_name());
        V_arg_name[c] = stoRString(param_map.at(p).show_distr_arg_names());
        V_arg_value[c] = stoRString(param_map.at(p).show_distr_arg_values());

        V_rule[c] = stoRString(param_map.at(p).show_prop_name());

        const std::string rule_scal = param_map.at(p).show_prop_scale();
        V_rule_scale[c] = tools::is_na(rule_scal) ? NA_REAL : tools::stod(rule_scal);

        const std::string min_val = param_map.at(p).show_min();
        V_min[c] = tools::is_na(min_val) ? NA_REAL : tools::stod(min_val);

        const std::string max_val = param_map.at(p).show_max();
        V_max[c] = tools::is_na(max_val) ? NA_REAL : tools::stod(max_val);
    }

    return DataFrame::create(_["name"] = wrap(V_names), _["value"] = V_value,
                             _["constant"] = V_constant, _["discrete"] = V_discrete,
                             _["prior"] = V_prior, _["arg_name"] = V_arg_name, _["arg_value"] = V_arg_value,
                             _["rule"] = V_rule, _["rule_scale"] = V_rule_scale,
                             _["min"] = V_min, _["max"] = V_max);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Monte_Carlo_Build : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_dataset() const
{
    using namespace Rcpp;

    CharacterVector df_index(n_individuals);
    CharacterVector df_house_index(n_individuals);
    IntegerVector df_status(n_individuals);
    std::vector<NumericVector> df_dbl;
    for (size_t s = 0; s < data_headers.size() + augm_list_full.size(); ++s)
        df_dbl.emplace_back(n_individuals);

    size_t r = 0, c = 0;
    for (const Household &house : house_vector)
        for (indiv i = 0; i < house.size(); ++i)
        {
            df_index[r] = house.id(i);
            df_house_index[r] = house.house_id();
            df_status[r] = house.status(i);

            c = 0;
            for (const dataID &d : data_headers)
                df_dbl[c++][r] = house.data(i, d);

            for (const augmID &a : augm_list_full)
                df_dbl[c++][r] = house.augm(i, a);

            ++r;
        }

    DataFrame output = DataFrame::create(_["id"] = df_index,
                                         _["house_id"] = df_house_index,
                                         _["status"] = df_status);
    c = 0;
    for (const dataID &d : data_headers)
        output.push_back(df_dbl[c++], d);

    for (const augmID &a : augm_list_full)
        output.push_back(df_dbl[c++], a);

    return output;
}

/* Export MCMC traces of parameter posterior samples. */
Rcpp::NumericMatrix MCMC::Impl::export_traces() const
{
    using namespace Rcpp;

    if ((setup_n_chains == 0) || (is_power_post()))
        tools::error(
            "No traces to export were found.");

    //----- Collect the column names -----//
    str_vec headers = {"loglikelihood"};
    headers += param_standard_infer_list + param_extra_infer_list;

    // augmented variable of id trace
    if (!dgn_augm_value.empty())
        for (const augmID &a : augm_infer_list)
            headers += {a + "_id_" + setup_id_trace};

    // multiple chains
    if (setup_n_chains > 1)
        headers += {"CHAIN__"};

    //----- NumericMatrix with traces -----//
    NumericMatrix output(setup_n_chains * setup_c_sample_size,
                         headers.size());
    size_t col = 0;

    // loglikelihood
    output(_, col++) = NumericVector(wrap(dgn_glob_loglike));

    // parameters
    for (const parID &p : param_standard_infer_list + param_extra_infer_list)
        output(_, col++) = NumericVector(wrap(dgn_param_value.at(p)));

    // augmented variables
    if (!dgn_augm_value.empty())
        for (const augmID &a : augm_infer_list)
            output(_, col++) = NumericVector(wrap(dgn_augm_value.at(a)));

    // retrieve chain identities
    if (setup_n_chains > 1)
    {
        IntegerVector df_chain(setup_c_sample_size * setup_n_chains);
        for (size_t c = 0; c < setup_n_chains; ++c)
            std::fill(df_chain.begin() + c * setup_c_sample_size,
                      df_chain.begin() + (c + 1) * setup_c_sample_size, c + 1);
        output(_, col++) = df_chain;
    }

    //----- associate column names -----//
    colnames(output) = CharacterVector(wrap(headers));

    //----- associate row names -----//
    IntegerVector sweep(setup_c_sample_size * setup_n_chains);
    for (size_t c = 0; c < setup_n_chains; ++c)
        for (size_t k = 0; k < setup_c_sample_size; k++)
            sweep[k + c * setup_c_sample_size] = setup_n_burnin + (k + 1) * setup_thin;
    rownames(output) = CharacterVector(wrap(sweep));

    return output;
}

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix MCMC::Impl::export_household_loglike_matrix() const
{
    using namespace Rcpp;

    if ((setup_n_chains == 0) || (is_power_post()))
        tools::error(
            "No traces to export were found.");
            
    // Create a NumericMatrix with n_house columns
    ulong sample_size = setup_n_chains * setup_c_sample_size;
    NumericMatrix output(sample_size, n_households);
    CharacterVector headers(n_households);

    // Assign each NumericVector to a column of the matrix
    for (size_t H = 0; H < n_households; ++H)
    {
        output(_, H) = NumericVector(wrap(house_vector[H].get_dgn_loglike()));
        headers[H] = house_vector[H].house_id();
    }

    //----- associate column names -----//
    colnames(output) = headers;
    
    //----- associate row names -----//
    IntegerVector sweep(setup_c_sample_size * setup_n_chains);
    for (size_t c = 0; c < setup_n_chains; ++c)
        for (size_t k = 0; k < setup_c_sample_size; k++)
            sweep[k + c * setup_c_sample_size] = setup_n_burnin + (k + 1) * setup_thin;
    rownames(output) = CharacterVector(wrap(sweep));
    
    return output;
}

/* Export a showcase of the selected infectivity profile. */
Rcpp::DataFrame MCMC::Impl::export_infectivity_profile(const Rcpp::List &opts,
                                                 const Rcpp::List &xlim,
                                                 int n_bins) const
{
    using namespace Rcpp;

    // check the call makes sense
    if (_infec_profile_ == nullptr)
        tools::print_stderr("No infectivity profile to export.");

    // load options
    dbl_vec opts_ = {};
    for (int k = 0; k < opts.size(); ++k)
        opts_.push_back(opts.at(k));

    // load x limits
    dbl_vec xlim_ = {-20., 20.};
    if (xlim.size() > 0)
    {
        std::map<std::string, std::string> load_user_xlim = {};
        rlist_to_map(xlim, load_user_xlim);
        if (load_user_xlim.count("min"))
            if (!tools::is_na(load_user_xlim.at("min")))
                xlim_[0] = tools::stod(load_user_xlim.at("min"));
        if (load_user_xlim.count("max"))
            if (!tools::is_na(load_user_xlim.at("max")))
                xlim_[1] = tools::stod(load_user_xlim.at("max"));
    }

    // produce a showcase of the distribution profile
    std::map<double, dbl_vec> values = _infec_profile_->showcase(opts_, xlim_, n_bins);
    size_t c = 0, n = values.size();
    NumericVector t(n), inst(n), cumu(n);

    for (const auto &pair : values)
    {
        t[c] = pair.first;
        inst[c] = pair.second[0];
        cumu[c] = pair.second[1];
        ++c;
    }

    return DataFrame::create(_["generation_time"] = t,
                             _["instantaneous"] = inst,
                             _["cumulative"] = cumu);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : exports to Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::List MCMC::export_setup() const { return pImpl->export_setup(); }
Rcpp::DataFrame MCMC::export_parameter() const { return pImpl->export_parameter(); }
Rcpp::DataFrame MCMC::export_dataset() const { return pImpl->export_dataset(); }
Rcpp::NumericMatrix MCMC::export_traces() const { return pImpl->export_traces(); }
Rcpp::NumericMatrix MCMC::export_household_loglike_matrix() const { return pImpl->export_household_loglike_matrix(); }
Rcpp::DataFrame MCMC::export_infectivity_profile(const Rcpp::List &opts,
                                           const Rcpp::List &xlim,
                                           int n_bins)
{
    init_infec_profile();
    return pImpl->export_infectivity_profile(opts, xlim, n_bins);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(const Rcpp::List &x) { pImpl->import_setup(x); }
void MCMC::import_infectivity_profile(const Rcpp::List &x) { pImpl->import_infectivity_profile(x); }
void MCMC::import_parameter(const Rcpp::DataFrame &x) { pImpl->import_parameter(x); }
void MCMC::import_augmented(const Rcpp::DataFrame &x) { pImpl->import_augmented(x); }
void MCMC::import_dataset(const Rcpp::DataFrame &x) { pImpl->import_dataset(x); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#else

#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

/* Read Setup file and produce a definition list. */
void MCMC::Impl::read_setup(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_list(pathtofile, _setup_deflist_);
    }
    catch (const std::exception &e)
    {
        tools::error("While reading setup definitions: " + std::string(e.what()));
    }
    build_setup();
    _report_pathtofile_["setup"] = pathtofile;
}

/* Read Setup file and produce a definition list. */
void MCMC::Impl::read_infec_profile(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_list(pathtofile, _infec_profile_deflist_);
    }
    catch (const std::exception &e)
    {
        tools::error("While reading infectivity profile definitions: " + std::string(e.what()));
    }
    build_infec_profile();
    _report_pathtofile_["infec_profile"] = pathtofile;
}

/* Read parameter definitions and run assignment. */
void MCMC::Impl::read_params(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_table(pathtofile, _param_deftable_);
    }
    catch (const std::exception &e)
    {
        tools::error("While reading parameter definitions: " + std::string(e.what()));
    }
    build_params();
    _report_pathtofile_["params"] = pathtofile;
}

/* Read Model Augmented Variable file. */
void MCMC::Impl::read_augm(const std::string &pathtofile)
{
    try
    {
        tools::read_csv_table(pathtofile, _augm_deftable_);
    }
    catch (const std::exception &e)
    {
        tools::error(
            "While reading augmented definitions: " + std::string(e.what()));
    }
    build_augm();
    _report_pathtofile_["augm"] = pathtofile;
}

/* Read data file and then run assignment. */
void MCMC::Impl::read_data(const std::string &pathtofile)
{
    dict_str_vec datatable = {};
    try
    {
        tools::read_csv_table(pathtofile, datatable);
    }
    catch (const std::exception &e)
    {
        tools::error("While reading dataset: " + std::string(e.what()));
    }
    build_data(datatable);
    _report_pathtofile_["data"] = pathtofile;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show summary  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary function to call percentuage. */
inline const std::string _aux_percent_(int num, int den)
{
    return "(" + tools::percentage(static_cast<double>(num) / den, 1) + ")";
}

/* Print the report of the database and of the model group rates. */
void MCMC::Impl::show_summary() const
{
    std::ostringstream oss;

    oss << "n. of individuals: " << n_individuals
        << "\n    * n. of cases: " << _report_data_.at("n.cases")
        << " " << _aux_percent_(_report_data_.at("n.cases"), n_individuals);

    if (_report_data_.at("n.miss.stat") > 0)
        oss << "\n    * n. of missing status : " << _report_data_.at("n.miss.stat")
            << " " << _aux_percent_(_report_data_.at("n.miss.stat"), n_individuals);

    oss << "\nn. of households: " << n_households;

    if (_report_data_.at("no.HH.case") > 0)
        oss << "\n    * with no cases: " << _report_data_.at("no.HH.case")
            << " " << _aux_percent_(_report_data_.at("no.HH.case"), n_households);

    oss << "\n    * with at least 1 case: " << _report_data_.at("lst.1.case")
        << " " << _aux_percent_(_report_data_.at("lst.1.case"), n_households);

    if (_report_data_.at("lst.1.case") > 0)
        oss << "\n    * with at least 2 cases: " << _report_data_.at("lst.2.case")
            << " " << _aux_percent_(_report_data_.at("lst.2.case"), n_households);

    oss << "\nmax n. house cases: " << _report_data_.at("max.n.HH.cases")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.max.cases"), n_households)
        << "\nmin n. house members: " << _report_data_.at("min.HH.size")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.min.size"), n_households)
        << "\nmax n. house members: " << _report_data_.at("max.HH.size")
        << " " << _aux_percent_(_report_data_.at("n.HH.with.max.size"), n_households);

    tools::print_stdout(oss.str() + "\n");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : show setup  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the complete choices for the setup. */
void MCMC::Impl::show_setup() const
{
    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << setup_seed
        << "\nn_threads" << C_DELIM << setup_n_threads
        << "\nn_chains" << C_DELIM << setup_n_chains
        << "\nsample_size" << C_DELIM << setup_c_sample_size
        << "\nn_steps" << C_DELIM << setup_n_steps
        << "\nn_adapt" << C_DELIM << setup_n_adapt
        << "\nn_burnin" << C_DELIM << setup_n_burnin
        << "\nthin" << C_DELIM << setup_thin
        << "\nrand_init" << C_DELIM << tools::truth(setup_rand_init_params)
        << "\nprogress" << C_DELIM << tools::truth(setup_is_progbar)
        << "\nid_trace" << C_DELIM << setup_id_trace
        << "\ninv_temper" << C_DELIM << setup_inv_T
        << "\n";

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void MCMC::Impl::show_all_params() const
{
    std::ostringstream oss;
    oss << param_file_headers[0];

    for (size_t k = 1; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];

    std::vector<parID> all_params = param_standard_infer_list + param_extra_infer_list + param_const_list;

    for (const parID &p : all_params)
    {
        const Variable &aux = param_map.at(p);
        oss << aux.show_name()
            << C_DELIM << tools::natod(aux.get_value())
            << C_DELIM << tools::truth(aux.is_constant())
            << C_DELIM << tools::truth(aux.is_discrete())
            << C_DELIM << aux.show_distr_name()
            << C_DELIM << aux.show_distr_arg_names()
            << C_DELIM << aux.show_distr_arg_values()
            << C_DELIM << aux.show_prop_name()
            << C_DELIM << tools::natod(aux.show_prop_scale())
            << C_DELIM << tools::natod(aux.show_min())
            << C_DELIM << tools::natod(aux.show_max())
            << "\n";
    }

    tools::print_stdout(oss.str());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : show full build  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::Impl::show_full_build() const
{
    //----- SETUP -----//
    std::string add_file_info = NULL_STRING;
    if (_report_pathtofile_.count("setup"))
        if (!tools::is_na(_report_pathtofile_.at("setup")))
            add_file_info = "# import: " + _report_pathtofile_.at("setup");

    tools::print_decorator("Setup", add_file_info);
    show_setup();

    //----- PARAMETERS -----//
    add_file_info = NULL_STRING;
    if (_report_pathtofile_.count("params"))
        if (!tools::is_na(_report_pathtofile_.at("params")))
            add_file_info = "# import: " + _report_pathtofile_.at("params");

    tools::print_decorator("Parameters", add_file_info);
    show_all_params();

    //----- DATABASE -----//
    add_file_info = NULL_STRING;
    if (_report_pathtofile_.count("data"))
        if (!tools::is_na(_report_pathtofile_.at("data")))
            add_file_info = "# import: " + _report_pathtofile_.at("data");

    tools::print_decorator("Database", add_file_info);
    show_summary();

    //----- AUGMENTED -----//
    add_file_info = NULL_STRING;
    if (_report_pathtofile_.count("augm"))
        if (!tools::is_na(_report_pathtofile_.at("augm")))
            add_file_info = "# import: " + _report_pathtofile_.at("augm");

    if (!augm_list_full.empty())
        tools::print_decorator("Augmented", add_file_info);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void MCMC::Impl::show_results() const
{
    std::ostringstream oss;
    oss << "name" << C_DELIM << "value" << C_DELIM << "value_sd";

    for (size_t k = 2; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];

    oss << "\n";

    for (const parID &p : param_standard_infer_list + param_extra_infer_list)
    {
        const Variable &this_param = param_map.at(p);
        double avg = 0, var = 0;
        tools::compute_stats(dgn_param_value.at(p), avg, var);

        oss << p
            << C_DELIM << avg
            << C_DELIM << std::sqrt(var)
            << C_DELIM << tools::truth(this_param.is_constant())
            << C_DELIM << tools::truth(this_param.is_discrete())
            << C_DELIM << this_param.show_distr_name()
            << C_DELIM << this_param.show_distr_arg_names()
            << C_DELIM << this_param.show_distr_arg_values()
            << C_DELIM << this_param.show_prop_name()
            << C_DELIM << this_param.show_prop_scale()
            << C_DELIM << tools::natod(this_param.show_min())
            << C_DELIM << tools::natod(this_param.show_max())
            << "\n";
    }

    tools::print_stdout(oss.str());
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the output file with inferred parameters with column names and initial values */
void MCMC::Impl::print_infer_state(const std::string &output_filename) const
{
    if (!tools::is_na(output_filename))
    {
        std::ofstream inferstate_ofstream;
        inferstate_ofstream.open(output_filename.c_str());

        if (inferstate_ofstream.is_open())
        {
            //----- HEADERS -----//
            // header loglikelihood
            inferstate_ofstream << "iteration"
                                << C_DELIM << "loglikelihood";

            // header parameters
            for (const parID &p : param_standard_infer_list + param_extra_infer_list)
                inferstate_ofstream << C_DELIM << p;

            // augmented variables
            if (is_augm_sample() && (!tools::is_na(setup_id_trace)))
                for (const augmID &a : augm_infer_list)
                    inferstate_ofstream << C_DELIM << a + "_id_" + setup_id_trace;

            // header chain
            if (setup_n_chains > 1)
                inferstate_ofstream << C_DELIM << "CHAIN__";

            inferstate_ofstream << std::endl;

            //----- VALUES -----//

            for (size_t c = 0; c < setup_n_chains; ++c)
                for (size_t k = 0; k < setup_c_sample_size; k++)
                {
                    size_t i = k + c * setup_c_sample_size;

                    // print iteration
                    ulong sweep = setup_n_burnin + (k + 1) * setup_thin;
                    inferstate_ofstream << sweep;

                    // print loglikelihood
                    inferstate_ofstream << C_DELIM << std::setprecision(2)
                                        << std::fixed << dgn_glob_loglike[i];

                    // print parameter values
                    for (const parID &p : param_standard_infer_list + param_extra_infer_list)
                        inferstate_ofstream << C_DELIM << std::setprecision(5)
                                            << std::fixed << dgn_param_value.at(p)[i];

                    // print selected augmented variables
                    if (is_augm_sample() && (!tools::is_na(setup_id_trace)))
                        for (const augmID &a : augm_infer_list)
                            inferstate_ofstream << C_DELIM << std::setprecision(5)
                                                << std::fixed << dgn_augm_value.at(a)[i];

                    // print which chain
                    if (setup_n_chains > 1)
                        inferstate_ofstream << C_DELIM << (c + 1);

                    inferstate_ofstream << std::endl;
                }

            inferstate_ofstream.close();
        }
        else
            throw std::runtime_error(
                "Unable to open file: " + output_filename);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC::Impl : print augm variables  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Print all augmented variables of each individual */
void MCMC::Impl::print_augm_variables(const std::string &pathtofile) const
{
    if (!tools::is_na(pathtofile))
    {
        std::ofstream augm_ofstream;
        augm_ofstream.open(pathtofile.c_str());

        if (augm_ofstream.is_open())
        {
            //===== Print column headers =====//
            // default headers
            augm_ofstream << "id" << C_DELIM << "house_id" << C_DELIM << "status";

            // data headers
            for (const dataID &d : data_headers)
                augm_ofstream << C_DELIM << d;

            // augmented variables headers
            for (const augmID &a : augm_list_full)
                augm_ofstream << C_DELIM << a;

            augm_ofstream << std::endl;

            //===== Print data and augmented variables =====//
            for (const Household &house : house_vector)
                for (indiv i = 0; i < house.size(); ++i)
                {
                    augm_ofstream << house.id(i)
                                  << C_DELIM << house.house_id()
                                  << C_DELIM << house.status(i);
                    for (const dataID &d : data_headers)
                        augm_ofstream << C_DELIM << house.data(i, d);

                    for (const augmID &a : augm_list_full)
                        if (house.status(i) == INFECTED_STATUS)
                            augm_ofstream << C_DELIM << house.augm(i, a);
                        else
                            augm_ofstream << C_DELIM << std::numeric_limits<double>::quiet_NaN();

                    augm_ofstream << std::endl;
                }

            augm_ofstream.close();
        }
        else
        {
            throw std::runtime_error(
                "Unable to open file: " + pathtofile);
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : read  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::read_setup(const std::string &pathtofile) { pImpl->read_setup(pathtofile); }
void MCMC::read_params(const std::string &pathtofile) { pImpl->read_params(pathtofile); }
void MCMC::read_augm(const std::string &pathtofile) { pImpl->read_augm(pathtofile); }
void MCMC::read_data(const std::string &pathtofile) { pImpl->read_data(pathtofile); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : show  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::show_full_build() const { pImpl->show_full_build(); }
void MCMC::show_results() const { pImpl->show_results(); }
void MCMC::show_summary() const { pImpl->show_summary(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MCMC : C++ Print  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::print_infer_state(const std::string &pathtofile) const { pImpl->print_infer_state(pathtofile); }
void MCMC::print_augm_variables(const std::string &pathtofile) const { pImpl->print_augm_variables(pathtofile); }

#endif
