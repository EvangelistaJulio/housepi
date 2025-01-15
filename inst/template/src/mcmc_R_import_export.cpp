/*
 * mcmc_R_import_export.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc_build.h"
#include "mcmc.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef _HOUSEPI_RCPP

/* Convert R List into map<string, string>. */
void rlist_to_map(const Rcpp::List &list, std::map<std::string, std::string> &output)
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
void rdataframe_to_map(const Rcpp::DataFrame &df, dict_str_vec &output)
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

/* Import setup from R list and run build. */
void Monte_Carlo_Build::import_setup(const Rcpp::List &setup)
{
    rlist_to_map(setup, _setup_deflist_);
    build_setup();
    if (_report_pathtofile_.find("setup") != _report_pathtofile_.end())
        _report_pathtofile_.erase("setup");
}

/* Import params from R DataFrame and run build. */
void Monte_Carlo_Build::import_params(const Rcpp::DataFrame &params)
{
    rdataframe_to_map(params, _param_deftable_);
    build_params();
    if (_report_pathtofile_.find("params") != _report_pathtofile_.end())
        _report_pathtofile_.erase("params");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  Monte_Carlo_Build : imports from Rcpp ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Import data from R DataFrame then run build. */
void Monte_Carlo_Build::import_augm(const Rcpp::DataFrame &augm)
{
    rdataframe_to_map(augm, _augm_deftable_);
    build_augm();
    if (_report_pathtofile_.find("augm") != _report_pathtofile_.end())
        _report_pathtofile_.erase("augm");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current setup configurations to R List. */
Rcpp::List Monte_Carlo_Build::export_setup() const
{
    using namespace Rcpp;

    return List::create(
        _["seed"] = setup_seed,
        _["n_threads"] = setup_n_threads,
        _["n_chains"] = setup_n_chains,
        _["n_adapt"] = setup_n_adapt,
        _["n_burnin"] = setup_n_burnin,
        _["n_steps"] = setup_n_steps,
        _["thin"] = setup_thin,
        _["rand_init"] = setup_rand_init_params,
        _["inv_temper"] = setup_inv_T);
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame Monte_Carlo_Build::export_params() const
{
    using namespace Rcpp;

    std::vector<parID> V_names = {};

    if (!param_to_be_inferred.empty())
        V_names.insert(V_names.end(), param_to_be_inferred.begin(), param_to_be_inferred.end());

    if (!param_to_be_constant.empty())
        V_names.insert(V_names.end(), param_to_be_constant.begin(), param_to_be_constant.end());

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

        V_prior[c] = stoRString(param_map.at(p).get_distr_name());
        V_arg_name[c] = stoRString(param_map.at(p).show_distr_arg_names());
        V_arg_value[c] = stoRString(param_map.at(p).show_distr_arg_values());

        V_rule[c] = stoRString(param_map.at(p).get_prop_name());

        const std::string rule_scal = param_map.at(p).show_prop_scale();
        V_rule_scale[c] = tools::is_na(rule_scal) ? NA_REAL : tools::stod(rule_scal);

        const std::string min_val = param_map.at(p).show_range_min();
        V_min[c] = tools::is_na(min_val) ? NA_REAL : tools::stod(min_val);

        const std::string max_val = param_map.at(p).show_range_max();
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
Rcpp::DataFrame Monte_Carlo_Build::export_data() const
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
    {
        for (indiv i = 0; i < house.size(); ++i)
        {
            df_index[r] = house.ID(i);
            df_house_index[r] = house.house_ID();
            df_status[r] = house.status(i);

            c = 0;
            for (const dataID &d : data_headers)
                df_dbl[c++][r] = house.data(i, d);

            for (const augmID &a : augm_list_full)
                df_dbl[c++][r] = house.augm(i, a);

            ++r;
        }
    }

    DataFrame output;
    output = DataFrame::create(_["id"] = df_index,
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
Rcpp::DataFrame Monte_Carlo_Build::export_traces() const
{
    using namespace Rcpp;

    DataFrame output = DataFrame::create(
        _["loglikelihood"] = wrap(dgn_glob_loglike));

    // parameters
    for (const parID &p : param_to_be_inferred)
        output.push_back(wrap(dgn_param_value.at(p)), p);

    // augmented variables
    if (is_augm_infer() && (!tools::is_na(setup_id_trace)))
        for (const augmID &a : augm_list_to_be_inferred)
            output.push_back(wrap(dgn_augm_value.at(a)), a);

    // retrieve chain identities
    if (setup_n_chains > 1)
    {
        IntegerVector df_chain(setup_c_sample_size * setup_n_chains);
        for (ulong c = 0; c < setup_n_chains; ++c)
            std::fill(df_chain.begin() + c * setup_c_sample_size,
                      df_chain.begin() + (c + 1) * setup_c_sample_size, c + 1);

        output.push_back(df_chain, "CHAIN__");
    }

    return output;
}

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix Monte_Carlo_Build::export_household_loglike_traces() const
{
    using namespace Rcpp;

    // FIXME : this should be prevented if is_power_post

    // Create a NumericMatrix with n_house columns
    ulong sample_size = setup_n_chains * setup_c_sample_size;
    NumericMatrix output(sample_size, n_households);

    // Assign each NumericVector to a column of the matrix
    for (size_t H = 0; H < n_households; ++H)
    {
        NumericVector column = wrap(house_vector[H].get_dgn_loglike());
        output(_, H) = column;
    }

    return output;
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

Rcpp::List MCMC::export_setup() const { return pBuild->export_setup(); }
Rcpp::DataFrame MCMC::export_params() const { return pBuild->export_params(); }
Rcpp::DataFrame MCMC::export_data() const { return pBuild->export_data(); }
Rcpp::DataFrame MCMC::export_traces() const { return pBuild->export_traces(); }
Rcpp::NumericMatrix MCMC::export_household_loglike_traces() const { return pBuild->export_household_loglike_traces(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(const Rcpp::List &x) { pBuild->import_setup(x); }
void MCMC::import_params(const Rcpp::DataFrame &x) { pBuild->import_params(x); }
void MCMC::import_augm(const Rcpp::DataFrame &x) { pBuild->import_augm(x); }

/* Import data from R DataFrame and run build. */
void MCMC::import_data(const Rcpp::DataFrame &data)
{
    dict_str_vec datatable = {};
    rdataframe_to_map(data, datatable);
    build_data(datatable);
    if (pBuild->_report_pathtofile_.find("data") != pBuild->_report_pathtofile_.end())
        pBuild->_report_pathtofile_.erase("data");
}

#endif
