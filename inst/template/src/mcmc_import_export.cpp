/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_import_export.cpp
 *
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
    for (R_xlen_t i = 0; i < list.size(); ++i)
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
    for (R_xlen_t i = 0; i < df.size(); ++i)
    {
        CharacterVector col_str = as<CharacterVector>(df[i]);
        str_vec col_data(col_str.begin(), col_str.end());
        output[as<std::string>(keys[i])] = std::move(col_data);
    }
}

/* Format C++ string into Rcpp string. */
Rcpp::String stoRString(const std::string &x)
{
    if (tools::is_na(x))
        return HOUSEPI_NA_STRING;
    return Rcpp::String(x);
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
        rlist_to_map(input, _setup_deflist);
    }
    catch (const std::exception &e)
    {
        tools::error({"While importing setup definitions:",
                      std::string(e.what())});
    }
    build_setup();
}

/* Import infectivity profile from R list and run build. */
void MCMC::Impl::import_infec_profile(const Rcpp::List &input)
{
    try
    {
        rlist_to_map(input, _infec_profile.deflist);
    }
    catch (const std::exception &e)
    {
        tools::error({"While importing infectivity profile definitions:",
                      std::string(e.what())});
    }
    _infec_profile.build();
}

/* Import params from R DataFrame and run build. */
void MCMC::Impl::import_param(const Rcpp::DataFrame &input)
{
    try
    {
        rdataframe_to_map(input, _param_deftable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading parameter definitions:",
                      std::string(e.what())});
    }
    build_param();
}

/* Import data from R DataFrame then run build. */
void MCMC::Impl::import_augm(const Rcpp::DataFrame &input)
{
    try
    {
        rdataframe_to_map(input, _augm_deftable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading augmented definitions:",
                      std::string(e.what())});
    }
    build_augm();
}

/* Import data from R DataFrame and run build. */
void MCMC::Impl::import_dataset(const Rcpp::DataFrame &input)
{
    dict_str_vec datatable = {};
    try
    {
        rdataframe_to_map(input, datatable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading dataset:",
                      std::string(e.what())});
    }
    build_dataset(datatable);
}

/* Import multiplicity definitions from R DataFrame and run build. */
void MCMC::Impl::import_multiplicity(const Rcpp::DataFrame &input)
{
    dict_str_vec multipltable = {};
    try
    {
        rdataframe_to_map(input, multipltable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading multiplicity:",
                      std::string(e.what())});
    }
    build_multiplicity(multipltable);
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
                        _["progress"] = setup_progress,
                        _["inv_temper"] = setup_inv_T);
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_param() const
{
    using namespace Rcpp;

    try
    {
        std::vector<parID> V_names = param_standard_infer_list + param_extra_infer_list + param_const_list;

        // check if not empty
        if (V_names.empty())
            throw std::length_error("No parameters to be exported.");

        // pre-allocate memory
        size_t n_rows = V_names.size();
        NumericVector V_value(n_rows), V_min(n_rows), V_max(n_rows);
        LogicalVector V_constant(n_rows), V_discrete(n_rows);
        CharacterVector V_prior(n_rows), V_arg_name(n_rows), V_arg_value(n_rows);
        CharacterVector V_rule(n_rows);
        NumericVector V_rule_scale(n_rows);

        // fill datafrane columns
        for (size_t r = 0; r < n_rows; ++r)
        {
            parID p = V_names[r];
            V_value[r] = get_param(p);
            V_constant[r] = param_map.at(p).is_constant();
            V_discrete[r] = param_map.at(p).is_discrete();

            V_prior[r] = stoRString(param_map.at(p).show_distr_name());
            V_arg_name[r] = stoRString(param_map.at(p).show_distr_arg_names());
            V_arg_value[r] = stoRString(param_map.at(p).show_distr_arg_values());

            V_rule[r] = stoRString(param_map.at(p).show_prop_name());

            const std::string rule_scale = param_map.at(p).show_prop_scale();
            V_rule_scale[r] = tools::is_na(rule_scale) ? NA_REAL : tools::stod(rule_scale);

            const std::string min_val = param_map.at(p).show_min();
            V_min[r] = tools::is_na(min_val) ? NA_REAL : tools::stod(min_val);

            const std::string max_val = param_map.at(p).show_max();
            V_max[r] = tools::is_na(max_val) ? NA_REAL : tools::stod(max_val);
        }

        return DataFrame::create(_["name"] = std::move(wrap(V_names)),
                                 _["value"] = std::move(V_value),
                                 _["const"] = std::move(V_constant),
                                 _["discr"] = std::move(V_discrete),
                                 _["prior"] = std::move(V_prior),
                                 _["arg_key"] = std::move(V_arg_name),
                                 _["arg_val"] = std::move(V_arg_value),
                                 _["min"] = std::move(V_min),
                                 _["max"] = std::move(V_max),
                                 _["rule"] = std::move(V_rule),
                                 _["scale"] = std::move(V_rule_scale));
    }
    catch (const std::exception &e)
    {
        tools::error({"While exporting parameter definitions:",
                      std::string(e.what())});
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Build : dataset exports  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_dataset()
{
    using namespace Rcpp;

    try
    {
        if (n_individuals == 0)
            throw std::length_error("No dataset to export.");

        size_t n_rows = n_individuals;
        size_t r = 0;
        Rcpp::CharacterVector df_index = CharacterVector(n_rows);
        Rcpp::CharacterVector df_house_index = CharacterVector(n_rows);
        Rcpp::IntegerVector df_status = IntegerVector(n_rows);
        std::vector<Rcpp::NumericVector> df_dbl = {};

        // other data and augmetned variable columns
        for (size_t s = 0; s < data_headers.size(); ++s)
            df_dbl.emplace_back(n_rows);

        for (const Household &house : house_vector)
            for (indiv i = 0; i < house.size(); ++i)
            {
                df_index[r] = house.id(i);
                df_house_index[r] = house.house_id();
                df_status[r] = house.status(i);

                size_t col = 0;
                for (const dataID &d : data_headers)
                {
                    double x = house.data(i, d);
                    df_dbl[col++][r] = tools::is_na(x) ? NA_REAL : x;
                }

                ++r;
            }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(
            _["id"] = df_index,
            _["house_id"] = df_house_index,
            _["status"] = df_status);

        size_t col = 0;
        for (const dataID &d : data_headers)
            output.push_back(df_dbl[col++], d);

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error({"While exporting dataset:",
                      std::string(e.what())});
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Build : report exports  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_data_report() const
{
    using namespace Rcpp;

    try
    {
        if (n_individuals == 0)
            throw std::length_error(
                "No report available for empty dataset.");

        size_t n_rows = _report_indiv_.size() + _report_house_.size();
        CharacterVector type(n_rows), comment(n_rows);
        IntegerVector counts(n_rows);
        NumericVector pct(n_rows);

        size_t k = 0;
        for (const auto &pair : _report_indiv_)
        {
            type[k] = "individual";
            comment[k] = pair.first;
            counts[k] = pair.second;
            pct[k] = 100. * static_cast<double>(pair.second) / n_individuals;
            ++k;
        }
        for (const auto &pair : _report_house_)
        {
            type[k] = "household";
            comment[k] = pair.first;
            counts[k] = pair.second;
            pct[k] = 100. * static_cast<double>(pair.second) / n_households;
            ++k;
        }

        return DataFrame::create(_["type"] = std::move(type),
                                 _["comment"] = std::move(comment),
                                 _["counts"] = std::move(counts),
                                 _["pct"] = std::move(pct));
    }
    catch (const std::exception &e)
    {
        tools::error({"While exporting dataset report:",
                      std::string(e.what())});
    }
    return DataFrame();
}

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_house_size_report() const
{
    using namespace Rcpp;

    try
    {
        size_t n_rows = _report_hh_size_.size();
        CharacterVector comment(n_rows);
        IntegerVector value(n_rows), counts(n_rows);
        NumericVector pct(n_rows);

        size_t k = 0;
        for (const auto &pair : _report_hh_size_)
        {
            comment[k] = pair.first;
            value[k] = pair.second.first;
            counts[k] = pair.second.second;
            pct[k] = 100. * static_cast<double>(pair.second.second) / n_individuals;
            ++k;
        }

        return DataFrame::create(_["comment"] = std::move(comment),
                                 _["value"] = std::move(value),
                                 _["counts"] = std::move(counts),
                                 _["pct"] = std::move(pct));
    }
    catch (const std::exception &e)
    {
        tools::error({"While exporting household size report:",
                      std::string(e.what())});
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Monte_Carlo_Build : traces exports  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export MCMC traces of parameter posterior samples. */
Rcpp::NumericMatrix MCMC::Impl::export_traces() const
{
    using namespace Rcpp;

    try
    {
        // UPDATE: I shall allow for export even if is_power_post
        if (is_power_post() || smpl_full_loglike.empty())
            throw std::length_error("No traces to export were found.");

        //----- Collect the column names -----//
        str_vec headers = {"loglikelihood"};
        headers += param_standard_infer_list + param_extra_infer_list;

        // augmented variable of id trace
        if (!smpl_augm_value.empty())
            for (const augmID &a : augm_list_infer + augm_list_func)
                headers += {a + "_id_" + setup_id_trace};

        // multiple chains
        if (setup_n_chains > 1)
            headers += {"CHAIN__"};

        //----- NumericMatrix with traces -----//
        NumericMatrix output(setup_n_chains * setup_c_sample_size,
                             headers.size());
        size_t col = 0;

        // loglikelihood
        output(_, col++) = NumericVector(wrap(smpl_full_loglike));

        // parameters
        for (const parID &p : param_standard_infer_list + param_extra_infer_list)
            output(_, col++) = NumericVector(wrap(smpl_param_value.at(p)));

        // augmented variables
        if (!smpl_augm_value.empty())
            for (const augmID &a : augm_list_infer + augm_list_func)
                output(_, col++) = NumericVector(wrap(smpl_augm_value.at(a)));

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
    catch (const std::exception &e)
    {
        tools::error({"While exporting MCMC traces:",
                      std::string(e.what())});
    }
    return NumericMatrix();
}

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix MCMC::Impl::export_house_loglike_mtx() const
{
    using namespace Rcpp;

    try
    {
        // UPDATE: I shall allow for export even if is_power_post
        if (is_power_post() || smpl_full_loglike.empty())
            tools::error("No traces to export were found.");

        // create a NumericMatrix with n_house columns
        ulong sample_size = setup_n_chains * setup_c_sample_size;
        NumericMatrix output(sample_size, n_households);
        CharacterVector headers(n_households);

        // assign each NumericVector to a column of the matrix
        for (size_t H = 0; H < n_households; ++H)
        {
            size_t H0 = house_vector[H].multiplicity_class_index();
            output(_, H) = NumericVector(wrap(house_vector[H0].get_smpl_loglike()));
            headers[H] = house_vector[H].house_id();
        }

        // associate column names
        colnames(output) = headers;

        // associate row names
        IntegerVector sweep(setup_c_sample_size * setup_n_chains);
        for (size_t c = 0; c < setup_n_chains; ++c)
            for (size_t k = 0; k < setup_c_sample_size; ++k)
                sweep[k + c * setup_c_sample_size] = setup_n_burnin + (k + 1) * setup_thin;
        rownames(output) = CharacterVector(wrap(sweep));

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error({"While exporting household loglikelihood matrix:",
                      std::string(e.what())});
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  Monte_Carlo_Build : infectivity profile export ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export a showcase of the selected infectivity profile. */
Rcpp::List MCMC::Impl::showcase_infec_profile(const Rcpp::List &opts, int n_bins) const
{
    using namespace Rcpp;

    try
    {
        // check the call makes sense
        if (_infec_profile.profile == nullptr)
            throw std::length_error("No infectivity profile to export.");

        // load options
        dbl_vec opts_ = {};
        for (int k = 0; k < opts.size(); ++k)
            opts_.push_back(opts.at(k));

        // produce a showcase of the distribution profile
        std::vector<dbl_vec> values = _infec_profile.profile->showcase(opts_, n_bins);

        // compute mean and variance
        const dbl_vec &x = values[0];
        const dbl_vec &prob = values[1];
        const size_t N = x.size();
        const double Dx = x[1] - x[0];

        double m = x[0] * prob[0] + x[N - 1] * prob[N - 1];
        double m2 = x[0] * x[0] * prob[0] + x[N - 1] * x[N - 1] * prob[N - 1];
        // integration wiht trapezoidal rule for continuos case
        if (!_infec_profile.profile->is_discrete())
        {
            m *= 0.5;
            m2 *= 0.5;
        }
        for (size_t i = 1; i < N - 1; ++i)
        {
            m += x[i] * prob[i];
            m2 += x[i] * x[i] * prob[i];
        }

        m *= Dx;
        m2 *= Dx;

        return List::create(
            _["name"] = _infec_profile.profile->show_name(),
            _["discr"] = _infec_profile.profile->is_discrete(),
            _["mean"] = m,
            _["var"] = m2 - m * m,
            _["gentime"] = wrap(values[0]),
            _["inst"] = wrap(values[1]),
            _["cumul"] = wrap(values[2]));
    }
    catch (const std::exception &e)
    {
        tools::error({"While exporting infectivity profile:",
                      std::string(e.what())});
    }
    return List();
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
Rcpp::DataFrame MCMC::export_param() const { return pImpl->export_param(); }
Rcpp::DataFrame MCMC::export_dataset() const { return pImpl->export_dataset(); }
Rcpp::DataFrame MCMC::export_data_report() const { return pImpl->export_data_report(); }
Rcpp::DataFrame MCMC::export_house_size_report() const { return pImpl->export_house_size_report(); }
Rcpp::NumericMatrix MCMC::export_traces() const { return pImpl->export_traces(); }
Rcpp::NumericMatrix MCMC::export_house_loglike_mtx() const { return pImpl->export_house_loglike_mtx(); }

Rcpp::List MCMC::showcase_infec_profile(const Rcpp::List &opts, int n_bins)
{
    init_infec_profile();
    return pImpl->showcase_infec_profile(opts, n_bins);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : imports from Rcpp ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::import_setup(const Rcpp::List &x) { pImpl->import_setup(x); }
void MCMC::import_infec_profile(const Rcpp::List &x) { pImpl->import_infec_profile(x); }
void MCMC::import_param(const Rcpp::DataFrame &x) { pImpl->import_param(x); }
void MCMC::import_augm(const Rcpp::DataFrame &x) { pImpl->import_augm(x); }
void MCMC::import_dataset(const Rcpp::DataFrame &x) { pImpl->import_dataset(x); }
void MCMC::import_multiplicity(const Rcpp::DataFrame &x) { pImpl->import_multiplicity(x); }

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
void MCMC::Impl::read_setup(const std::string &f)
{
    try
    {
        tools::read_csv_list(f, _setup_deflist);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading setup definitions:",
                      std::string(e.what())});
    }
    build_setup();
}

/* Read Setup file and produce a definition list. */
void MCMC::Impl::read_infec_profile(const std::string &f)
{
    try
    {
        tools::read_csv_list(f, _infec_profile.deflist);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading infectivity profile definitions:",
                      std::string(e.what())});
    }
    _infec_profile.build();
}

/* Read parameter definitions and run assignment. */
void MCMC::Impl::read_param(const std::string &f)
{
    try
    {
        tools::read_csv_table(f, _param_deftable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading parameter definitions:",
                      std::string(e.what())});
    }
    build_param();
}

/* Read Model Augmented Variable file. */
void MCMC::Impl::read_augm(const std::string &f)
{
    try
    {
        tools::read_csv_table(f, _augm_deftable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading augmented definitions:",
                      std::string(e.what())});
    }
    build_augm();
}

/* Read dataset file and then run assignment. */
void MCMC::Impl::read_dataset(const std::string &f)
{
    dict_str_vec datatable = {};
    try
    {
        tools::read_csv_table(f, datatable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading dataset:",
                      std::string(e.what())});
    }
    build_dataset(datatable);
}

/* Read multiplicity file and then run assignment. */
void MCMC::Impl::read_multiplicity(const std::string &f)
{
    dict_str_vec datatable = {};
    try
    {
        tools::read_csv_table(f, datatable);
    }
    catch (const std::exception &e)
    {
        tools::error({"While reading multiplicity:",
                      std::string(e.what())});
    }
    build_multiplicity(datatable);
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
    return " (" + tools::percentage(static_cast<double>(num) / den, 1) + ")";
}

/* Auxiliary function to call percentuage. */
inline const std::string _trail_percent_(int num, int den)
{
    return std::to_string(num) + _aux_percent_(num, den);
}

/* Print the report of the database and of the model group rates. */
void MCMC::Impl::show_summary() const
{
    std::ostringstream oss;

    oss << "n. of individuals: " << n_individuals
        << "\n    * n. of cases: "
        << _trail_percent_(get_report_indiv("infected"), n_individuals)
        << "\n    * n. of missing status : "
        << _trail_percent_(get_report_indiv("missing"), n_individuals);

    oss << "\nn. of households: " << n_households
        << "\n    * with no cases: "
        << _trail_percent_(get_report_house("no cases"), n_households)
        << "\n    * with at least 1 case: "
        << _trail_percent_(get_report_house("at least 1 case"), n_households)
        << "\n    * with at least 2 cases: "
        << _trail_percent_(get_report_house("at least 2 cases"), n_households);

    oss << "\nmax n. house cases: " << get_report_hh_size("max cases").first
        << _aux_percent_(get_report_hh_size("max cases").second, n_households)
        << "\nmin n. house members: " << get_report_hh_size("min size").first
        << _aux_percent_(get_report_hh_size("min size").second, n_households)
        << "\nmax n. house members: " << get_report_hh_size("max size").first
        << _aux_percent_(get_report_hh_size("max size").second, n_households);

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
        << "\nprogress" << C_DELIM << tools::truth(setup_progress)
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
            << C_DELIM << aux.show_value()
            << C_DELIM << tools::truth(aux.is_constant())
            << C_DELIM << tools::truth(aux.is_discrete())
            << C_DELIM << aux.show_distr_name()
            << C_DELIM << aux.show_distr_arg_names()
            << C_DELIM << aux.show_distr_arg_values()
            << C_DELIM << aux.show_prop_name()
            << C_DELIM << aux.show_prop_scale()
            << C_DELIM << aux.show_min()
            << C_DELIM << aux.show_max()
            << "\n";
    }

    tools::print_stdout(oss.str());
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
        const std::pair<double, double> stats = tools::compute_stats(smpl_param_value.at(p));

        oss << p
            << C_DELIM << stats.first
            << C_DELIM << std::sqrt(stats.second)
            << C_DELIM << tools::truth(this_param.is_constant())
            << C_DELIM << tools::truth(this_param.is_discrete())
            << C_DELIM << this_param.show_distr_name()
            << C_DELIM << this_param.show_distr_arg_names()
            << C_DELIM << this_param.show_distr_arg_values()
            << C_DELIM << this_param.show_min()
            << C_DELIM << this_param.show_max()
            << C_DELIM << this_param.show_prop_name()
            << C_DELIM << this_param.show_prop_scale()
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
                for (const augmID &a : augm_list_infer + augm_list_func)
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
                                        << std::fixed << smpl_full_loglike[i];

                    // print parameter values
                    for (const parID &p : param_standard_infer_list + param_extra_infer_list)
                        inferstate_ofstream << C_DELIM << std::setprecision(5)
                                            << std::fixed << smpl_param_value.at(p)[i];

                    // print selected augmented variables
                    if (is_augm_sample() && (!tools::is_na(setup_id_trace)))
                        for (const augmID &a : augm_list_infer + augm_list_func)
                            inferstate_ofstream << C_DELIM << std::setprecision(5)
                                                << std::fixed << smpl_augm_value.at(a)[i];

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
void MCMC::Impl::_aux_print_dataset_init_(const std::string &f)
{
    aux_data_ofstream.open(f.c_str());
    if (aux_data_ofstream.is_open())
    {
        // print default data headers
        aux_data_ofstream << "id" << C_DELIM << "house_id" << C_DELIM << "status";

        // print other data headers
        for (const dataID &d : data_headers)
            aux_data_ofstream << C_DELIM << d;

        // print augmented variables headers
        for (const augmID &a : augm_list_full)
            aux_data_ofstream << C_DELIM << a;

        aux_data_ofstream << std::endl;
    }
    else
        throw std::runtime_error(
            "Unable to open file: " + f);
}

/* Print all augmented variables of each individual */
void MCMC::Impl::_aux_print_dataset_update_(const Household &house)
{
    if (aux_data_ofstream.is_open())
    {
        for (indiv i = 0; i < house.size(); ++i)
        {
            aux_data_ofstream << house.id(i)
                              << C_DELIM << house.house_id()
                              << C_DELIM << house.status(i);
            for (const dataID &d : data_headers)
                aux_data_ofstream << C_DELIM << house.data(i, d);

            for (const augmID &a : augm_list_full)
                if (house.status(i) == INFECTED_STATUS)
                {
                    double x = house.augm(i, a);
                    aux_data_ofstream << C_DELIM
                                      << (tools::is_na(x) ? HOUSEPI_NA_DOUBLE : x);
                }
                else
                    aux_data_ofstream << C_DELIM << HOUSEPI_NA_DOUBLE;

            aux_data_ofstream << std::endl;
        }
    }
    else
        throw std::runtime_error(
            "Output file closed unexpectedly.");
}

/* Print all augmented variables of each individual */
void MCMC::Impl::_aux_print_dataset_finalize_()
{
    if (aux_data_ofstream.is_open())
        aux_data_ofstream.close();
    else
        throw std::runtime_error(
            "Output file closed unexpectedly.");
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

void MCMC::read_dataset(const std::string &f) { pImpl->read_dataset(f); }
void MCMC::read_setup(const std::string &f) { pImpl->read_setup(f); }
void MCMC::read_param(const std::string &f) { pImpl->read_param(f); }
void MCMC::read_infec_profile(const std::string &f) { pImpl->read_infec_profile(f); }
void MCMC::read_augm(const std::string &f) { pImpl->read_augm(f); }
void MCMC::read_multiplicity(const std::string &f) { pImpl->read_multiplicity(f); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : show, C++ print  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::show_setup() const { pImpl->show_setup(); }
void MCMC::show_all_params() const { pImpl->show_all_params(); }
void MCMC::show_results() const { pImpl->show_results(); }
void MCMC::show_summary() const { pImpl->show_summary(); }
void MCMC::print_infer_state(const std::string &f) const { pImpl->print_infer_state(f); }

#endif
