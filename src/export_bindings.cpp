/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/atomic.h"
#include "tools/global.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/read.h"
#include "tools/types.h"

#include "function.h"
#include "distribution.h"
#include "data.h"
#include "data_reports.h"
#include "parameter.h"
#include "group_rate.h"
#include "infec_profile.h"
#include "augmented.h"
#include "mcmc_impl.h"
#include "mcmc.h"

#include <Rcpp.h>

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : imports from Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using std::string;
    using std::vector;

    /// Auxiliary converter to R NA real values.
    double hsp_to_rcpp(double x) { return hsp::na::is(x) ? NA_REAL : x; }

    /// Convert R dataframe map of string vectors.
    hsp::str_table rdataframe_to_strmap(SEXP input)
    {
        using namespace Rcpp;

        if (Rf_isNull(input))
            return {};

        DataFrame df(input);
        CharacterVector colnames = df.names();
        hsp::str_table output;

        for (R_xlen_t i = 0; i < df.size(); ++i)
        {
            CharacterVector col = df[i];
            vector<string> col_data(col.begin(), col.end());
            output[as<string>(colnames[i])] = std::move(col_data);
        }

        return output;
    }

    /// Convert R dataframe map of double vectors.
    hsp::dbl_table rdataframe_to_dblmap(SEXP input)
    {
        using namespace Rcpp;

        if (Rf_isNull(input))
            return {};

        DataFrame df(input);
        CharacterVector colnames = df.names();
        hsp::dbl_table output;

        for (R_xlen_t i = 0; i < df.size(); ++i)
        {
            NumericVector col = df[i];
            vector<double> col_data(col.begin(), col.end());
            output[as<string>(colnames[i])] = std::move(col_data);
        }

        return output;
    }

    /// Export Function argument list to R.
    Rcpp::NumericVector rcpp_export_args(const std::shared_ptr<Function> f)
    {
        using namespace Rcpp;

        size_t n_args = f->n_args();
        NumericVector vals(n_args);
        CharacterVector keys(n_args);
        for (size_t k = 0; k < n_args; ++k)
        {
            vals[k] = f->arg(k);
            keys[k] = f->arg_name(k);
        }
        vals.names() = std::move(keys);
        return vals;
    }

    /// Standaridized percentage formatting.
    double format_percentage(double value, double total) { return std::round(10000. * value / total) / 100.; }

    // Copies src into dest. If destroy_src is true, releases src's storage after copy.
    template <typename DestNumeric>
    void helper_vector_flush(DestNumeric dest, vector<double> &src)
    {
        std::copy(src.begin(), src.end(), dest.begin());
        vector<double>().swap(src);
    }
}

#define MAKE_RCPP_IMPORT_METHOD(module, section_name)           \
    void MCMC::import_##module(SEXP input)                      \
    {                                                           \
        try                                                     \
        {                                                       \
            pImpl_->module->build(rdataframe_to_strmap(input)); \
        }                                                       \
        catch (const std::exception &e)                         \
        {                                                       \
            hsp::msg::error(string("While importing ") +        \
                            section_name + ":\n" + e.what());   \
        }                                                       \
    }
MAKE_RCPP_IMPORT_METHOD(setup, "setup")
MAKE_RCPP_IMPORT_METHOD(data_wide, "data-wide")
MAKE_RCPP_IMPORT_METHOD(data_long, "data-long")
MAKE_RCPP_IMPORT_METHOD(param, "parameters")
MAKE_RCPP_IMPORT_METHOD(group_rate, "group rates")
MAKE_RCPP_IMPORT_METHOD(infec_profile, "infectivity profile")
MAKE_RCPP_IMPORT_METHOD(augm, "augmented")

// UPDATE: this will be removed once automatized in C++ reflection
void MCMC::import_multiplicity(SEXP input)
{
    try
    {
        pImpl_->data_wide->multiplicity_build(rdataframe_to_strmap(input), pImpl_->augm->is_any());
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(string("While importing ") + "multiplicity" + ":\n" + e.what());
    }
}

void MCMC::import_param_for_prediction(SEXP input)
{
    try
    {
        pImpl_->param->build_for_expec(rdataframe_to_dblmap(input));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While importanting parameter table:\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::CharacterVector MCMC::export_param_list_const() const { return Rcpp::wrap(pImpl_->param->list_const); }
Rcpp::CharacterVector MCMC::export_param_list_sample_standard() const { return Rcpp::wrap(pImpl_->param->list_sample_standard); }
Rcpp::CharacterVector MCMC::export_param_list_sample_cache() const { return Rcpp::wrap(pImpl_->param->list_sample_cached); }
Rcpp::CharacterVector MCMC::export_augm_list_draws() const { return Rcpp::wrap(pImpl_->augm->list_draws); }
Rcpp::CharacterVector MCMC::export_augm_list_full() const { return Rcpp::wrap(pImpl_->augm->list_full); }
Rcpp::CharacterVector MCMC::export_augm_list_export() const { return Rcpp::wrap(pImpl_->augm->list_export); }
Rcpp::CharacterVector MCMC::export_augm_list_discrete() const { return Rcpp::wrap(pImpl_->augm->user_list_discrete); }

/// Export current setup configurations to R List.
Rcpp::List MCMC::export_setup() const
{
    using namespace Rcpp;
    const Setup &setup = *pImpl_->setup;

    return List::create(_["steps"] = setup.n_steps,
                        _["burnin"] = setup.n_burnin,
                        _["thin"] = setup.n_thin,
                        _["adapt"] = setup.n_adapt,
                        _["chains"] = setup.n_chains,
                        _["threads"] = setup.n_threads,
                        _["id_draws"] = CharacterVector(wrap(setup.id_draws)),
                        _["seed"] = setup.seed,
                        _["progress"] = setup.progress);
}

/// Export to R definitions of infectivity profile.
/// FIXME: this is not working with new argument encoding.
Rcpp::List MCMC::export_infec_profile() const
{
    using namespace Rcpp;
    const InfecProfile &infec_profile = *pImpl_->infec_profile;

    if (infec_profile.is_augmented())
        hsp::msg::error(
            "Impossible to export infectivity profile when if definition is not generic.");

    if (infec_profile.generic_profile == nullptr ||
        infec_profile.generic_profile->is_empty())
        return List::create();

    return List::create(_["func"] = infec_profile.generic_profile->name(),
                        // NOTE: arguments should not be polymorphic when called here
                        _["func_args"] = rcpp_export_args(infec_profile.generic_profile),
                        _["min"] = infec_profile.generic_profile->lower(),
                        _["max"] = infec_profile.generic_profile->upper(),
                        _["discr"] = infec_profile.generic_profile->is_discrete());
}

/// Export to R all augmented variables at the current state.
Rcpp::DataFrame MCMC::export_augmented_data_wide() const
{
    using namespace Rcpp;
    const auto &augm = *pImpl_->augm;

    const size_t n_cols = pImpl_->data_wide->report_indiv->total;

    vector<int> status_col(n_cols);
    vector<vector<double>> augm_cols(augm.list_export.size());
    for (size_t i = 0; i < augm.list_export.size(); ++i)
        augm_cols[i] = vector<double>(n_cols);

    size_t row = 0;
    for (auto &house : *pImpl_->data_wide)
        for (const auto &member : house)
        {
            status_col[row] = member.status();
            for (size_t i = 0; i < augm.list_export.size(); ++i)
                augm_cols[i][row] = hsp_to_rcpp(member.aug(augm.list_export[i]));
            ++row;
        }

    // create the ouput List
    List output = List::create(_["status"] = IntegerVector(wrap(std::move(status_col))));
    for (size_t i = 0; i < augm.list_export.size(); ++i)
        output[augm.list_export[i]] = NumericVector(wrap(std::move(augm_cols[i])));
    return DataFrame(output);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : DataWide : exports  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame MCMC::helper_export_data_wide_(const vector<string> &columns, bool is_status) const
{
    using namespace Rcpp;
    const auto &data_wide = *pImpl_->data_wide;

    try
    {
        if (data_wide.is_empty())
            return DataFrame();

        const size_t n_rows = data_wide.report_indiv->total;
        CharacterVector df_id = CharacterVector(n_rows);
        CharacterVector df_house_id = CharacterVector(n_rows);
        IntegerVector df_status = IntegerVector(n_rows);

        // other non standard data columns
        vector<NumericVector> df_double;
        df_double.reserve(columns.size());
        for (size_t c = 0; c < columns.size(); ++c)
            df_double.emplace_back(n_rows);

        size_t r = 0;
        for (const auto &house : *pImpl_->data_wide)
            for (const auto &member : house)
            {
                df_id[r] = member.id();
                df_house_id[r] = house.house_id();
                df_status[r] = member.status();

                for (size_t c = 0; c < columns.size(); ++c)
                    df_double[c][r] = hsp_to_rcpp(member.data(columns[c]));
                ++r;
            }

        DataFrame output = DataFrame::create(_["id"] = df_id,
                                             _["house_id"] = df_house_id);
        if (is_status)
            output.push_back(df_status, "status");
        for (size_t c = 0; c < columns.size(); ++c)
            output.push_back(df_double[c], columns[c]);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data (wide):\n" +
                        string(e.what()));
    }
}
Rcpp::DataFrame MCMC::export_data_wide() const
{
    return helper_export_data_wide_(pImpl_->data_wide->headers, true);
}
Rcpp::DataFrame MCMC::data_wide_for_simulation() const
{
    return helper_export_data_wide_(pImpl_->data_wide->headers - pImpl_->augm->list_export, false);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : DataLong : exports  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::export_data_long() const
{
    using namespace Rcpp;
    const auto &data_long = *pImpl_->data_long;

    try
    {
        if (data_long.indiv_datalist.size() == 0)
            return DataFrame();

        // retrieve total number of rows
        size_t n_rows = 0;
        const auto &headers = data_long.headers;
        for (size_t i = 0; i < data_long.indiv_datalist.size(); ++i)
            n_rows += data_long.indiv_datalist[i].at(headers[0]).size();

        CharacterVector datalist_index(n_rows);

        // other non standard data-long columns
        vector<NumericVector> datalist_dbl = {};
        datalist_dbl.reserve(headers.size());
        for (size_t c = 0; c < headers.size(); ++c)
            datalist_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const auto &house : *pImpl_->data_wide)
            for (const auto &member : house)
                for (size_t k = 0; k < member.data_long(headers[0]).size(); ++k, ++r)
                {
                    datalist_index[r] = member.id();
                    for (size_t c = 0; c < headers.size(); ++c)
                        datalist_dbl[c][r] = hsp_to_rcpp(member.data_long(headers[c])[k]);
                }

        // assign vectors to dataframe
        auto output = DataFrame::create(_["id"] = datalist_index);
        for (size_t c = 0; c < headers.size(); ++c)
            output.push_back(datalist_dbl[c], headers[c]);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data (long):\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : report exports  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::summarise_cases() const
{
    using namespace Rcpp;
    const auto &report_indiv = *pImpl_->data_wide->report_indiv;
    const auto &report_house = *pImpl_->data_wide->report_house;

    try
    {
        if (pImpl_->data_wide->is_empty())
            hsp::msg::warning("Called 'summarise_cases' on empty data.");

        CharacterVector type = wrap(vector<string>(4, "individual") + vector<string>(4, "household"));
        const size_t n_rows = type.size();
        CharacterVector comment = {"total", "missing", "not infected", "infected",
                                   "total", "no cases", "at least 1 case", "at least 2 cases"};
        IntegerVector counts(n_rows);
        NumericVector pct(n_rows);

        counts[0] = report_indiv.total;
        pct[0] = 100.;

        counts[1] = report_indiv.missing;
        pct[1] = format_percentage(report_indiv.missing, report_indiv.total);

        counts[2] = report_indiv.susceptible;
        pct[2] = format_percentage(report_indiv.susceptible, report_indiv.total);

        counts[3] = report_indiv.infectious;
        pct[3] = format_percentage(report_indiv.infectious, report_indiv.total);

        counts[4] = report_house.total;
        pct[4] = 100.;
        counts[5] = report_house.no_cases;
        pct[5] = format_percentage(report_house.no_cases, report_house.total);

        counts[6] = report_house.at_least_1_case;
        pct[6] = format_percentage(report_house.at_least_1_case, report_house.total);

        counts[7] = report_house.at_least_2_cases;
        pct[7] = format_percentage(report_house.at_least_2_cases, report_house.total);

        return DataFrame::create(_["type"] = std::move(type),
                                 _["comment"] = std::move(comment),
                                 _["counts"] = std::move(counts),
                                 _["pct"] = std::move(pct));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data report:\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : export house size report  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::summarise_house_sizes() const
{
    using namespace Rcpp;
    const auto &report_indiv = *pImpl_->data_wide->report_indiv;
    const auto &report_hh_size = *pImpl_->data_wide->report_hh_size;

    try
    {
        CharacterVector comment({"max cases", "max miss", "min size", "max size"});
        const size_t n_rows = comment.size();
        IntegerVector value(n_rows), counts(n_rows);
        NumericVector pct(n_rows);

        value[0] = report_hh_size.max_cases.first;
        counts[0] = report_hh_size.max_cases.second;
        pct[0] = format_percentage(counts[0], report_indiv.total);

        value[1] = report_hh_size.max_miss.first;
        counts[1] = report_hh_size.max_miss.second;
        pct[1] = format_percentage(counts[1], report_indiv.total);

        value[2] = report_hh_size.min_size.first;
        counts[2] = report_hh_size.min_size.second;
        pct[2] = format_percentage(counts[2], report_indiv.total);

        value[3] = report_hh_size.max_size.first;
        counts[3] = report_hh_size.max_size.second;
        pct[3] = format_percentage(counts[3], report_indiv.total);

        return DataFrame::create(_["comment"] = std::move(comment),
                                 _["value"] = std::move(value),
                                 _["counts"] = std::move(counts),
                                 _["pct"] = std::move(pct));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting household size report:\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : export draws  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export MCMC draws of total log-likelihood posterior samples.
Rcpp::NumericMatrix MCMC::flush_total_log_lik_draws()
{
    using namespace Rcpp;

    try
    {
        NumericMatrix output(pImpl_->smpl_count, 1);
        helper_vector_flush(output(_, 0), pImpl_->smpl_total_log_lik);
        colnames(output) = CharacterVector({"log_lik"});
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting total log-likelihood draws:\n" +
                        string(e.what()));
    }
}

/// Export MCMC draws of parameter values posterior samples.
Rcpp::NumericMatrix MCMC::flush_param_draws()
{
    using namespace Rcpp;
    Parameter &param = *pImpl_->param;

    try
    {
        if (param.list_sample.empty())
            return NumericMatrix();

        const size_t n_cols = param.list_sample.size();
        CharacterVector headers(n_cols);
        NumericMatrix output(pImpl_->smpl_count, n_cols);

        size_t col = 0;
        for (const string &p : param.list_sample)
        {
            headers[col] = p;
            helper_vector_flush(output(_, col++), param.smpl_values.at(p));
        }
        colnames(output) = std::move(headers);
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting parameters draws:\n" +
                        string(e.what()));
    }
}

/// Export MCMC draws of parameter log-priors posterior samples.
Rcpp::NumericMatrix MCMC::flush_param_lprior_draws()
{
    using namespace Rcpp;
    auto &param = *pImpl_->param;

    try
    {
        if (param.list_sample.empty())
            return NumericMatrix();

        const size_t n_cols = param.list_sample.size();
        CharacterVector headers(n_cols);
        NumericMatrix output(pImpl_->smpl_count, n_cols);

        size_t col = 0;
        for (const string &p : param.list_sample)
        {
            headers[col] = "lprior_" + p;
            helper_vector_flush(output(_, col++), param.smpl_priors.at(p));
        }
        colnames(output) = std::move(headers);
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting parameter log-priors draws:\n" +
                        string(e.what()));
    }
}

/// Export MCMC draws of augmented posterior samples.
Rcpp::NumericMatrix MCMC::flush_augm_draws()
{
    using namespace Rcpp;
    auto &augm = *pImpl_->augm;

    try
    {
        if (augm.smpl_values.empty())
            return NumericMatrix();

        const size_t n_cols = augm.smpl_individuals.size() * augm.list_draws.size();
        CharacterVector headers(n_cols);
        NumericMatrix output(pImpl_->smpl_count, n_cols);

        size_t col = 0;
        for (size_t k = 0; k < augm.smpl_values.size(); ++k)
            for (const string &a : augm.list_draws)
            {
                headers[col] = augm.smpl_individuals[k]->id() + "_" + a;
                helper_vector_flush(output(_, col++), augm.smpl_values.at(k).at(a));
            }
        colnames(output) = std::move(headers);
        augm.smpl_values.clear();
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting augmented draws:\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : export household log-likelihood matrix  ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export all sampled values of the log-likelihood for all household.
Rcpp::NumericMatrix MCMC::flush_house_log_lik_draws()
{
    using namespace Rcpp;
    auto &data_wide = *pImpl_->data_wide;

    try
    {
        if (data_wide.smpl_values.empty())
            return NumericMatrix();

        // create a NumericMatrix with n_house columns
        auto n_houses = data_wide.report_house->total;
        NumericMatrix output(pImpl_->smpl_count, n_houses);
        CharacterVector headers(n_houses);

        // retrieve restricted index by general index of multiplicity head class
        const auto n_classes = data_wide.head_class_house_indexes.size();
        vector<vector<size_t>> house_classes(n_classes);

        std::unordered_map<size_t, size_t> restr_idx_map = {};
        restr_idx_map.reserve(n_classes);
        for (size_t R = 0; R < n_classes; ++R)
            restr_idx_map.emplace(data_wide.head_class_house_indexes[R], R);

        for (const auto &house : data_wide)
        {
            const auto R = restr_idx_map.at(house.multiplicity_class_index());
            const auto H = house.house_index();
            house_classes[R].push_back(H);
            headers[H] = house.house_id();
        }

        for (size_t R = 0; R < n_classes; ++R)
        {
            auto &src = data_wide.smpl_values[R];

            for (auto H : house_classes[R])
                std::copy(src.begin(), src.end(), output(_, H).begin());
            vector<double>().swap(src);
        }
        colnames(output) = std::move(headers);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting household log-likelihood draws:\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : house state snapshot /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export to R a snapshot of current state of house data.
Rcpp::NumericVector MCMC::helper_param_rcpp_snapshot_() const
{
    using namespace Rcpp;
    const auto &param = *pImpl_->param;
    const auto list_full = param.list_sample + param.list_const;

    if (list_full.empty())
        return NumericVector::create();

    NumericVector output(list_full.size());
    size_t c = 0;
    for (const string &p : list_full)
        output[c++] = param.get(p);
    output.names() = CharacterVector(wrap(list_full));

    return output;
}

/// Export to R a snapshot of current state of house data.
Rcpp::DataFrame MCMC::helper_data_wide_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const auto &data_wide = *pImpl_->data_wide;
    const auto &headers = data_wide.headers;

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);
    CharacterVector df_house_id = CharacterVector(n_rows);
    IntegerVector df_status = IntegerVector(n_rows);

    vector<NumericVector> df_double;
    df_double.reserve(headers.size());
    for (size_t c = 0; c < headers.size(); ++c)
        df_double.emplace_back(NumericVector(n_rows));

    for (const auto &member : house)
    {
        df_id[member] = member.id();
        df_house_id[member] = house.house_id();
        df_status[member] = member.status();

        for (size_t c = 0; c < headers.size(); ++c)
            df_double[c][member] = hsp_to_rcpp(member.data(headers[c]));
    }

    auto output = DataFrame::create(_["id"] = std::move(df_id),
                                    _["house_id"] = std::move(df_house_id),
                                    _["status"] = std::move(df_status));
    for (size_t c = 0; c < headers.size(); ++c)
        output.push_back(df_double[c], headers[c]);

    return output;
}

/// Export to R a snapshot of current state of house data-long.
Rcpp::DataFrame MCMC::helper_data_long_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const auto &data_long = *pImpl_->data_long;
    const auto &headers = data_long.headers;

    if (data_long.indiv_datalist.empty())
        return DataFrame::create();

    size_t n_rows = 0;

    for (const auto &member : house)
        n_rows += member.data_long(headers[0]).size();

    CharacterVector datalist_id(n_rows);

    vector<NumericVector> datalist_dbl = {};
    datalist_dbl.reserve(headers.size());
    for (size_t i = 0; i < headers.size(); ++i)
        datalist_dbl.emplace_back(n_rows);

    size_t r = 0;
    for (const auto &member : house)
        for (size_t i = 0; i < member.data_long(headers[0]).size(); ++i)
        {
            datalist_id[r] = member.id();
            for (size_t c = 0; c < headers.size(); ++c)
                datalist_dbl[c][r] = hsp_to_rcpp(member.data_long(headers[c])[i]);
            ++r;
        }

    // assign vectors to dataframe
    auto output = DataFrame::create(_["id"] = std::move(datalist_id));
    for (size_t c = 0; c < headers.size(); ++c)
        output.push_back(datalist_dbl[c], headers[c]);

    return output;
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::helper_group_rates_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;

    const auto &headers = pImpl_->group_rate->headers;

    if (headers.empty())
        return DataFrame::create();

    const auto n_rows = house.size();
    CharacterVector df_id(n_rows);

    vector<NumericVector> group_rate_dbl;
    group_rate_dbl.reserve(headers.size());
    for (size_t c = 0; c < headers.size(); ++c)
        group_rate_dbl.emplace_back(NumericVector(n_rows));

    for (const auto &member : house)
    {
        df_id[member] = member.id();

        for (size_t c = 0; c < headers.size(); ++c)
            group_rate_dbl[c][member] = hsp_to_rcpp(member.group_rate(headers[c]));
    }

    auto output = DataFrame::create(_["id"] = df_id);
    for (size_t c = 0; c < headers.size(); ++c)
        output.push_back(group_rate_dbl[c], headers[c]);

    return output;
}

/// Export to R a snapshot of current state of house infectivity profile.
Rcpp::List MCMC::helper_infec_profile_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const auto &infec_profile = *pImpl_->infec_profile;

    if (!infec_profile.is_augmented() && infec_profile.generic_profile->is_empty())
        return List::create();

    const size_t n_rows = house.size();
    CharacterVector df_id(n_rows);

    CharacterVector infprof_name(n_rows, NA_STRING);
    LogicalVector infprof_discr(n_rows, NA_LOGICAL);
    List infprof_args(n_rows, NumericVector(0));
    NumericVector infprof_min(n_rows, NA_REAL);
    NumericVector infprof_max(n_rows, NA_REAL);

    for (const auto &member : house)
        df_id[member] = member.id();

    for (const auto &member : house.infected())
        if (infec_profile.is_augmented())
        {
            const auto &f = member.get_augm_variable(hsp::infprof::INSTANT).get_function();

            infprof_name[member] = f->name();
            infprof_discr[member] = f->is_discrete();
            infprof_args[member] = rcpp_export_args(f);
            infprof_min[member] = f->lower();
            infprof_max[member] = f->upper();
        }
        else
        {
            List generic = export_infec_profile();
            infprof_name[member] = wrap(generic["func"]);
            infprof_discr[member] = generic["discr"];
            infprof_args[member] = generic["func_args"];
            infprof_min[member] = generic["min"];
            infprof_max[member] = generic["max"];
        }

    return List::create(_["id"] = df_id,
                        _["discr"] = infprof_discr,
                        _["func"] = infprof_name,
                        _["func_args"] = infprof_args,
                        _["min"] = infprof_min,
                        _["max"] = infprof_max);
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::helper_augm_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const auto &list_export = pImpl_->augm->list_export;

    if (list_export.empty())
        return DataFrame::create();

    const auto n_rows = house.size();
    CharacterVector df_id(n_rows);

    vector<NumericVector> augm_dbl;
    augm_dbl.reserve(list_export.size());
    for (size_t c = 0; c < list_export.size(); ++c)
        augm_dbl.emplace_back(NumericVector(n_rows));

    for (const auto &member : house)
        df_id[member] = member.id();

    for (const auto &member : house)
        for (size_t c = 0; c < list_export.size(); ++c)
            augm_dbl[c][member] = hsp_to_rcpp(member.aug(list_export[c]));

    auto output = DataFrame::create(_["id"] = df_id);
    for (size_t c = 0; c < list_export.size(); ++c)
        output.push_back(augm_dbl[c], list_export[c]);

    return output;
}

/// Export a list containing all relevant information for a debug.
Rcpp::List MCMC::export_house_snapshot(const string &house_id) const
{
    using namespace Rcpp;

    try
    {
        if (pImpl_->data_wide->is_empty())
            return List::create();

        const auto &house = (*pImpl_->data_wide)(house_id);

        return List::create(_["param"] = helper_param_rcpp_snapshot_(),
                            _["data_wide"] = helper_data_wide_rcpp_snapshot_(house),
                            _["data_long"] = helper_data_long_rcpp_snapshot_(house),
                            _["group_rate"] = helper_group_rates_rcpp_snapshot_(house),
                            _["infec_profile"] = helper_infec_profile_rcpp_snapshot_(house),
                            _["augm"] = helper_augm_rcpp_snapshot_(house));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting snapshot of household '" + house_id + "':\n" +
                        string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : simulate ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Perform simulation and return simulated data.
Rcpp::List MCMC::helper_simulate_(size_t n, bool is_expec_measure)
{
    using namespace Rcpp;
    auto &param = *pImpl_->param;

    pImpl_->atomic_error->restore();

    try
    {
        size_t n_rep = 1;
        // change execute mode
        // UPDATE: this could later enter in the loop
        if (is_expec_measure)
            param.load_sample_values(n);
        else
            n_rep = n;

        // initialization
        pImpl_->helper_init_simulation();
        List output(n_rep);

        // run simulations
#pragma omp parallel for
        for (size_t r = 0; r < n_rep; ++r)
        {
            if (pImpl_->atomic_error->is_pending())
                continue;
            output[r] = helper_simulate_single_repetition_();
        }
        pImpl_->atomic_error->throw_exception();

        return output;
    }
    catch (std::exception &e)
    {
        hsp::msg::error(e.what());
    }
}

/// Perform a unique repetition.
Rcpp::DataFrame MCMC::helper_simulate_single_repetition_()
{
    using namespace Rcpp;
    DataFrame output;
    try
    {
        const auto &augm_list_export = pImpl_->augm->list_export;
        const auto n_rows = pImpl_->data_wide->report_indiv->total;

        IntegerVector df_status(n_rows);

        vector<NumericVector> df_augm;
        df_augm.reserve(augm_list_export.size());
        for (size_t c = 0; c < augm_list_export.size(); ++c)
            df_augm.emplace_back(n_rows);

        size_t r = 0;
        for (auto house : *pImpl_->data_wide) // NOTE: deep copy of households
            try
            {
                pImpl_->simulate_single_house_in_place(house);
                for (const auto &member : house)
                {
                    df_status[r] = member.status();
                    for (size_t c = 0; c < augm_list_export.size(); ++c)
                        df_augm[c][r] = hsp_to_rcpp(member.aug(augm_list_export[c]));
                    ++r;
                }
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating " + house.err_name() + ":\n" + string(e.what()));
            }

        output = DataFrame::create(_["status"] = df_status);
        for (size_t c = 0; c < augm_list_export.size(); ++c)
            output.push_back(df_augm[c], augm_list_export[c]);
    }
    catch (const std::exception &e)
    {
        pImpl_->atomic_error->set_if_first(e);
    }
    return output;
}

/// Estimate log_evidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
Rcpp::List MCMC::log_evidence(const string &scheme, unsigned long n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const auto pp = pImpl_->power_posteriors(scheme, n_rungs);
        const auto &inv_tem = pp[0];
        const auto &avg_vec = pp[1];
        const auto &var_vec = pp[2];

        for (unsigned long k = 0; k < pImpl_->setup->n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(_["estimate"] = output,
                            _["scheme"] = scheme,
                            _["n_rungs"] = pImpl_->setup->c_sample_size,
                            _["inv_temper"] = wrap(inv_tem),
                            _["mean_log_lik"] = wrap(avg_vec),
                            _["var_log_lik"] = wrap(var_vec));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While estimating log_evidence with power posteriors:\n" +
                        string(e.what()));
    }
}
