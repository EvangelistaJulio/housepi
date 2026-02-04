/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "tools/atomic.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/operators.h"
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
            std::vector<std::string> col_data(col.begin(), col.end());
            output[as<std::string>(colnames[i])] = std::move(col_data);
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
            std::vector<double> col_data(col.begin(), col.end());
            output[as<std::string>(colnames[i])] = std::move(col_data);
        }

        return output;
    }

    /// Export Function argument list to R.
    Rcpp::NumericVector rcpp_export_args(const Function *f)
    {
        using namespace Rcpp;

        size_t n_args = f->n_args();
        NumericVector vals(n_args);
        CharacterVector keys(n_args);
        for (size_t k = 0; k < n_args; ++k)
        {
            vals[k] = f->arg(k);
            keys[k] = std::string(f->arg_name(k));
        }
        vals.names() = std::move(keys);
        return vals;
    }

    /// Standaridized percentage formatting.
    double format_percentage(size_t value, size_t total)
    {
        return std::round(10000. * static_cast<double>(value) / static_cast<double>(total)) / 100.;
    }

    // Copies src into dest. If destroy_src is true, releases src's storage after copy.
    template <typename DestNumeric>
    void helper_vector_flush(DestNumeric dest, std::vector<double> &src)
    {
        std::copy(src.begin(), src.end(), dest.begin());
        std::vector<double>().swap(src);
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
            hsp::msg::error(std::string("While importing ") +   \
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
        pImpl_->data_wide->multiplicity_build(rdataframe_to_strmap(input),
                                              pImpl_->augm->pRegistry->size() > 0);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(std::string("While importing ") + "multiplicity" + ":\n" + e.what());
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
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    /// Export to R definitions of infectivity profile.
    Rcpp::List helper_export_infec_profile(const InfecProfile *infec_profile)
    {
        using namespace Rcpp;

        if (infec_profile->is_augmented())
            hsp::msg::error(
                "Impossible to export infectivity profile when if definition is not generic.");

        if (infec_profile->generic == nullptr)
            return List::create();

        return List::create(_["func"] = std::string(infec_profile->generic->name()),
                            // NOTE: arguments should not be polymorphic when called here
                            _["func_args"] = rcpp_export_args(infec_profile->generic.get()),
                            _["min"] = infec_profile->generic->lower(),
                            _["max"] = infec_profile->generic->upper(),
                            _["discr"] = infec_profile->generic->is_discrete());
    }
}

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

Rcpp::List MCMC::export_infec_profile() const { return helper_export_infec_profile(pImpl_->infec_profile.get()); }

/// Export to R all augmented variables at the current state.
Rcpp::DataFrame MCMC::export_augmented_data_wide() const
{
    using namespace Rcpp;
    const Augmented &augm = *pImpl_->augm;

    size_t n_cols = augm.list_export.size();
    size_t n_rows = pImpl_->data_wide->report_indiv->total;

    IntegerVector status_col(n_rows);
    std::vector<NumericVector> augm_cols(n_cols);
    for (size_t i = 0; i < n_cols; ++i)
        augm_cols[i] = NumericVector(n_rows);

    size_t row = 0;
    for (Household &house : *pImpl_->data_wide)
        for (Individual &member : house)
        {
            status_col[row] = member.status();
            size_t i = 0;
            for (size_t a : augm.list_export)
                augm_cols[i++][row] = hsp_to_rcpp(member.get_augm_variable(a).get());
            ++row;
        }

    // create the ouput List
    List output = List::create(_["status"] = status_col);
    size_t i = 0;
    for (size_t a : augm.list_export)
        output[augm.pRegistry->as_string(a)] = augm_cols[i++];
    return DataFrame(output);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : DataWide : exports  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame helper_export_data_wide(const DataWide &data_wide,
                                        const std::vector<std::string> &columns, bool is_status)
{
    using namespace Rcpp;

    try
    {
        if (data_wide.is_empty())
            return DataFrame();

        size_t n_rows = data_wide.report_indiv->total;
        CharacterVector df_id = CharacterVector(n_rows);
        CharacterVector df_house_id = CharacterVector(n_rows);
        IntegerVector df_status = IntegerVector(n_rows);

        // other non standard data columns
        std::vector<NumericVector> df_double(columns.size());
        for (size_t d = 0; d < columns.size(); ++d)
            df_double[d] = NumericVector(n_rows);

        size_t r = 0;
        for (const Household &house : data_wide)
            for (const Individual &member : house)
            {
                df_id[r] = member.id();
                df_house_id[r] = house.house_id();
                df_status[r] = member.status();

                const std::vector<double> &data_vec = member.get_data_wide_vec();
                for (size_t d = 0; d < columns.size(); ++d)
                    df_double[d][r] = hsp_to_rcpp(data_vec[d]);
                ++r;
            }

        DataFrame output = DataFrame::create(_["id"] = df_id,
                                             _["house_id"] = df_house_id);
        if (is_status)
            output.push_back(df_status, "status");
        for (size_t d = 0; d < columns.size(); ++d)
            output.push_back(df_double[d], columns[d]);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data (wide):\n" +
                        std::string(e.what()));
    }
}
Rcpp::DataFrame MCMC::export_data_wide() const
{
    const DataWide &data_wide = *pImpl_->data_wide;
    return helper_export_data_wide(data_wide, data_wide.pRegistry->get_all(), true);
}
Rcpp::DataFrame MCMC::data_wide_for_simulation() const
{
    const DataWide &data_wide = *pImpl_->data_wide;
    const Augmented &augm = *pImpl_->augm;

    std::vector<std::string> export_names(augm.list_export.size());
    size_t i = 0;
    for (size_t a : augm.list_export)
        export_names[i++] = augm.pRegistry->as_string(a);

    return helper_export_data_wide(data_wide, data_wide.pRegistry->get_all() - export_names, false);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : DataLong : exports  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::export_data_long() const
{
    using namespace Rcpp;
    const DataLong &data_long = *pImpl_->data_long;

    try
    {
        if (data_long.indiv_datalist.size() == 0)
            return DataFrame();

        // retrieve total number of rows
        size_t n_rows = 0;
        for (const auto &tbl : data_long.indiv_datalist)
            if (!tbl.empty())
                n_rows += tbl.begin()->size();

        CharacterVector datalist_index(n_rows);

        // other non standard data-long columns
        std::vector<NumericVector> datalist_dbl(data_long.pRegistry->size());
        for (size_t d = 0; d < datalist_dbl.size(); ++d)
            datalist_dbl[d] = NumericVector(n_rows);

        size_t r = 0;
        for (const Household &house : *pImpl_->data_wide)
            for (const Individual &member : house)
            {
                const std::vector<std::vector<double>> *ptr = member.get_data_long_ptr();
                if (ptr != nullptr && !ptr->empty())
                    for (size_t k = 0; k < ptr->begin()->size(); ++k)
                    {
                        datalist_index[r] = member.id();
                        for (size_t d = 0; d < datalist_dbl.size(); ++d)
                            datalist_dbl[d][r] = hsp_to_rcpp(ptr->at(d)[k]);
                        ++r;
                    }
            }

        // assign vectors to dataframe
        auto output = DataFrame::create(_["id"] = datalist_index);
        for (size_t d = 0; d < datalist_dbl.size(); ++d)
            output.push_back(datalist_dbl[d], data_long.pRegistry->as_string(d));

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data (long):\n" + std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : report exports  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::summarise_cases() const
{
    using namespace Rcpp;
    const ReportIndividual &report_indiv = *pImpl_->data_wide->report_indiv;
    const ReportHousehold &report_house = *pImpl_->data_wide->report_house;

    try
    {
        if (pImpl_->data_wide->is_empty())
            hsp::msg::warning("Called 'summarise_cases' on empty data.");

        CharacterVector type = wrap(std::vector<std::string>(4, "individual") +
                                    std::vector<std::string>(4, "household"));
        size_t n_rows = type.size();
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
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : export house size report  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::summarise_house_sizes() const
{
    using namespace Rcpp;
    const ReportIndividual &report_indiv = *pImpl_->data_wide->report_indiv;
    const ReportHouseSizes &report_hh_size = *pImpl_->data_wide->report_hh_size;

    try
    {
        CharacterVector comment({"max cases", "max miss", "min size", "max size"});
        size_t n_rows = comment.size();
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
                        std::string(e.what()));
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
        NumericMatrix output(pImpl_->tot_log_lik.count(), 1);
        helper_vector_flush(output(_, 0), pImpl_->tot_log_lik.smpl_values);
        colnames(output) = CharacterVector({"log_lik"});
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting total log-likelihood draws:\n" +
                        std::string(e.what()));
    }
}

/// Export MCMC draws of parameter values posterior samples.
Rcpp::NumericMatrix MCMC::flush_param_draws()
{
    using namespace Rcpp;
    Parameter &param = *pImpl_->param;

    if (param.list_sample.empty())
        return NumericMatrix();

    try
    {
        size_t n_cols = param.list_sample.size();
        CharacterVector headers = wrap(param.pRegistry->as_string(param.list_sample));
        NumericMatrix output(pImpl_->tot_log_lik.count(), n_cols);

        for (size_t c = 0; c < n_cols; ++c)
            helper_vector_flush(output(_, c), param.smpl_values[c]);

        colnames(output) = std::move(headers);
        param.smpl_values.clear();
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting parameters draws:\n" +
                        std::string(e.what()));
    }
}

/// Export MCMC draws of parameter log-priors posterior samples.
Rcpp::NumericMatrix MCMC::flush_param_lprior_draws()
{
    using namespace Rcpp;
    Parameter &param = *pImpl_->param;

    if (param.list_sample.empty())
        return NumericMatrix();

    try
    {
        size_t n_cols = param.list_sample.size();
        CharacterVector headers(n_cols);
        NumericMatrix output(pImpl_->tot_log_lik.count(), n_cols);

        size_t c = 0;
        for (size_t p : param.list_sample)
        {
            headers[c] = "lprior_" + param.pRegistry->as_string(p);
            helper_vector_flush(output(_, c), param.smpl_priors[c]);
            ++c;
        }
        colnames(output) = std::move(headers);
        param.smpl_priors.clear();
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting parameter log-priors draws:\n" +
                        std::string(e.what()));
    }
}

/// Export MCMC draws of augmented posterior samples.
Rcpp::NumericMatrix MCMC::flush_augm_draws()
{
    using namespace Rcpp;
    Augmented &augm = *pImpl_->augm;

    if (augm.smpl_values.empty())
        return NumericMatrix();

    try
    {
        size_t n_cols = augm.smpl_individuals.size() * augm.list_draws.size();
        CharacterVector headers(n_cols);
        NumericMatrix output(pImpl_->tot_log_lik.count(), n_cols);

        size_t col = 0;
        for (size_t i = 0; i < augm.smpl_values.size(); ++i)
        {
            const Individual *member_ptr = augm.smpl_individuals[i];
            auto &member_matrix = augm.smpl_values[i];

            size_t j = 0;
            for (size_t a : augm.list_draws)
            {
                headers[col] = member_ptr->id() + "_" + augm.pRegistry->as_string(a);
                helper_vector_flush(output(_, col++), member_matrix[j++]);
            }
        }

        colnames(output) = headers;
        augm.smpl_values.clear();
        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting augmented draws:\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : export household log-likelihood matrix  ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export all sampled values of the log-likelihood for all household.
Rcpp::NumericMatrix MCMC::flush_house_log_lik_draws()
{
    using namespace Rcpp;
    DataWide &data_wide = *pImpl_->data_wide;

    try
    {
        if (data_wide.smpl_values.empty())
            return NumericMatrix();

        // create a NumericMatrix with n_house columns
        size_t n_houses = data_wide.report_house->total;
        NumericMatrix output(pImpl_->tot_log_lik.count(), n_houses);
        CharacterVector headers(n_houses);

        // retrieve restricted index by general index of multiplicity head class
        size_t n_classes = data_wide.head_class_house_indexes.size();
        std::vector<std::vector<size_t>> house_classes(n_classes);

        std::unordered_map<size_t, size_t> restr_idx_map = {};
        restr_idx_map.reserve(n_classes);
        for (size_t R = 0; R < n_classes; ++R)
            restr_idx_map.emplace(data_wide.head_class_house_indexes[R], R);

        for (const Household &house : data_wide)
        {
            size_t R = restr_idx_map.at(house.multiplicity_class_index());
            size_t H = house.house_index();
            house_classes[R].push_back(H);
            headers[H] = house.house_id();
        }

        for (size_t R = 0; R < n_classes; ++R)
        {
            std::vector<double> &src = data_wide.smpl_values[R];

            for (size_t H : house_classes[R])
                std::copy(src.begin(), src.end(), output(_, H).begin());
            std::vector<double>().swap(src);
        }
        colnames(output) = std::move(headers);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting household log-likelihood draws:\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : house state snapshot /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{

    /// Export to R a snapshot of current state of house data.
    /// NOTE: there are not private parameters so far.
    Rcpp::NumericVector helper_param_rcpp_snapshot(const Parameter *param)
    {
        using namespace Rcpp;
        size_t n = param->pRegistry->size();

        if (n == 0)
            return NumericVector::create();

        NumericVector output(n);
        CharacterVector headers(n);
        for (size_t p = 0; p < n; ++p)
        {
            output[p] = param->at(p).get();
            headers[p] = std::string(param->at(p).name());
        }
        output.names() = CharacterVector(headers);

        return output;
    }

    /// Export to R a snapshot of current state of house data.
    Rcpp::DataFrame helper_data_wide_rcpp_snapshot(const DataWide *data_wide,
                                                   const Household &house)
    {
        using namespace Rcpp;
        const std::vector<std::string> &headers = data_wide->pRegistry->get_all();

        size_t n_rows = house.size();
        CharacterVector df_id = CharacterVector(n_rows);
        CharacterVector df_house_id = CharacterVector(n_rows);
        IntegerVector df_status = IntegerVector(n_rows);

        std::vector<NumericVector> df_double(headers.size());
        for (size_t d = 0; d < headers.size(); ++d)
            df_double[d] = NumericVector(n_rows);

        for (const Individual &member : house)
        {
            size_t i = member.hhindex();
            df_id[i] = member.id();
            df_house_id[i] = house.house_id();
            df_status[i] = member.status();

            const std::vector<double> &data_vec = member.get_data_wide_vec();
            for (size_t d = 0; d < headers.size(); ++d)
                df_double[d][i] = hsp_to_rcpp(data_vec[d]);
        }

        auto output = DataFrame::create(_["id"] = std::move(df_id),
                                        _["house_id"] = std::move(df_house_id),
                                        _["status"] = std::move(df_status));
        for (size_t d = 0; d < headers.size(); ++d)
            output.push_back(df_double[d], headers[d]);

        return output;
    }

    /// Export to R a snapshot of current state of house data-long.
    Rcpp::DataFrame helper_data_long_rcpp_snapshot(const DataLong *data_long,
                                                   const Household &house)
    {
        using namespace Rcpp;

        if (data_long->indiv_datalist.empty())
            return DataFrame::create();

        size_t n_rows = 0;
        for (const Individual &member : house)
        {
            const std::vector<std::vector<double>> *ptr = member.get_data_long_ptr();
            if (ptr != nullptr && !ptr->empty())
                n_rows += ptr->begin()->size();
        }

        CharacterVector datalist_id(n_rows);
        std::vector<NumericVector> df_double(data_long->pRegistry->size());
        for (size_t d = 0; d < df_double.size(); ++d)
            df_double[d] = NumericVector(n_rows);

        size_t r = 0;
        for (const Individual &member : house)
        {
            const std::vector<std::vector<double>> *ptr = member.get_data_long_ptr();
            if (ptr != nullptr && !ptr->empty())
                for (size_t k = 0; k < ptr->begin()->size(); ++k)
                {
                    for (size_t d = 0; d < df_double.size(); ++d)
                        df_double[d][r] = hsp_to_rcpp(ptr->at(d)[k]);
                    datalist_id[r++] = member.id();
                }
        }

        // assign vectors to dataframe
        auto output = DataFrame::create(_["id"] = std::move(datalist_id));
        for (size_t d = 0; d < df_double.size(); ++d)
            output.push_back(df_double[d], data_long->pRegistry->as_string(d));

        return output;
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_group_rates_rcpp_snapshot(const GroupRate *group_rate,
                                                     const Household &house)
    {
        using namespace Rcpp;

        const std::vector<std::string> &headers = group_rate->pRegistry->get_all();

        if (headers.empty())
            return DataFrame::create();

        size_t n_rows = house.size();
        CharacterVector df_id(n_rows);

        std::vector<NumericVector> df_double(headers.size());
        for (size_t c = 0; c < headers.size(); ++c)
            df_double[c] = NumericVector(n_rows);

        for (const Individual &member : house)
        {
            size_t i = member.hhindex();
            df_id[i] = member.id();

            for (size_t c = 0; c < headers.size(); ++c)
                df_double[c][i] = hsp_to_rcpp(member.get_group_rate_attribute(c).get());
        }

        auto output = DataFrame::create(_["id"] = df_id);
        for (size_t c = 0; c < headers.size(); ++c)
            output.push_back(df_double[c], headers[c]);

        return output;
    }

    /// Export to R a snapshot of current state of house infectivity profile.
    Rcpp::List helper_infec_profile_rcpp_snapshot(const InfecProfile *infec_profile,
                                                  const Household &house)
    {
        using namespace Rcpp;

        if (infec_profile->generic == nullptr && !infec_profile->is_augmented())
            return List::create();

        size_t n_rows = house.size();
        CharacterVector df_id(n_rows);

        CharacterVector infprof_name(n_rows, NA_STRING);
        LogicalVector infprof_discr(n_rows, NA_LOGICAL);
        List infprof_args(n_rows, NumericVector(0));
        NumericVector infprof_min(n_rows, NA_REAL);
        NumericVector infprof_max(n_rows, NA_REAL);

        for (size_t i = 0; i < house.size(); ++i)
            df_id[i] = house.member(i).id();

        for (size_t i : house.infected_hhindexes())
            if (infec_profile->is_augmented())
            {
                const auto &f = house.member(i).get_augm_variable(infec_profile->augm_key).get_function();

                infprof_name[i] = std::string(f->name());
                infprof_discr[i] = f->is_discrete();
                infprof_args[i] = rcpp_export_args(f.get());
                infprof_min[i] = f->lower();
                infprof_max[i] = f->upper();
            }
            else
            {
                List generic = helper_export_infec_profile(infec_profile);
                infprof_name[i] = as<std::string>(generic["func"]);
                infprof_discr[i] = as<bool>(generic["discr"]);
                infprof_args[i] = generic["func_args"];
                infprof_min[i] = as<double>(generic["min"]);
                infprof_max[i] = as<double>(generic["max"]);
            }

        return List::create(_["id"] = df_id,
                            _["discr"] = infprof_discr,
                            _["func"] = infprof_name,
                            _["func_args"] = infprof_args,
                            _["min"] = infprof_min,
                            _["max"] = infprof_max);
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_augm_rcpp_snapshot(const Augmented *augm, const Household &house)
    {
        using namespace Rcpp;

        size_t n_cols = augm->list_export.size();

        if (n_cols == 0)
            return DataFrame::create();

        size_t n_rows = house.size();
        CharacterVector df_id(n_rows);

        std::vector<NumericVector> df_double(n_cols);
        for (size_t c = 0; c < n_cols; ++c)
            df_double[c] = NumericVector(n_rows);

        for (size_t i = 0; i < house.size(); ++i)
        {
            const Individual &member_i = house.member(i);
            df_id[i] = member_i.id();
            size_t c = 0;
            for (size_t a : augm->list_export)
                df_double[c++][i] = hsp_to_rcpp(member_i.get_augm_variable(a).get());
        }

        auto output = DataFrame::create(_["id"] = df_id);
        size_t i = 0;
        for (size_t a : augm->list_export)
            output.push_back(df_double[i++], augm->pRegistry->as_string(a));

        return output;
    }
}

/// Export a list containing all relevant information for a debug.
Rcpp::List MCMC::export_house_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    try
    {
        if (pImpl_->data_wide->is_empty())
            return List::create();

        const Household &house = pImpl_->data_wide->at(house_id);

        return List::create(
            _["param"] = helper_param_rcpp_snapshot(pImpl_->param.get()),
            _["data_wide"] = helper_data_wide_rcpp_snapshot(pImpl_->data_wide.get(), house),
            _["data_long"] = helper_data_long_rcpp_snapshot(pImpl_->data_long.get(), house),
            _["group_rate"] = helper_group_rates_rcpp_snapshot(pImpl_->group_rate.get(), house),
            _["infec_profile"] = helper_infec_profile_rcpp_snapshot(pImpl_->infec_profile.get(), house),
            _["augm"] = helper_augm_rcpp_snapshot(pImpl_->augm.get(), house));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting snapshot of household '" + house_id + "':\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MCMC : simulate ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Perform simulation and return simulated data.
Rcpp::List MCMC::helper_simulate(size_t n, bool is_expec_measure)
{
    using namespace Rcpp;
    Parameter &param = *pImpl_->param;

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
        init_model_for_simulation();
        List output(n_rep);

        // run simulations
#pragma omp parallel for
        for (size_t r = 0; r < n_rep; ++r)
        {
            if (pImpl_->atomic_error->is_pending())
                continue;
            output[r] = helper_simulate_single_repetition();
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
Rcpp::DataFrame MCMC::helper_simulate_single_repetition()
{
    using namespace Rcpp;
    const Augmented &augm = *pImpl_->augm;
    DataWide &data_wide = *pImpl_->data_wide;
    DataFrame output;

    try
    {
        size_t n_cols = augm.list_export.size();
        size_t n_rows = data_wide.report_indiv->total;

        IntegerVector df_status(n_rows);

        std::vector<NumericVector> df_double(n_cols);
        for (size_t c = 0; c < n_cols; ++c)
            df_double[c] = NumericVector(n_rows);

        size_t r = 0;
        for (Household house : data_wide) // NOTE: deep copy of households
            try
            {
                data_wide.init_house_cache(house);
                simulate_single_house(house);
                for (const Individual &member : house)
                {
                    df_status[r] = member.status();
                    size_t c = 0;
                    for (size_t a : augm.list_export)
                        df_double[c++][r] = hsp_to_rcpp(member.get_augm_variable(a).get());
                    ++r;
                }
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating " + house.err_name() + ":\n" +
                    std::string(e.what()));
            }

        output = DataFrame::create(_["status"] = df_status);
        size_t c = 0;
        for (size_t a : augm.list_export)
            output.push_back(df_double[c++], augm.pRegistry->as_string(a));
    }
    catch (const std::exception &e)
    {
        pImpl_->atomic_error->set_if_first(e);
    }
    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  MCMC : log-evidence /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Estimate log_evidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
Rcpp::List MCMC::log_evidence(const std::string &scheme, unsigned long n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const auto pp = pImpl_->power_posteriors(scheme, n_rungs);
        const std::vector<double> &inv_tem = pp[0];
        const std::vector<double> &avg_vec = pp[1];
        const std::vector<double> &var_vec = pp[2];

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
                        std::string(e.what()));
    }
}
