/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_input_output.cpp
 *
 */

#include "tools/error.h"
#include "tools/read.h"
#include "mcmc_impl.h"

#ifdef _HOUSEPI_RCPP

#include "tools/params.h"
#include "tools/infprof.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : imports from Rcpp /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Convert R List into map<string, string>.
u_dict<std::string> rlist_to_map(SEXP input)
{
    using namespace Rcpp;

    if (Rf_isNull(input))
        return {};

    List list(input);
    CharacterVector rownames = list.names();
    u_dict<std::string> output = {};

    for (R_xlen_t i = 0; i < list.size(); ++i)
    {
        CharacterVector element = list[i];
        std::string value = as<std::string>(wrap(element[0]));
        output[as<std::string>(rownames[i])] = value;
    }

    return output;
}

/// Convert R dataframe into unordered_map<string, vector<string>>.
u_dict<str_vec> rdataframe_to_map(SEXP input)
{
    using namespace Rcpp;

    if (Rf_isNull(input))
        return {};

    DataFrame df(input);
    CharacterVector colnames = df.names();
    u_dict<str_vec> output = {};

    for (R_xlen_t i = 0; i < df.size(); ++i)
    {
        CharacterVector col_str = as<CharacterVector>(df[i]);
        str_vec col_data(col_str.begin(), col_str.end());
        output[as<std::string>(colnames[i])] = std::move(col_data);
    }

    return output;
}

#define MAKE_RCPP_IMPORT_METHOD(ImplSection, section_name, r_to_map) \
    void MCMC::Impl::ImplSection::rcpp_import(SEXP input)            \
    {                                                                \
        try                                                          \
        {                                                            \
            build(r_to_map(input));                                  \
        }                                                            \
        catch (const std::exception &e)                              \
        {                                                            \
            hsp::msg::error(                                         \
                "While reading " + section_name_ + ":\n" +           \
                std::string(e.what()));                              \
        }                                                            \
    }                                                                \
    void MCMC::import_##section_name(SEXP input)                     \
    {                                                                \
        pImpl_->section_name.rcpp_import(input);                     \
    }

MAKE_RCPP_IMPORT_METHOD(Setup, setup, rlist_to_map)
MAKE_RCPP_IMPORT_METHOD(InfecProfile, infec_profile, rlist_to_map)
MAKE_RCPP_IMPORT_METHOD(Parameter, param, rdataframe_to_map)
MAKE_RCPP_IMPORT_METHOD(DataWide, data_wide, rdataframe_to_map)
MAKE_RCPP_IMPORT_METHOD(Multiplicity, multiplicity, rdataframe_to_map)
MAKE_RCPP_IMPORT_METHOD(DataLong, data_long, rdataframe_to_map)
MAKE_RCPP_IMPORT_METHOD(GroupRates, group_rates, rdataframe_to_map)
MAKE_RCPP_IMPORT_METHOD(Augmented, augm, rdataframe_to_map)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export Function argument list to R.
Rcpp::NumericVector rcpp_export_args(const Function &f)
{
    using namespace Rcpp;

    size_t n_args = f.n_args();
    NumericVector vals(n_args);
    CharacterVector keys(n_args);
    for (size_t k = 0; k < n_args; ++k)
    {
        vals[k] = f.arg(k);
        keys[k] = f.arg_name(k);
    }
    vals.names() = std::move(keys);
    return vals;
}
/// Export arguments of a Variable which is a function to R.
Rcpp::NumericVector rcpp_export_func_args(const Variable &v) { return rcpp_export_args(*v.get_function()); }
/// Export arguments of a Variable model to R.
Rcpp::NumericVector rcpp_export_distr_args(const Variable &v) { return rcpp_export_args(*v.get_distribution()); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : exports to Rcpp ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::CharacterVector MCMC::export_param_list_const() const { return Rcpp::wrap(pImpl_->param.list_const); }
Rcpp::CharacterVector MCMC::export_param_list_sample_standard() const { return Rcpp::wrap(pImpl_->param.list_sample_standard); }
Rcpp::CharacterVector MCMC::export_param_list_sample_cache() const { return Rcpp::wrap(pImpl_->param.list_sample_cached); }
Rcpp::CharacterVector MCMC::export_augm_list_draws() const { return Rcpp::wrap(pImpl_->augm.list_draws); }
Rcpp::CharacterVector MCMC::export_augm_list_full() const { return Rcpp::wrap(pImpl_->augm.list_full); }
Rcpp::CharacterVector MCMC::export_augm_list_export() const { return Rcpp::wrap(pImpl_->augm.list_export); }
Rcpp::CharacterVector MCMC::export_augm_list_discrete() const { return Rcpp::wrap(pImpl_->augm.user_list_discrete); }

/// Export current setup configurations to R List.
Rcpp::List MCMC::Impl::Setup::rcpp_export() const
{
    using namespace Rcpp;

    return List::create(_["steps"] = n_steps,
                        _["burnin"] = n_burnin,
                        _["thin"] = n_thin,
                        _["adapt"] = n_adapt,
                        _["chains"] = n_chains,
                        _["threads"] = n_threads,
                        _["id_draws"] = CharacterVector(wrap(id_draws)),
                        _["seed"] = seed,
                        _["progress"] = progress);
}
Rcpp::List MCMC::export_setup() const { return pImpl_->setup.rcpp_export(); }

/// Export to R definitions of infectivity profile.
/// FIXME: this is not working with new argument encoding.
Rcpp::List MCMC::Impl::InfecProfile::rcpp_export() const
{
    using namespace Rcpp;

    if (is_augmented())
        hsp::msg::error(
            "Impossible to export infectivity profile when if definition is not generic.");

    if (profile_->is_empty())
        return List::create();

    return List::create(_["name"] = profile_->name(),
                        _["args"] = rcpp_export_args(*profile_),
                        _["min"] = profile_->lower(),
                        _["max"] = profile_->upper(),
                        _["discr"] = profile_->is_discrete());
}
Rcpp::List MCMC::export_infec_profile() const { return pImpl_->infec_profile.rcpp_export(); }

/// Export to R all augmented variables at the current state.
Rcpp::DataFrame MCMC::Impl::Augmented::rcpp_export_data_wide() const
{
    using namespace Rcpp;

    const size_t n_cols = pImpl_.data_wide.report_indiv.total;

    int_vec status_col(n_cols);
    std::vector<dbl_vec> augm_cols(list_export.size());
    for (size_t i = 0; i < list_export.size(); ++i)
        augm_cols[i] = dbl_vec(n_cols);

    size_t row = 0;
    for (Household &house : pImpl_.data_wide)
        for (const Individual &member : house)
        {
            status_col[row] = member.status();
            for (size_t i = 0; i < list_export.size(); ++i)
                augm_cols[i][row] = hsp_to_rcpp(member.augm(list_export[i]));
            ++row;
        }

    // create the ouput List
    List output = List::create(_["status"] = IntegerVector(wrap(std::move(status_col))));
    for (size_t i = 0; i < list_export.size(); ++i)
        output[list_export[i]] = NumericVector(wrap(std::move(augm_cols[i])));
    return DataFrame(output);
}
Rcpp::DataFrame MCMC::export_augmented_data_wide() const { return pImpl_->augm.rcpp_export_data_wide(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : DataWide : exports  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_export(const str_vec &columns, bool is_status) const
{
    using namespace Rcpp;

    try
    {
        if (is_empty())
            return DataFrame();

        const size_t n_rows = report_indiv.total;
        CharacterVector df_id = CharacterVector(n_rows);
        CharacterVector df_house_id = CharacterVector(n_rows);
        IntegerVector df_status = IntegerVector(n_rows);

        // other non standard data columns
        std::vector<NumericVector> df_double;
        df_double.reserve(columns.size());
        for (size_t c = 0; c < columns.size(); ++c)
            df_double.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : house_vector_)
            for (const Individual &member : house)
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
                        std::string(e.what()));
    }
}
Rcpp::DataFrame MCMC::export_data_wide() const
{
    return pImpl_->data_wide.rcpp_export(pImpl_->data_wide.headers(), true);
}
Rcpp::DataFrame MCMC::data_wide_for_simulation() const
{
    return pImpl_->data_wide.rcpp_export(pImpl_->data_wide.headers() - pImpl_->augm.list_export, false);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : DataLong : exports  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::Impl::DataLong::rcpp_export() const
{
    using namespace Rcpp;

    try
    {
        if (indiv_datalist_.size() == 0)
            return DataFrame();

        // retrieve total number of rows
        size_t n_rows = 0;
        for (size_t i = 0; i < pImpl_.data_wide.report_indiv.total; ++i)
            n_rows += indiv_datalist_[i].at(headers_[0]).size();

        CharacterVector datalist_index = CharacterVector(n_rows);

        // other non standard data-long columns
        std::vector<NumericVector> datalist_dbl = {};
        datalist_dbl.reserve(headers_.size());
        for (size_t c = 0; c < headers_.size(); ++c)
            datalist_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : pImpl_.data_wide)
            for (const Individual &member : house)
                for (size_t k = 0; k < member.data_long(headers_[0]).size(); ++k, ++r)
                {
                    datalist_index[r] = member.id();
                    for (size_t c = 0; c < headers_.size(); ++c)
                        datalist_dbl[c][r] = hsp_to_rcpp(member.data_long(headers_[c])[k]);
                }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(_["id"] = datalist_index);
        for (size_t c = 0; c < headers_.size(); ++c)
            output.push_back(datalist_dbl[c], headers_[c]);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data (long):\n" +
                        std::string(e.what()));
    }
}
Rcpp::DataFrame MCMC::export_data_long() const { return pImpl_->data_long.rcpp_export(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : report exports  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double format_percentage(double value, double total) { return std::round(10000. * value / total) / 100.; }

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_summarise_cases() const
{
    using namespace Rcpp;

    try
    {
        if (is_empty())
            hsp::msg::warning("Called 'summarise_cases' on empty data.");

        CharacterVector type = wrap(str_vec(4, "individual") + str_vec(4, "household"));
        const size_t n_rows = type.size();
        CharacterVector comment = {"total", "missing", "healthy", "infected",
                                   "total", "no cases", "at least 1 case", "at least 2 cases"};
        IntegerVector counts(n_rows);
        NumericVector pct(n_rows);

        counts[0] = report_indiv.total;
        pct[0] = 100.;

        counts[1] = report_indiv.missing;
        pct[1] = format_percentage(report_indiv.missing, report_indiv.total);

        counts[2] = report_indiv.healthy;
        pct[2] = format_percentage(report_indiv.healthy, report_indiv.total);

        counts[3] = report_indiv.infected;
        pct[3] = format_percentage(report_indiv.infected, report_indiv.total);

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
Rcpp::DataFrame MCMC::summarise_cases() const { return pImpl_->data_wide.rcpp_summarise_cases(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : export house size report  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_summarise_house_sizes() const
{
    using namespace Rcpp;

    try
    {
        CharacterVector comment({"max cases", "max miss", "min size", "max size"});
        const size_t n_rows = comment.size();
        IntegerVector value(n_rows), counts(n_rows);
        NumericVector pct(n_rows);

        value[0] = report_hh_size.max_cases.first;
        counts[0] = report_hh_size.max_cases.second;
        pct[0] = format_percentage(report_hh_size.max_cases.second, report_indiv.total);

        value[1] = report_hh_size.max_miss.first;
        counts[1] = report_hh_size.max_miss.second;
        pct[1] = format_percentage(report_hh_size.max_miss.second, report_indiv.total);

        value[2] = report_hh_size.min_size.first;
        counts[2] = report_hh_size.min_size.second;
        pct[2] = format_percentage(report_hh_size.min_size.second, report_indiv.total);

        value[3] = report_hh_size.max_size.first;
        counts[3] = report_hh_size.max_size.second;
        pct[3] = format_percentage(report_hh_size.max_size.second, report_indiv.total);

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
Rcpp::DataFrame MCMC::summarise_house_sizes() const { return pImpl_->data_wide.rcpp_summarise_house_sizes(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : export draws  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Chain tag of each trace entry.
Rcpp::IntegerVector MCMC::Impl::Setup::draws_chain_tags() const
{
    using namespace Rcpp;

    IntegerVector output(c_sample_size * n_chains);
    for (size_t c = 0; c < n_chains; ++c)
        std::fill(output.begin() + c * c_sample_size,
                  output.begin() + (c + 1) * c_sample_size, c + 1);
    return output;
}

/// Iteration name of each trace entry.
Rcpp::CharacterVector MCMC::Impl::Setup::draws_iter_names() const
{
    using namespace Rcpp;

    std::vector<unsigned long> iteration(c_sample_size * n_chains);
    for (size_t c = 0; c < n_chains; ++c)
        for (size_t k = 0; k < c_sample_size; k++)
            iteration[k + c * c_sample_size] = n_burnin + (k + 1) * n_thin;
    return CharacterVector(wrap(std::move(iteration)));
}

/// Export MCMC draws of parameter posterior samples.
Rcpp::NumericMatrix MCMC::Impl::rcpp_export_draws() const
{
    using namespace Rcpp;

    try
    {
        // UPDATE: I shall allow for export even if _is_power_post
        if (smpl_total_log_lik_.empty())
            throw std::length_error("No draws to export were found.");

        // define matrix headers
        str_vec headers = {"log_lik"};
        headers += param.list_sample;
        // prior of sampled parameters
        for (const std::string &p : param.list_sample)
            headers.push_back("lprior_" + p);
        // augmented variable of id trace
        if (!augm.smpl_values.empty())
            for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                for (const std::string &a : augm.list_draws)
                    headers.push_back(augm.smpl_individuals[k]->id() + "_" + a);
        // multiple chains
        headers.push_back(hsp::read::CHAIN_COLUMN);

        // matrix with draws
        NumericMatrix output(setup.n_chains * setup.c_sample_size, headers.size());

        size_t col = 0;
        // log-likelihood
        output(_, col++) = NumericVector(wrap(smpl_total_log_lik_));
        // parameter values
        for (const std::string &p : param.list_sample)
            output(_, col++) = NumericVector(wrap(param.smpl_values.at(p)));
        // parameter values
        for (const std::string &p : param.list_sample)
            output(_, col++) = NumericVector(wrap(param.smpl_priors.at(p)));

        // augmented variables
        if (!augm.smpl_values.empty())
            for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                for (const std::string &a : augm.list_draws)
                    output(_, col++) = NumericVector(wrap(augm.smpl_values.at(k).at(a)));

        // asssociate chain tag
        output(_, col++) = setup.draws_chain_tags();
        // associate column names
        colnames(output) = CharacterVector(wrap(headers));
        // associate row names
        rownames(output) = setup.draws_iter_names();

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting MCMC draws:\n" +
                        std::string(e.what()));
    }
}
Rcpp::NumericMatrix MCMC::export_draws() const { return pImpl_->rcpp_export_draws(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  MCMC : Impl : export household log-likelihood matrix  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export all sampled values of the log-likelihood for all household.
Rcpp::NumericMatrix MCMC::Impl::DataWide::rcpp_export_draws() const
{
    using namespace Rcpp;

    try
    {
        if (smpl_values_.empty())
            throw std::length_error("No draws to export were found.");

        // create a NumericMatrix with n_house columns
        const unsigned long sample_size = pImpl_.setup.n_chains * pImpl_.setup.c_sample_size;
        NumericMatrix output(sample_size, report_house.total + 1);
        CharacterVector columns(report_house.total + 1);
        columns[report_house.total] = hsp::read::CHAIN_COLUMN;

        // retrieve restricted index by general index of multiplicity head class
        std::unordered_map<size_t, size_t> restr_idx_map = {};
        restr_idx_map.reserve(pImpl_.multiplicity.indexes.size());
        for (size_t R = 0; R < pImpl_.multiplicity.indexes.size(); ++R)
            restr_idx_map.emplace(pImpl_.multiplicity.indexes[R], R);

        // assign each NumericVector to a column of the matrix
        for (size_t H = 0; H < report_house.total; ++H)
        {
            size_t R = restr_idx_map.at(house_vector_[H].multiplicity_class_index());
            output(_, H) = NumericVector(wrap(smpl_values_[R]));
            columns[H] = house_vector_[H].house_id();
        }

        // asssociate chain tag
        output(_, report_house.total) = pImpl_.setup.draws_chain_tags();
        // associate column names
        colnames(output) = columns;
        // associate row names
        rownames(output) = pImpl_.setup.draws_iter_names();

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting household log-likelihood draws:\n" +
                        std::string(e.what()));
    }
}
Rcpp::NumericMatrix MCMC::export_house_log_lik_draws() const { return pImpl_->data_wide.rcpp_export_draws(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : house state snapshot //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export to R a snapshot of current state of house data.
Rcpp::NumericVector MCMC::Impl::Parameter::rcpp_snapshot() const
{
    using namespace Rcpp;

    if (dictionary_.empty())
        return NumericVector::create();

    NumericVector output = NumericVector(dictionary_.size());
    CharacterVector names(dictionary_.size());

    size_t c = 0;
    for (const auto &this_param : dictionary_)
    {
        names[c] = this_param.first;
        output[c++] = this_param.second.get();
    }
    output.names() = names;

    return output;
}
Rcpp::List MCMC::export_house_snapshot(const std::string &house_id) const { return pImpl_->rcpp_export_house_snapshot(house_id); }

/// Export to R a snapshot of current state of house data.
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_snapshot(const Household &house) const
{
    using namespace Rcpp;

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);
    CharacterVector df_house_id = CharacterVector(n_rows);
    IntegerVector df_status = IntegerVector(n_rows);

    std::vector<NumericVector> df_double;
    df_double.reserve(headers_.size());
    for (size_t c = 0; c < headers_.size(); ++c)
        df_double.emplace_back(NumericVector(n_rows));

    for (const Individual &member : house)
    {
        df_id[member] = member.id();
        df_house_id[member] = house.house_id();
        df_status[member] = member.status();

        for (size_t c = 0; c < headers_.size(); ++c)
            df_double[c][member] = hsp_to_rcpp(member.data(headers_[c]));
    }

    DataFrame output = DataFrame::create(_["id"] = std::move(df_id),
                                         _["house_id"] = std::move(df_house_id),
                                         _["status"] = std::move(df_status));
    for (size_t c = 0; c < headers_.size(); ++c)
        output.push_back(df_double[c], headers_[c]);

    return output;
}

/// Export to R a snapshot of current state of house data-long.
Rcpp::DataFrame MCMC::Impl::DataLong::rcpp_snapshot(const Household &house) const
{
    using namespace Rcpp;

    if (is_empty())
        return DataFrame::create();

    size_t n_rows = 0;

    for (const Individual &member : house)
        n_rows += member.data_long(headers_[0]).size();

    CharacterVector datalist_id = CharacterVector(n_rows);

    std::vector<NumericVector> datalist_dbl = {};
    datalist_dbl.reserve(headers_.size());
    for (size_t i = 0; i < headers_.size(); ++i)
        datalist_dbl.emplace_back(n_rows);

    size_t r = 0;
    for (const Individual &member : house)
        for (size_t i = 0; i < member.data_long(headers_[0]).size(); ++i)
        {
            datalist_id[r] = member.id();
            for (size_t c = 0; c < headers_.size(); ++c)
                datalist_dbl[c][r] = hsp_to_rcpp(member.data_long(headers_[c])[i]);
            ++r;
        }

    // assign vectors to dataframe
    DataFrame output = DataFrame::create(_["id"] = std::move(datalist_id));
    for (size_t c = 0; c < headers_.size(); ++c)
        output.push_back(datalist_dbl[c], headers_[c]);

    return output;
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::Impl::GroupRates::rcpp_snapshot(const Household &house) const
{
    using namespace Rcpp;

    if (headers_.empty())
        return DataFrame::create();

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);

    std::vector<NumericVector> group_rate_dbl;
    group_rate_dbl.reserve(headers_.size());
    for (size_t c = 0; c < headers_.size(); ++c)
        group_rate_dbl.emplace_back(NumericVector(n_rows));

    for (const Individual &member : house)
    {
        df_id[member] = member.id();

        for (size_t c = 0; c < headers_.size(); ++c)
            group_rate_dbl[c][member] = hsp_to_rcpp(member.group_rate(headers_[c]));
    }

    DataFrame output = DataFrame::create(_["id"] = df_id);
    for (size_t c = 0; c < headers_.size(); ++c)
        output.push_back(group_rate_dbl[c], headers_[c]);

    return output;
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::Impl::Augmented::rcpp_snapshot(const Household &house) const
{
    using namespace Rcpp;

    if (list_export.empty())
        return DataFrame::create();

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);

    std::vector<NumericVector> augm_dbl;
    augm_dbl.reserve(list_export.size());
    for (size_t c = 0; c < list_export.size(); ++c)
        augm_dbl.emplace_back(NumericVector(n_rows, NA_REAL));

    for (const Individual &member : house)
        df_id[member] = member.id();

    for (const Individual &member : house.infected())
        for (size_t c = 0; c < list_export.size(); ++c)
            augm_dbl[c][member] = hsp_to_rcpp(member.augm(list_export[c]));

    DataFrame output = DataFrame::create(_["id"] = df_id);
    for (size_t c = 0; c < list_export.size(); ++c)
        output.push_back(augm_dbl[c], list_export[c]);

    return output;
}

/// Export to R a snapshot of current state of house infectivity profile.
Rcpp::List MCMC::Impl::InfecProfile::rcpp_snapshot(const Household &house) const
{
    using namespace Rcpp;

    if (!is_augmented() && profile_->is_empty())
        return List::create();

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);

    CharacterVector infprof_name(n_rows, NA_STRING);
    LogicalVector infprof_discr(n_rows, NA_LOGICAL);
    List infprof_args(n_rows, NumericVector(0));
    NumericVector infprof_min(n_rows, NA_REAL);
    NumericVector infprof_max(n_rows, NA_REAL);

    for (const Individual &member : house)
        df_id[member] = member.id();

    for (const Individual &member : house.infected())
        if (is_augmented())
        {
            const Variable &var = member.get_augm_variable(hsp::infprof::INSTANTANEOUS);
            infprof_name[member] = var.func_name();
            infprof_discr[member] = var.is_discrete();
            infprof_args[member] = rcpp_export_func_args(var);
            infprof_min[member] = var.get_distribution()->lower();
            infprof_max[member] = var.get_distribution()->upper();
        }
        else
        {
            infprof_name[member] = profile_->name();
            infprof_discr[member] = profile_->is_discrete();
            infprof_args[member] = rcpp_export_args(*profile_);
            infprof_min[member] = profile_->lower();
            infprof_max[member] = profile_->upper();
        }

    return List::create(_["id"] = df_id,
                        _["discr"] = infprof_discr,
                        _["name"] = infprof_name,
                        _["args"] = infprof_args,
                        _["min"] = infprof_min,
                        _["max"] = infprof_max);
}

/// Export a list containing all relevant information for a debug.
Rcpp::List MCMC::Impl::rcpp_export_house_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    try
    {
        if (data_wide.is_empty())
            return List::create();

        const Household &house = data_wide(house_id);

        return List::create(_["param"] = param.rcpp_snapshot(),
                            _["data_wide"] = data_wide.rcpp_snapshot(house),
                            _["data_long"] = data_long.rcpp_snapshot(house),
                            _["group_rates"] = group_rates.rcpp_snapshot(house),
                            _["infec_profile"] = infec_profile.rcpp_snapshot(house),
                            _["augm"] = augm.rcpp_snapshot(house));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting snaposhot of household '" + house_id + "':\n" +
                        std::string(e.what()));
    }
}

#else

#include <fstream>
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : read  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_READ_METHOD(ImplSection, section_name, csv_object)        \
    void MCMC::Impl::ImplSection::read(const std::string &f)           \
    {                                                                  \
        try                                                            \
        {                                                              \
            build(hsp::read::csv_object(f));                           \
        }                                                              \
        catch (const std::exception &e)                                \
        {                                                              \
            hsp::msg::error("While reading " + section_name_ + ":\n" + \
                            std::string(e.what()));                    \
        }                                                              \
    }                                                                  \
    void MCMC::read_##section_name(const std::string &f)               \
    {                                                                  \
        pImpl_->section_name.read(f);                                  \
    }

MAKE_READ_METHOD(Setup, setup, csv_list)
MAKE_READ_METHOD(InfecProfile, infec_profile, csv_list)
MAKE_READ_METHOD(Parameter, param, csv_table)
MAKE_READ_METHOD(DataWide, data_wide, csv_table)
MAKE_READ_METHOD(Multiplicity, multiplicity, csv_table)
MAKE_READ_METHOD(DataLong, data_long, csv_table)
MAKE_READ_METHOD(GroupRates, group_rates, csv_table)
MAKE_READ_METHOD(Augmented, augm, csv_table)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize the output file with inferred parameters with column names and initial values.
void MCMC::Impl::print_infer_state(const std::string &output_filename) const
{
    if (hsp::na::is(output_filename))
        return;

    std::ofstream ofs;
    ofs.open(output_filename.c_str());

    if (!ofs.is_open())
        throw std::runtime_error(
            "Unable to open file: " + output_filename);

    // Headers
    // header log-likelihood
    ofs << "iteration" << hsp::read::DELIMITER << "log_lik";

    // header parameters
    for (const std::string &p : param.list_sample)
        ofs << hsp::read::DELIMITER << p;

    // header augmented variable of id trace
    if (augm.is_any_draw())
        for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
            for (const std::string &a : augm.list_draws)
                ofs << hsp::read::DELIMITER << augm.smpl_individuals[k]->id() + "_" + a;

    // header chain
    ofs << hsp::read::DELIMITER << hsp::read::CHAIN_COLUMN << std::endl;

    // Values
    for (size_t c = 0; c < setup.n_chains; ++c)
        for (size_t k = 0; k < setup.c_sample_size; k++)
        {
            size_t i = k + c * setup.c_sample_size;

            // print iteration
            unsigned long sweep = setup.n_burnin + (k + 1) * setup.n_thin;
            ofs << sweep;

            // print log-likelihood
            ofs << hsp::read::DELIMITER << std::setprecision(2)
                << std::fixed << smpl_total_log_lik_[i];

            // print parameter values
            for (const std::string &p : param.list_sample)
                ofs << hsp::read::DELIMITER << std::setprecision(5)
                    << std::fixed << param.smpl_values.find(p)->second[i];

            // print augmented variable of id trace
            if (!augm.smpl_values.empty())
                for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                    for (const std::string &a : augm.list_draws)
                        ofs << hsp::read::DELIMITER << std::setprecision(5)
                            << std::fixed << augm.smpl_values[k].find(a)->second[i];

            // print which chain
            ofs << hsp::read::DELIMITER << (c + 1) << std::endl;
        }

    ofs.close();
}
void MCMC::print_infer_state(const std::string &f) const { pImpl_->print_infer_state(f); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : print simulation  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Print all augmented variables of each individual.
std::ofstream open_ofstream_datawide(const std::string &f,
                                     const str_vec &data_headers,
                                     const str_vec &augm_headers)
{
    std::ofstream ofs;
    ofs.open(f.c_str());

    if (!ofs.is_open())
        throw std::runtime_error("Unable to open file: '" + f + "'");

    // print default data headers
    ofs << "id" << hsp::read::DELIMITER << "house_id" << hsp::read::DELIMITER << "status";

    // print other data headers
    for (const std::string &d : data_headers)
        ofs << hsp::read::DELIMITER << d;

    // print augmented variables headers
    for (const std::string &a : augm_headers)
        ofs << hsp::read::DELIMITER << a;

    ofs << std::endl;

    return ofs;
}

/// Print all augmented variables of each individual.
void update_ofstream_datawide(std::ofstream &ofs,
                              Household &house,
                              const str_vec &data_headers,
                              const str_vec &augm_headers)
{
    if (!ofs.is_open())
        throw std::runtime_error("Output file closed unexpectedly.");

    for (const Individual &member : house)
    {
        ofs << member.id()
            << hsp::read::DELIMITER << house.house_id()
            << hsp::read::DELIMITER << member.status();
        for (const std::string &d : data_headers)
            ofs << hsp::read::DELIMITER << member.data(d);

        for (const std::string &a : augm_headers)
            if (member.is_infected_status())
                ofs << hsp::read::DELIMITER << member.augm(a);
            else
                ofs << hsp::read::DELIMITER << hsp::na::DOUBLE;

        ofs << std::endl;
    }
}

/// Print the simulated data to the provided file.
void MCMC::Impl::print_simulated_data_wide(const std::string &pathtofile)
{
    try
    {
        // initialization
        setup.prog_declare("simulate");
        std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                                   data_wide.headers(),
                                                   augm.list_export);
        aux_init_simulate_();
        // open progress bar
        setup.prog_open(data_wide.report_house.total);

        for (Household house : data_wide) // NOTE: deep copy of household
            try
            {
                // inplace simulation
                simulate_single_household_(house);
                // print the simulated data
                update_ofstream_datawide(ofs, house,
                                         data_wide.headers(),
                                         augm.list_export);
                // update progress bar
                setup.prog_update();
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating " + house.err_name() + ":\n" + e.what());
            }

        // close progress bar
        setup.prog_close();

        // close ofstream
        if (!ofs.is_open())
            throw std::runtime_error(
                "Output file '" + pathtofile + "' closed unexpectedly.");
        ofs.close();
    }
    catch (const std::exception &e)
    {
        setup.prog_crash();
        hsp::msg::error(e.what());
    }
}
void MCMC::print_simulated_data_wide(const std::string &f) { pImpl_->print_simulated_data_wide(f); }
#endif
