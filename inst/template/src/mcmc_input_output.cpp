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
 * mcmc_i_o_rcpp.cpp
 *
 */

#include "tools/error.h"
#include "tools/read.h"
#include "mcmc_impl.h"

#ifdef _HOUSEPI_RCPP

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
    u_dict<std::string> output = {};

    if (!Rf_isNull(input))
    {
        List list(input);
        CharacterVector rownames = list.names();

        for (R_xlen_t i = 0; i < list.size(); ++i)
        {
            CharacterVector element = list[i];
            std::string value = as<std::string>(wrap(element[0]));
            output[as<std::string>(rownames[i])] = value;
        }
    }

    return output;
}

/// Convert R dataframe into unordered_map<string, vector<string>>.
u_dict<str_vec> rdataframe_to_map(SEXP input)
{
    using namespace Rcpp;
    u_dict<str_vec> output = {};

    if (Rf_isNull(input))
        return output;

    DataFrame df(input);
    CharacterVector colnames = df.names();

    for (R_xlen_t i = 0; i < df.size(); ++i)
    {
        CharacterVector col_str = as<CharacterVector>(df[i]);
        str_vec col_data(col_str.begin(), col_str.end());
        output[as<std::string>(colnames[i])] = std::move(col_data);
    }

    return output;
}

#define MAKE_RCPP_IMPORT_METHOD(ImplSection, r_to_map)     \
    void MCMC::Impl::ImplSection::rcpp_import(SEXP input)  \
    {                                                      \
        try                                                \
        {                                                  \
            build(r_to_map(input));                        \
        }                                                  \
        catch (const std::exception &e)                    \
        {                                                  \
            hsp::msg::error(                               \
                "While reading " + section_name_ + ":\n" + \
                std::string(e.what()));                    \
        }                                                  \
    }

#define MAKE_RCPP_IMPORT_LIST_METHOD(ImplSection) MAKE_RCPP_IMPORT_METHOD(ImplSection, rlist_to_map)
MAKE_RCPP_IMPORT_LIST_METHOD(Setup)
MAKE_RCPP_IMPORT_LIST_METHOD(InfecProfile)

#define MAKE_RCPP_IMPORT_TABLE_METHOD(ImplSection) MAKE_RCPP_IMPORT_METHOD(ImplSection, rdataframe_to_map)
MAKE_RCPP_IMPORT_TABLE_METHOD(Parameter);
MAKE_RCPP_IMPORT_TABLE_METHOD(DataWide);
MAKE_RCPP_IMPORT_TABLE_METHOD(DataLong);
MAKE_RCPP_IMPORT_TABLE_METHOD(GroupRates);

/// Import multiplicity definitions from R DataFrame and run build.
void MCMC::Impl::DataWide::rcpp_import_multiplicity(SEXP input)
{
    try
    {
        const u_dict<str_vec> deftable = rdataframe_to_map(input);
        build_multiplicity(deftable);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While reading multiplicity:\n" +
                        std::string(e.what()));
    }
}

/// Import data from R DataFrame then run build.
void MCMC::Impl::Augmented::rcpp_import(SEXP input)
{
    try
    {
        deftable = rdataframe_to_map(input);
        build();
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While reading " + section_name_ + ":\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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

/// Export to R definitions of infectivity profile.
Rcpp::List MCMC::Impl::InfecProfile::rcpp_export() const
{
    using namespace Rcpp;

    // FIXME: this is not working with new argument encoding

    if (profile_ != nullptr)
        return List::create(_["name"] = profile_->name(),
                            _["discr"] = profile_->is_discrete(),
                            _["args"] = profile_->rcpp_export_args(),
                            _["range"] = profile_->rcpp_export_range());
    return List::create();
}

/// Export to R all augmented variables at the current state.
Rcpp::DataFrame MCMC::Impl::Augmented::rcpp_export_data_wide() const
{
    using namespace Rcpp;

    const size_t n_cols = pImpl_.dataset.n_individuals();

    int_vec status_col(n_cols);
    std::vector<dbl_vec> augm_cols(list_export.size());
    for (size_t i = 0; i < list_export.size(); ++i)
        augm_cols[i] = dbl_vec(n_cols);

    size_t row = 0;
    for (Household &house : pImpl_.dataset)
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

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : Impl : dataset exports  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_export(const str_vec &columns, bool is_status) const
{
    using namespace Rcpp;

    try
    {
        if (is_empty())
            return DataFrame();

        const size_t n_rows = n_individuals_;
        CharacterVector df_id = CharacterVector(n_rows);
        CharacterVector df_house_id = CharacterVector(n_rows);
        IntegerVector df_status = IntegerVector(n_rows);

        // other non standard dataset columns
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
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : DataLong : export datalong  //////////////////////////
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
        for (size_t i = 0; i < pImpl_.dataset.n_individuals(); ++i)
            n_rows += indiv_datalist_[i].at(headers_[0]).size();

        CharacterVector datalist_index = CharacterVector(n_rows);

        // other non standard datalong columns
        std::vector<NumericVector> datalist_dbl = {};
        datalist_dbl.reserve(headers_.size());
        for (size_t c = 0; c < headers_.size(); ++c)
            datalist_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : pImpl_.dataset)
            for (const Individual &member : house)
                for (size_t k = 0; k < member.datalong(headers_[0]).size(); ++k, ++r)
                {
                    datalist_index[r] = member.id();
                    for (size_t c = 0; c < headers_.size(); ++c)
                        datalist_dbl[c][r] = hsp_to_rcpp(member.datalong(headers_[c])[k]);
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
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : report exports  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_summarise_cases() const
{
    using namespace Rcpp;

    try
    {
        if (is_empty())
            hsp::msg::warning("Calling a summarise_cases on empty data.");

        CharacterVector type = wrap(str_vec(4, "individual") + str_vec(4, "household"));
        const size_t n_rows = type.size();
        CharacterVector comment = {"total", "missing", "healthy", "infected",
                                   "total", "no cases", "at least 1 case", "at least 2 cases"};
        IntegerVector counts(n_rows);
        NumericVector pct(n_rows);

        counts[0] = report_indiv.total;
        pct[0] = 100.;

        counts[1] = report_indiv.missing;
        pct[1] = 100. * static_cast<double>(report_indiv.missing) / report_indiv.total;

        counts[2] = report_indiv.healthy;
        pct[2] = 100. * static_cast<double>(report_indiv.healthy) / report_indiv.total;

        counts[3] = report_indiv.infected;
        pct[3] = 100. * static_cast<double>(report_indiv.infected) / report_indiv.total;

        counts[4] = report_house.total;
        pct[4] = 100.;
        counts[5] = report_house.no_cases;
        pct[5] = 100. * static_cast<double>(report_house.no_cases) / report_house.total;

        counts[6] = report_house.at_least_1_case;
        pct[6] = 100. * static_cast<double>(report_house.at_least_1_case) / report_house.total;

        counts[7] = report_house.at_least_2_cases;
        pct[7] = 100. * static_cast<double>(report_house.at_least_2_cases) / report_house.total;

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
    return DataFrame();
}

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
        pct[0] = 100. * static_cast<double>(report_hh_size.max_cases.second) / report_indiv.total;

        value[1] = report_hh_size.max_miss.first;
        counts[1] = report_hh_size.max_miss.second;
        pct[1] = 100. * static_cast<double>(report_hh_size.max_miss.second) / report_indiv.total;

        value[2] = report_hh_size.min_size.first;
        counts[2] = report_hh_size.min_size.second;
        pct[2] = 100. * static_cast<double>(report_hh_size.min_size.second) / report_indiv.total;

        value[3] = report_hh_size.max_size.first;
        counts[3] = report_hh_size.max_size.second;
        pct[3] = 100. * static_cast<double>(report_hh_size.max_size.second) / report_indiv.total;

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
    return DataFrame();
}

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
        if (smpl_full_loglike_.empty())
            throw std::length_error("No draws to export were found.");

        // cllect the column names
        str_vec headers = {"loglikelihood"};

        headers += param.list_sample;

        // augmented variable of id trace
        if (!augm.smpl_values.empty())
            for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                for (const std::string &a : augm.list_draws)
                    headers.push_back(augm.smpl_individuals[k]->id() + "_" + a);

        // multiple chains
        headers.push_back(hsp::read::CHAIN_COLUMN);

        // matrix with draws
        NumericMatrix output(setup.n_chains * setup.c_sample_size,
                             headers.size());
        size_t col = 0;

        // loglikelihood
        output(_, col++) = NumericVector(wrap(std::move(smpl_full_loglike_)));

        // parameters
        for (const std::string &p : param.list_sample)
            output(_, col++) = NumericVector(wrap(std::move(param.smpl_values.at(p))));

        // augmented variables
        if (!augm.smpl_values.empty())
            for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
            {
                const auto &smpl_augm_value_k = augm.smpl_values[k];
                for (const std::string &a : augm.list_draws)
                    output(_, col++) = NumericVector(wrap(std::move(smpl_augm_value_k.at(a))));
            }

        // asssociate chain tag
        output(_, col++) = setup.draws_chain_tags();
        // associate column names
        colnames(output) = CharacterVector(wrap(std::move(headers)));
        // associate row names
        rownames(output) = setup.draws_iter_names();

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting MCMC draws:\n" +
                        std::string(e.what()));
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  MCMC : Impl : export household loglikelihood matrix  //////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export all sampled values of the loglikelihood for all household.
Rcpp::NumericMatrix MCMC::Impl::DataWide::rcpp_export_draws() const
{
    using namespace Rcpp;

    try
    {
        if (smpl_values_.empty())
            throw std::length_error("No draws to export were found.");

        // create a NumericMatrix with n_house columns
        const unsigned long sample_size = pImpl_.setup.n_chains * pImpl_.setup.c_sample_size;
        NumericMatrix output(sample_size, n_houses() + 1);
        CharacterVector columns(n_houses() + 1);
        columns[n_houses()] = hsp::read::CHAIN_COLUMN;

        // retrieve restricted index by general index of multiplicity head class
        std::unordered_map<size_t, size_t> restr_idx_map = {};
        restr_idx_map.reserve(restrict_indexes_.size());
        for (size_t R = 0; R < restrict_indexes_.size(); ++R)
            restr_idx_map.emplace(restrict_indexes_[R], R);

        // assign each NumericVector to a column of the matrix
        for (size_t H = 0; H < n_houses(); ++H)
        {
            size_t R = restr_idx_map.at(house_vector_[H].multiplicity_class_index());
            output(_, H) = NumericVector(wrap(std::move(smpl_values_[R])));
            columns[H] = house_vector_[H].house_id();
        }

        // asssociate chain tag
        output(_, n_houses()) = pImpl_.setup.draws_chain_tags();
        // associate column names
        colnames(output) = columns;
        // associate row names
        rownames(output) = pImpl_.setup.draws_iter_names();

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting household loglikelihood draws:\n" +
                        std::string(e.what()));
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : house state snapshot //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export to R a snapshot of current state of house dataset.
Rcpp::NumericVector MCMC::Impl::Parameter::rcpp_snapshot() const
{
    using namespace Rcpp;

    if (!dictionary_.empty())
    {
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

    return NumericVector::create();
}

/// Export to R a snapshot of current state of house dataset.
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    const Household &house = this->operator()(house_id);
    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);
    CharacterVector df_house_id = CharacterVector(n_rows);
    IntegerVector df_status = IntegerVector(n_rows);

    std::vector<NumericVector> df_double;
    df_double.reserve(headers_.size());
    for (size_t c = 0; c < headers_.size(); ++c)
        df_double.emplace_back(NumericVector(n_rows));

    size_t r = 0;
    for (const Individual &member : house)
    {
        df_id[r] = member.id();
        df_house_id[r] = house.house_id();
        df_status[r] = member.status();

        for (size_t c = 0; c < headers_.size(); ++c)
            df_double[c][r] = hsp_to_rcpp(member.data(headers_[c]));
        ++r;
    }

    DataFrame output = DataFrame::create(_["id"] = std::move(df_id),
                                         _["house_id"] = std::move(df_house_id),
                                         _["status"] = std::move(df_status));
    for (size_t c = 0; c < headers_.size(); ++c)
        output.push_back(df_double[c], headers_[c]);

    return output;
}

/// Export to R a snapshot of current state of house datalong.
Rcpp::DataFrame MCMC::Impl::DataLong::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (!is_empty())
    {
        const Household &house = pImpl_.dataset(house_id);
        size_t n_rows = 0;

        for (const Individual &member : house)
            n_rows += member.datalong(headers_[0]).size();

        CharacterVector datalist_id = CharacterVector(n_rows);

        std::vector<NumericVector> datalist_dbl = {};
        datalist_dbl.reserve(headers_.size());
        for (size_t i = 0; i < headers_.size(); ++i)
            datalist_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Individual &member : house)
            for (size_t i = 0; i < member.datalong(headers_[0]).size(); ++i)
            {
                datalist_id[r] = member.id();
                for (size_t c = 0; c < headers_.size(); ++c)
                    datalist_dbl[c][r] = hsp_to_rcpp(member.datalong(headers_[c])[i]);
                ++r;
            }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(_["id"] = std::move(datalist_id));
        for (size_t c = 0; c < headers_.size(); ++c)
            output.push_back(datalist_dbl[c], headers_[c]);

        return output;
    }
    return DataFrame::create();
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::Impl::GroupRates::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (!headers_.empty())
    {
        const Household &house = pImpl_.dataset(house_id);
        const size_t n_rows = house.size();
        CharacterVector df_id = CharacterVector(n_rows);

        std::vector<NumericVector> group_rate_dbl;
        group_rate_dbl.reserve(headers_.size());
        for (size_t c = 0; c < headers_.size(); ++c)
            group_rate_dbl.emplace_back(NumericVector(n_rows));

        size_t r = 0;
        for (const Individual &member : house)
        {
            df_id[r] = member.id();

            for (size_t c = 0; c < headers_.size(); ++c)
                group_rate_dbl[c][r] = hsp_to_rcpp(member.group_rate(headers_[c]));
            ++r;
        }

        DataFrame output = DataFrame::create(_["id"] = df_id);
        for (size_t c = 0; c < headers_.size(); ++c)
            output.push_back(group_rate_dbl[c], headers_[c]);

        return output;
    }
    return DataFrame::create();
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::Impl::Augmented::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (!list_export.empty())
    {
        const Household &house = pImpl_.dataset(house_id);
        const size_t n_rows = house.size();
        CharacterVector df_id = CharacterVector(n_rows);

        std::vector<NumericVector> augm_dbl;
        augm_dbl.reserve(list_export.size());
        for (size_t c = 0; c < list_export.size(); ++c)
            augm_dbl.emplace_back(NumericVector(n_rows));

        size_t r = 0;
        for (const Individual &member : house)
        {
            df_id[r] = member.id();
            for (size_t c = 0; c < list_export.size(); ++c)
                augm_dbl[c][r] = hsp_to_rcpp(member.augm(list_export[c]));
            ++r;
        }

        DataFrame output = DataFrame::create(_["id"] = df_id);
        for (size_t c = 0; c < list_export.size(); ++c)
            output.push_back(augm_dbl[c], list_export[c]);

        return output;
    }
    return DataFrame::create();
}

/// Export to R a snapshot of current state of house infectivity profile.
Rcpp::List MCMC::Impl::InfecProfile::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (is_augmented() || (profile_ != nullptr))
    {
        const Household &house = pImpl_.dataset(house_id);
        const size_t n_rows = house.size();
        CharacterVector df_id = CharacterVector(n_rows);

        CharacterVector infprof_name(n_rows, NA_STRING);
        LogicalVector infprof_discr(n_rows, NA_LOGICAL);
        List infprof_range(n_rows);
        List infprof_args(n_rows);

        size_t r = 0;
        for (const Individual &member : house)
        {
            df_id[r] = member.id();
            if (member.is_infected_status())
            {
                if (is_augmented())
                {
                    const Variable &var = member.get_augm_variable(hsp::infprof::INSTANTANEOUS);
                    infprof_name[r] = var.func_name();
                    infprof_discr[r] = var.is_discrete();
                    infprof_args[r] = var.rcpp_export_func_args();
                    infprof_range[r] = var.rcpp_export_range();
                }
                else
                {
                    infprof_name[r] = profile_->name();
                    infprof_discr[r] = profile_->is_discrete();
                    infprof_args[r] = profile_->rcpp_export_args();
                    infprof_range[r] = profile_->rcpp_export_range();
                }
            }
            else
            {
                infprof_args[r] = NumericVector(0);
                infprof_range[r] = NumericVector(0);
            }
            ++r;
        }

        return List::create(_["id"] = df_id,
                            _["discr"] = infprof_discr,
                            _["name"] = infprof_name,
                            _["args"] = infprof_args,
                            _["range"] = infprof_range,
                            _["stringsAsFactors"] = false);
    }
    return List::create();
}

/// Export a list containing all relevant information for a debug.
Rcpp::List MCMC::Impl::rcpp_export_house_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    try
    {
        if (!dataset.is_empty())
            return List::create(_["param"] = param.rcpp_snapshot(),
                                _["data_wide"] = dataset.rcpp_snapshot(house_id),
                                _["data_long"] = datalong.rcpp_snapshot(house_id),
                                _["group_rates"] = group_rates.rcpp_snapshot(house_id),
                                _["infec_profile"] = infec_profile.rcpp_snapshot(house_id),
                                _["augm"] = augm.rcpp_snapshot(house_id));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting house \"" + house_id + "\" for debug:\n" +
                        std::string(e.what()));
    }
    return List::create();
}

#else

#include <fstream>
#include <iomanip>

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : read  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_READ_METHOD(ImplSection, csv_object)                      \
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
    }

#define MAKE_READ_LIST_METHOD(ImplSection) MAKE_READ_METHOD(ImplSection, csv_list)
MAKE_READ_LIST_METHOD(Setup)
MAKE_READ_LIST_METHOD(InfecProfile)

#define MAKE_READ_TABLE_METHOD(ImplSection) MAKE_READ_METHOD(ImplSection, csv_table)
MAKE_READ_TABLE_METHOD(Parameter);
MAKE_READ_TABLE_METHOD(DataWide);
MAKE_READ_TABLE_METHOD(DataLong);
MAKE_READ_TABLE_METHOD(GroupRates);

/// Read multiplicity file and then run assignment.
void MCMC::Impl::DataWide::read_multiplicity(const std::string &f)
{
    try
    {
        const u_dict<str_vec> deftable = hsp::read::csv_table(f);
        build_multiplicity(deftable);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While reading multiplicity:\n" +
                        std::string(e.what()));
    }
}

/// Read Model Augmented Variable file.
void MCMC::Impl::Augmented::read(const std::string &f)
{
    try
    {
        // UPDATE: deftable should not be stored internally in this format
        deftable = hsp::read::csv_table(f);
        build();
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While reading " + section_name_ + ":\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize the output file with inferred parameters with column names and initial values.
void MCMC::Impl::print_infer_state(const std::string &output_filename) const
{
    if (!hsp::na::is(output_filename))
    {
        std::ofstream ofs;
        ofs.open(output_filename.c_str());

        if (ofs.is_open())
        {
            // Headers
            // header loglikelihood
            ofs << "iteration" << hsp::read::DELIMITER << "loglikelihood";

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

                    // print loglikelihood
                    ofs << hsp::read::DELIMITER << std::setprecision(2)
                        << std::fixed << smpl_full_loglike_[i];

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
        else
            throw std::runtime_error(
                "Unable to open file: " + output_filename);
    }
}

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

    if (ofs.is_open())
    {
        // print default data headers
        ofs << "id" << hsp::read::DELIMITER << "house_id" << hsp::read::DELIMITER << "status";

        // print other data headers
        for (const std::string &d : data_headers)
            ofs << hsp::read::DELIMITER << d;

        // print augmented variables headers
        for (const std::string &a : augm_headers)
            ofs << hsp::read::DELIMITER << a;

        ofs << std::endl;
    }
    else
        throw std::runtime_error("Unable to open file: \"" + f + "\"");

    return ofs;
}

/// Print all augmented variables of each individual.
void update_ofstream_datawide(std::ofstream &ofs,
                              Household &house,
                              const str_vec &data_headers,
                              const str_vec &augm_headers)
{
    if (ofs.is_open())
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
    else
        throw std::runtime_error("Output file closed unexpectedly.");
}

/// Print the simulated dataset to the provided file.
void MCMC::Impl::print_simulated_data_wide(const std::string &pathtofile)
{
    try
    {
        // initialization
        setup.prog_declare("simulate");
        std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                                   dataset.headers(),
                                                   augm.list_export);
        aux_init_simulate_();
        // open progress bar
        setup.prog_open(dataset.n_houses());

        for (Household house : dataset) // NOTE: deep copy of household
            try
            {
                // inplace simulation
                simulate_single_household_(house);
                // print the simulated dataset
                update_ofstream_datawide(ofs, house,
                                         dataset.headers(),
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
        if (ofs.is_open())
            ofs.close();
        else
            throw std::runtime_error(
                "Output file \"" + pathtofile + "\" closed unexpectedly.");
    }
    catch (const std::exception &e)
    {
        setup.prog_crash();
        hsp::msg::error(e.what());
    }
}

#endif
