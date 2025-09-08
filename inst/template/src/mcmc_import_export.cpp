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
 * mcmc_import_export.cpp
 *
 */

#include "mcmc_impl.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifdef _HOUSEPI_RCPP
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : imports from Rcpp /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define DEFINE_RCPP_IMPORT_METHOD(ImplSection, r_to_map)                 \
    void MCMC::Impl::ImplSection::rcpp_import(SEXP input)                \
    {                                                                    \
        try                                                              \
        {                                                                \
            build(tools::r_to_map(input));                               \
        }                                                                \
        catch (const std::exception &e)                                  \
        {                                                                \
            tools::error("While initializing " + section_name_ + ":\n" + \
                         std::string(e.what()));                         \
        }                                                                \
    }

#define DEFINE_RCPP_IMPORT_LIST_METHOD(ImplSection) DEFINE_RCPP_IMPORT_METHOD(ImplSection, rlist_to_map)
DEFINE_RCPP_IMPORT_LIST_METHOD(Setup)
DEFINE_RCPP_IMPORT_LIST_METHOD(InfecProfile)

#define DEFINE_RCPP_IMPORT_TABLE_METHOD(ImplSection) DEFINE_RCPP_IMPORT_METHOD(ImplSection, rdataframe_to_map)
DEFINE_RCPP_IMPORT_TABLE_METHOD(Parameter);
DEFINE_RCPP_IMPORT_TABLE_METHOD(DataWide);
DEFINE_RCPP_IMPORT_TABLE_METHOD(DataLong);
DEFINE_RCPP_IMPORT_TABLE_METHOD(GroupRates);

/* Import multiplicity definitions from R DataFrame and run build. */
void MCMC::Impl::DataWide::rcpp_import_multiplicity(SEXP input)
{
    try
    {
        const dict_str_vec deftable = tools::rdataframe_to_map(input);
        build_multiplicity(deftable);
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing multiplicity:\n" +
                     std::string(e.what()));
    }
}

/* Import data from R DataFrame then run build. */
void MCMC::Impl::Augmented::rcpp_import(SEXP input)
{
    try
    {
        // UPDATE: deftable should not be stored internally in this format
        deftable = tools::rdataframe_to_map(input);
        build();
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing " + section_name_ + ":\n" +
                     std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current setup configurations to R List. */
Rcpp::List MCMC::Impl::Setup::rcpp_export() const
{
    using namespace Rcpp;

    return List::create(_["n_chains"] = n_chains,
                        _["sample_size"] = c_sample_size,
                        _["n_steps"] = n_steps,
                        _["n_adapt"] = n_adapt,
                        _["n_burnin"] = n_burnin,
                        _["thin"] = thin,
                        _["rand_init"] = rand_init_params,
                        _["seed"] = seed,
                        _["n_threads"] = n_threads,
                        _["id_draws"] = CharacterVector(wrap(id_draws)),
                        _["progress"] = progress);
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::Parameter::rcpp_export() const
{
    using namespace Rcpp;

    try
    {
        str_vec V_names = list_sample + list_const;

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
            const Variable &this_par = dictionary_.find(V_names[r])->second;
            V_value[r] = this_par.get();
            V_constant[r] = this_par.is_constant();
            V_discrete[r] = this_par.is_discrete();

            V_prior[r] = tools::stoRString(this_par.distr_name());
            V_arg_name[r] = tools::stoRString(this_par.distr_str_arg_names());
            V_arg_value[r] = tools::stoRString(this_par.distr_str_arg_values());

            V_rule[r] = tools::stoRString(this_par.prop_name());
            const std::string rule_scale = this_par.prop_str_scale();
            V_rule_scale[r] = tools::is_na(rule_scale) ? NA_REAL : tools::str::tod(rule_scale);

            const std::string min_val = this_par.str_min();
            V_min[r] = tools::is_na(min_val) ? NA_REAL : tools::str::tod(min_val);

            const std::string max_val = this_par.str_max();
            V_max[r] = tools::is_na(max_val) ? NA_REAL : tools::str::tod(max_val);
        }

        return DataFrame::create(_["name"] = CharacterVector(wrap(std::move(V_names))),
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
        tools::error("While exporting parameter definitions:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

/* Export to R definitions of infectivity profile. */
Rcpp::List MCMC::Impl::InfecProfile::rcpp_export() const
{
    using namespace Rcpp;

    if (profile_ != nullptr)
        return List::create(_["name"] = profile_->name(),
                            _["discr"] = profile_->is_discrete(),
                            _["args"] = profile_->rcpp_export_args(),
                            _["range"] = profile_->rcpp_export_range());
    return List::create();
}

/* Export to R all augmented variables at the current state. */
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
                augm_cols[i][row] = member.augm(list_export[i]);
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
        CharacterVector dataset_id = CharacterVector(n_rows);
        CharacterVector dataset_house_id = CharacterVector(n_rows);
        IntegerVector dataset_status = IntegerVector(n_rows);

        // other non standard dataset columns
        std::vector<NumericVector> dataset_dbl;
        dataset_dbl.reserve(headers_.size());
        for (size_t c = 0; c < headers_.size(); ++c)
            dataset_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : house_vector_)
            for (const Individual &member : house)
            {
                dataset_id[r] = member.id();
                dataset_house_id[r] = house.house_id();
                dataset_status[r] = member.status();

                for (size_t c = 0; c < columns.size(); ++c)
                    dataset_dbl[c][r] = member.data(columns[c]);
                ++r;
            }

        DataFrame output = DataFrame::create(_["id"] = dataset_id,
                                             _["house_id"] = dataset_house_id);
        if (is_status)
            output.push_back(dataset_status, "status");
        for (size_t c = 0; c < columns.size(); ++c)
            output.push_back(dataset_dbl[c], columns[c]);

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting data (wide):\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC : Impl  : export datalong  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
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
                        datalist_dbl[c][r] = member.datalong(headers_[c])[k];
                }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(_["id"] = datalist_index);
        for (size_t c = 0; c < headers_.size(); ++c)
            output.push_back(datalist_dbl[c], headers_[c]);

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting data (long):\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : report exports  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_summarise_cases() const
{
    using namespace Rcpp;

    try
    {
        if (is_empty())
            throw std::length_error(
                "No report available for empty data.");

        size_t n_rows = report_indiv_.size() + report_house_.size();
        CharacterVector type(n_rows), comment(n_rows);
        IntegerVector counts(n_rows);
        NumericVector pct(n_rows);

        size_t k = 0;
        for (const auto &pair : report_indiv_)
        {
            type[k] = "individual";
            comment[k] = pair.first;
            counts[k] = pair.second;
            pct[k] = 100. * static_cast<double>(pair.second) / n_individuals();
            ++k;
        }
        for (const auto &pair : report_house_)
        {
            type[k] = "household";
            comment[k] = pair.first;
            counts[k] = pair.second;
            pct[k] = 100. * static_cast<double>(pair.second) / n_houses();
            ++k;
        }

        return DataFrame::create(_["type"] = std::move(type),
                                 _["comment"] = std::move(comment),
                                 _["counts"] = std::move(counts),
                                 _["pct"] = std::move(pct));
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting data report:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : export house size report  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_summarise_house_sizes() const
{
    using namespace Rcpp;

    try
    {
        size_t n_rows = report_size_.size();
        CharacterVector comment(n_rows);
        IntegerVector value(n_rows), counts(n_rows);
        NumericVector pct(n_rows);

        size_t k = 0;
        for (const auto &pair : report_size_)
        {
            comment[k] = pair.first;
            value[k] = pair.second.first;
            counts[k] = pair.second.second;
            pct[k] = 100. * static_cast<double>(pair.second.second) / n_individuals_;
            ++k;
        }

        return DataFrame::create(_["comment"] = std::move(comment),
                                 _["value"] = std::move(value),
                                 _["counts"] = std::move(counts),
                                 _["pct"] = std::move(pct));
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting household size report:\n" +
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
            iteration[k + c * c_sample_size] = n_burnin + (k + 1) * thin;
    return CharacterVector(wrap(std::move(iteration)));
}

/* Export MCMC draws of parameter posterior samples. */
Rcpp::NumericMatrix MCMC::Impl::rcpp_export_draws() const
{
    using namespace Rcpp;

    try
    {
        // UPDATE: I shall allow for export even if _is_power_post
        if (smpl_full_loglike_.empty())
            throw std::length_error("No draws to export were found.");

        //----- Collect the column names -----//
        str_vec headers = {"loglikelihood"};

        headers += param.list_sample;

        // augmented variable of id trace
        if (!augm.smpl_values.empty())
            for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                for (const std::string &a : augm.list_draws)
                    headers.push_back(augm.smpl_individuals[k]->id() + "_" + a);

        // multiple chains
        headers.push_back(HSP_CHAIN_COL);

        //----- NumericMatrix with draws -----//
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
        tools::error("While exporting MCMC draws:\n" +
                     std::string(e.what()));
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  MCMC : Impl : export household loglikelihood matrix  //////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export all sampled values of the loglikelihood for all household. */
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
        columns[n_houses()] = HSP_CHAIN_COL;

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
        tools::error("While exporting household loglikelihood draws:\n" +
                     std::string(e.what()));
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : house state snapshot //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export to R a snapshot of current state of house dataset. */
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

/* Export to R a snapshot of current state of house dataset. */
Rcpp::DataFrame MCMC::Impl::DataWide::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    const Household &house = this->operator()(house_id);
    const size_t n_rows = house.size();
    CharacterVector dataset_id = CharacterVector(n_rows);
    CharacterVector dataset_house_id = CharacterVector(n_rows);
    IntegerVector dataset_status = IntegerVector(n_rows);

    std::vector<NumericVector> dataset_dbl;
    dataset_dbl.reserve(headers_.size());
    for (size_t c = 0; c < headers_.size(); ++c)
        dataset_dbl.emplace_back(NumericVector(n_rows));

    size_t r = 0;
    for (const Individual &member : house)
    {
        dataset_id[r] = member.id();
        dataset_house_id[r] = house.house_id();
        dataset_status[r] = member.status();

        for (size_t c = 0; c < headers_.size(); ++c)
            dataset_dbl[c][r] = member.data(headers_[c]);
        ++r;
    }

    DataFrame output = DataFrame::create(_["id"] = std::move(dataset_id),
                                         _["house_id"] = std::move(dataset_house_id),
                                         _["status"] = std::move(dataset_status));
    for (size_t c = 0; c < headers_.size(); ++c)
        output.push_back(dataset_dbl[c], headers_[c]);

    return output;
}

/* Export to R a snapshot of current state of house datalong. */
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
                    datalist_dbl[c][r] = member.datalong(headers_[c])[i];
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

/* Export to R a snapshot of current state of house augmented variables. */
Rcpp::DataFrame MCMC::Impl::GroupRates::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (!is_empty())
    {
        const Household &house = pImpl_.dataset(house_id);
        const size_t n_rows = house.size();
        CharacterVector dataset_id = CharacterVector(n_rows);

        std::vector<NumericVector> group_rate_dbl;
        group_rate_dbl.reserve(headers_.size());
        for (size_t c = 0; c < headers_.size(); ++c)
            group_rate_dbl.emplace_back(NumericVector(n_rows));

        size_t r = 0;
        for (const Individual &member : house)
        {
            dataset_id[r] = member.id();

            for (size_t c = 0; c < headers_.size(); ++c)
                group_rate_dbl[c][r] = member.group_rate(headers_[c]);
            ++r;
        }

        DataFrame output = DataFrame::create(_["id"] = dataset_id);
        for (size_t c = 0; c < headers_.size(); ++c)
            output.push_back(group_rate_dbl[c], headers_[c]);

        return output;
    }
    return DataFrame::create();
}

/* Export to R a snapshot of current state of house augmented variables. */
Rcpp::DataFrame MCMC::Impl::Augmented::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (!list_export.empty())
    {
        const Household &house = pImpl_.dataset(house_id);
        const size_t n_rows = house.size();
        CharacterVector dataset_id = CharacterVector(n_rows);

        std::vector<NumericVector> augm_dbl;
        augm_dbl.reserve(list_export.size());
        for (size_t c = 0; c < list_export.size(); ++c)
            augm_dbl.emplace_back(NumericVector(n_rows));

        size_t r = 0;
        for (const Individual &member : house)
        {
            dataset_id[r] = member.id();
            for (size_t c = 0; c < list_export.size(); ++c)
                augm_dbl[c][r] = member.augm(list_export[c]);
            ++r;
        }

        DataFrame output = DataFrame::create(_["id"] = dataset_id);
        for (size_t c = 0; c < list_export.size(); ++c)
            output.push_back(augm_dbl[c], list_export[c]);

        return output;
    }
    return DataFrame::create();
}

/* Export to R a snapshot of current state of house infectivity profile. */
Rcpp::List MCMC::Impl::InfecProfile::rcpp_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    if (is_augmented() || (profile_ != nullptr))
    {
        const Household &house = pImpl_.dataset(house_id);
        const size_t n_rows = house.size();
        CharacterVector dataset_id = CharacterVector(n_rows);

        CharacterVector infprof_name(n_rows, NA_STRING);
        LogicalVector infprof_discr(n_rows, NA_LOGICAL);
        List infprof_range(n_rows);
        List infprof_args(n_rows);

        size_t r = 0;
        for (const Individual &member : house)
        {
            dataset_id[r] = member.id();
            if (member.is_infected_status())
            {
                if (is_augmented())
                {
                    const Variable &var = member.get_augm_variable(HSP_PROFILE_INST);
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

        return List::create(_["id"] = dataset_id,
                            _["discr"] = infprof_discr,
                            _["name"] = infprof_name,
                            _["args"] = infprof_args,
                            _["range"] = infprof_range,
                            _["stringsAsFactors"] = false);
    }
    return List::create();
}

/* Export a list containing all relevant information for a debug. */
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
        tools::error("While exporting house \"" + house_id + "\" for debug:\n" +
                     std::string(e.what()));
    }
    return List::create();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#else
#define DEFINE_READ_METHOD(ImplSection, read_csv)                        \
    void MCMC::Impl::ImplSection::read(const std::string &f)             \
    {                                                                    \
        try                                                              \
        {                                                                \
            build(tools::read_csv(f));                                   \
        }                                                                \
        catch (const std::exception &e)                                  \
        {                                                                \
            tools::error("While initializing " + section_name_ + ":\n" + \
                         std::string(e.what()));                         \
        }                                                                \
    }

#define DEFINE_READ_LIST_METHOD(ImplSection) DEFINE_READ_METHOD(ImplSection, read_csv_list)
DEFINE_READ_LIST_METHOD(Setup)
DEFINE_READ_LIST_METHOD(InfecProfile)

#define DEFINE_READ_TABLE_METHOD(ImplSection) DEFINE_READ_METHOD(ImplSection, read_csv_table)
DEFINE_READ_TABLE_METHOD(Parameter);
DEFINE_READ_TABLE_METHOD(DataWide);
DEFINE_READ_TABLE_METHOD(DataLong);
DEFINE_READ_TABLE_METHOD(GroupRates);

/* Read multiplicity file and then run assignment. */
void MCMC::Impl::DataWide::read_multiplicity(const std::string &f)
{
    try
    {
        const dict_str_vec deftable = tools::read_csv_table(f);
        build_multiplicity(deftable);
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing multiplicity:\n" +
                     std::string(e.what()));
    }
}

/* Read Model Augmented Variable file. */
void MCMC::Impl::Augmented::read(const std::string &f)
{
    try
    {
        // UPDATE: deftable should not be stored internally in this format
        deftable = tools::read_csv_table(f);
        build();
    }
    catch (const std::exception &e)
    {
        tools::error("While initializing " + section_name_ + ":\n" +
                     std::string(e.what()));
    }
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
const std::string _percent_(int n, int d) { return std::string(std::round(100. * n / d), 1); }

/* Auxiliary function to call percentuage. */
const std::string _trail_percent_(int n, int d) { return std::to_string(n) + " (" + _percent_(n, d) + ")"; }

/* Print the report of the database and of the model group rates. */
void MCMC::Impl::DataWide::show_summary() const
{
    std::cout << "n. of individuals: " << n_individuals_
              << "\n    * n. of cases: "
              << _trail_percent_(get_report_indiv("infected"), n_individuals_)
              << "\n    * n. of missing status : "
              << _trail_percent_(get_report_indiv("missing"), n_individuals_);

    std::cout << "\nn. of households: " << n_houses()
              << "\n    * with no cases: "
              << _trail_percent_(get_report_house("no cases"), n_houses())
              << "\n    * with at least 1 case: "
              << _trail_percent_(get_report_house("at least 1 case"), n_houses())
              << "\n    * with at least 2 cases: "
              << _trail_percent_(get_report_house("at least 2 cases"), n_houses());

    std::cout << "\nmax n. house cases: " << get_report_hh_size("max cases").first
              << " (" << _percent_(get_report_hh_size("max cases").second, n_houses()) << ")"
              << "\nmin n. house members: " << get_report_hh_size("min size").first
              << " (" << _percent_(get_report_hh_size("min size").second, n_houses()) << ")"
              << "\nmax n. house members: " << get_report_hh_size("max size").first
              << " (" << _percent_(get_report_hh_size("max size").second, n_houses()) << ")";

    std::cout << std::endl;
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
void MCMC::Impl::Setup::show() const
{

    std::string id_draws_txt = HSP_NA_STRING;
    if (pImpl_.augm.is_any_draw())
    {
        id_draws_txt = pImpl_.augm.smpl_individuals[0]->id();
        for (size_t k = 1; k < pImpl_.augm.smpl_individuals.size(); ++k)
            id_draws_txt += HSP_LIST_DELIM + pImpl_.augm.smpl_individuals[k]->id();
    }

    // common setup variables
    std::cout << "seed" << HSP_DELIM << seed
              << "\nn_threads" << HSP_DELIM << n_threads
              << "\nn_chains" << HSP_DELIM << n_chains
              << "\nsample_size" << HSP_DELIM << c_sample_size
              << "\nn_steps" << HSP_DELIM << n_steps
              << "\nn_adapt" << HSP_DELIM << n_adapt
              << "\nn_burnin" << HSP_DELIM << n_burnin
              << "\nthin" << HSP_DELIM << thin
              << "\nrand_init" << HSP_DELIM << tools::str::format(rand_init_params)
              << "\nprogress" << HSP_DELIM << progress
              << "\nid_draws" << HSP_DELIM << id_draws_txt
              << "\n";
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void MCMC::Impl::Parameter::show() const
{
    // show the headers
    std::cout << file_headers[0];
    for (size_t k = 1; k < file_headers.size(); ++k)
        std::cout << HSP_DELIM << file_headers[k];
    std::cout << "\n";

    // show parameter definitions
    str_vec all_params = list_sample + list_const;
    for (const std::string &p : all_params)
    {
        const Variable &aux = dictionary_.find(p)->second;
        std::cout << aux.name()
                  << HSP_DELIM << tools::str::format(aux.get())
                  << HSP_DELIM << tools::str::format(aux.is_constant())
                  << HSP_DELIM << tools::str::format(aux.is_discrete())
                  << HSP_DELIM << aux.distr_name()
                  << HSP_DELIM << aux.distr_str_arg_names()
                  << HSP_DELIM << aux.distr_str_arg_values()
                  << HSP_DELIM << aux.prop_name()
                  << HSP_DELIM << aux.prop_str_scale()
                  << HSP_DELIM << aux.str_min()
                  << HSP_DELIM << aux.str_max()
                  << "\n";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : show parameter results  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the results of the paramters. */
void MCMC::Impl::Parameter::show_results() const
{
    std::cout << "name" << HSP_DELIM << "value" << HSP_DELIM << "value_sd";
    for (size_t k = 2; k < file_headers.size(); ++k)
        std::cout << HSP_DELIM << file_headers[k];
    std::cout << "\n";

    for (const std::string &p : list_sample)
    {
        const Variable &this_param = dictionary_.find(p)->second;
        const std::pair<double, double> stats = tools::compute_stats(smpl_values.find(p)->second);

        std::cout << p
                  << HSP_DELIM << stats.first
                  << HSP_DELIM << std::sqrt(stats.second)
                  << HSP_DELIM << tools::str::format(this_param.is_constant())
                  << HSP_DELIM << tools::str::format(this_param.is_discrete())
                  << HSP_DELIM << this_param.distr_name()
                  << HSP_DELIM << this_param.distr_str_arg_names()
                  << HSP_DELIM << this_param.distr_str_arg_values()
                  << HSP_DELIM << this_param.str_min()
                  << HSP_DELIM << this_param.str_max()
                  << HSP_DELIM << this_param.prop_name()
                  << HSP_DELIM << this_param.prop_str_scale()
                  << "\n";
    }
    std::cout << std::endl;
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
            inferstate_ofstream << "iteration" << HSP_DELIM << "loglikelihood";

            // header parameters
            for (const std::string &p : param.list_sample)
                inferstate_ofstream << HSP_DELIM << p;

            // augmented variable of id trace
            if (augm.is_any_draw())
                for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                    for (const std::string &a : augm.list_draws)
                        inferstate_ofstream << augm.smpl_individuals[k]->id() + "_" + a;

            // header chain
            inferstate_ofstream << HSP_DELIM << HSP_CHAIN_COL << std::endl;

            //----- VALUES -----//

            for (size_t c = 0; c < setup.n_chains; ++c)
                for (size_t k = 0; k < setup.c_sample_size; k++)
                {
                    size_t i = k + c * setup.c_sample_size;

                    // print iteration
                    unsigned long sweep = setup.n_burnin + (k + 1) * setup.thin;
                    inferstate_ofstream << sweep;

                    // print loglikelihood
                    inferstate_ofstream << HSP_DELIM << std::setprecision(2)
                                        << std::fixed << smpl_full_loglike_[i];

                    // print parameter values
                    for (const std::string &p : param.list_sample)
                        inferstate_ofstream << HSP_DELIM << std::setprecision(5)
                                            << std::fixed << param.smpl_values.find(p)->second[i];

                    // print augmented variable of id trace
                    if (!augm.smpl_values.empty())
                        for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                            for (const std::string &a : augm.list_draws)
                                inferstate_ofstream << HSP_DELIM << std::setprecision(5)
                                                    << std::fixed << augm.smpl_values[k].find(a)->second[i];

                    // print which chain
                    inferstate_ofstream << HSP_DELIM << (c + 1) << std::endl;
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
void MCMC::Impl::DataWide::print_open(const std::string &f)
{
    aux_ofstream_.open(f.c_str());

    if (aux_ofstream_.is_open())
    {
        // print default data headers
        aux_ofstream_ << "id" << HSP_DELIM << "house_id" << HSP_DELIM << "status";

        // print other data headers
        for (const std::string &d : headers_)
            aux_ofstream_ << HSP_DELIM << d;

        // print augmented variables headers
        for (const std::string &a : pImpl_.augm.list_export)
            aux_ofstream_ << HSP_DELIM << a;

        aux_ofstream_ << std::endl;
    }
    else
        throw std::runtime_error("Unable to open file: \"" + f + "\"");
}

/* Print all augmented variables of each individual */
void MCMC::Impl::DataWide::print_update(Household &house)
{
    if (aux_ofstream_.is_open())
        for (const Individual &member : house)
        {
            aux_ofstream_ << member.id()
                                  << HSP_DELIM << house.house_id()
                                  << HSP_DELIM << member.status();
            for (const std::string &d : headers_)
                aux_ofstream_ << HSP_DELIM << member.data(d);

            for (const std::string &a : pImpl_.augm.list_export)
                if (member.is_infected_status())
                    aux_ofstream_ << HSP_DELIM << tools::str::format(member.augm(a));
                else
                    aux_ofstream_ << HSP_DELIM << HSP_NA_DOUBLE;

            aux_ofstream_ << std::endl;
        }
    else
        throw std::runtime_error(
            "Output file closed unexpectedly.");
}

/* Print all augmented variables of each individual */
void MCMC::Impl::DataWide::print_close()
{
    if (aux_ofstream_.is_open())
        aux_ofstream_.close();
    else
        throw std::runtime_error(
            "Output file closed unexpectedly.");
}
#endif
