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

/* Import setup from R list and run build. */
void MCMC::Impl::import_setup(const Rcpp::List &input)
{
    tools::rlist_to_map(input, setup.deflist, "setup");
    setup.build();
}

/* Import infectivity profile from R list and run build. */
void MCMC::Impl::import_infec_profile(const Rcpp::List &input)
{
    tools::rlist_to_map(input, infec_profile.deflist, "infectivity profile");
    infec_profile.build();
}

/* Import data from R DataFrame and run build. */
void MCMC::Impl::import_dataset(const Rcpp::DataFrame &input)
{
    dict_str_vec datatable = {};
    tools::rdataframe_to_map(input, datatable, "dataset");
    build_dataset(datatable);
}

/* Import data from R DataFrame and run build. */
void MCMC::Impl::import_datalist(const Rcpp::DataFrame &input)
{
    dict_str_vec datatable = {};
    tools::rdataframe_to_map(input, datatable, "datalist");
    build_datalist(datatable);
}

/* Import multiplicity definitions from R DataFrame and run build. */
void MCMC::Impl::import_multiplicity(const Rcpp::DataFrame &input)
{
    dict_str_vec datatable = {};
    tools::rdataframe_to_map(input, datatable, "multiplicity");
    build_multiplicity(datatable);
}

/* Import params from R DataFrame and run build. */
void MCMC::Impl::import_param(const Rcpp::DataFrame &input)
{
    tools::rdataframe_to_map(input, _param_deftable, "parameter");
    build_param();
}

/* Import data from R DataFrame then run build. */
void MCMC::Impl::import_augm(const Rcpp::DataFrame &input)
{
    tools::rdataframe_to_map(input, _augm_deftable, "augmented");
    build_augm();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current setup configurations to R List. */
Rcpp::List MCMC::Impl::Setup::export_defs() const
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
                        _["id_trace"] = wrap(id_trace),
                        _["progress"] = progress);
}

/* Export current parameter configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_param() const
{
    using namespace Rcpp;

    try
    {
        std::vector<par_key> V_names = param_standard_infer_list + param_cache_infer_list + param_const_list;

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
            par_key p = V_names[r];
            V_value[r] = get_param(p);
            V_constant[r] = param_map.at(p).is_constant();
            V_discrete[r] = param_map.at(p).is_discrete();

            V_prior[r] = tools::stoRString(param_map.at(p).show_distr_name());
            V_arg_name[r] = tools::stoRString(param_map.at(p).show_distr_arg_names());
            V_arg_value[r] = tools::stoRString(param_map.at(p).show_distr_arg_values());

            V_rule[r] = tools::stoRString(param_map.at(p).show_prop_name());

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
        tools::error("While exporting parameter definitions:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC : Impl  : dataset exports  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_dataset()
{
    using namespace Rcpp;

    try
    {
        if (n_individuals == 0)
            return DataFrame();

        size_t n_rows = n_individuals;
        Rcpp::CharacterVector dataset_index = CharacterVector(n_rows);
        Rcpp::CharacterVector dataset_house_index = CharacterVector(n_rows);
        Rcpp::IntegerVector dataset_status = IntegerVector(n_rows);
        std::vector<Rcpp::NumericVector> dataset_dbl = {};

        // retrieve datalist headers
        str_vec headers = house_vector[0].get_member(0).keys_dataset();

        // other non standard dataset columns
        for (size_t s = 0; s < headers.size(); ++s)
            dataset_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : house_vector)
            for (const Individual &member : house)
            {
                dataset_index[r] = member.id();
                dataset_house_index[r] = house.house_id();
                dataset_status[r] = member.status();

                size_t col = 0;
                for (const data_key &d : headers)
                    dataset_dbl[col++][r] = member.data(d);
                ++r;
            }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(_["id"] = dataset_index,
                                             _["house_id"] = dataset_house_index,
                                             _["status"] = dataset_status);

        size_t col = 0;
        for (const data_key &d : headers)
            output.push_back(dataset_dbl[col++], d);

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting dataset:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC : Impl  : export datalist  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_datalist()
{
    using namespace Rcpp;

    try
    {
        if (indiv_datalist.size() == 0)
            return DataFrame();

        // retrieve datalist headers
        str_vec headers = house_vector[0].get_member(0).keys_datalist();

        // retrieve total number of rows
        size_t n_rows = 0;
        for (size_t i = 0; i < n_individuals; ++i)
            n_rows += indiv_datalist[i].at(headers[0]).size();

        Rcpp::CharacterVector datalist_index = CharacterVector(n_rows);
        std::vector<Rcpp::NumericVector> datalist_dbl = {};

        // other non standard datalist columns
        for (size_t s = 0; s < headers.size(); ++s)
            datalist_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : house_vector)
            for (const Individual &member : house)
                for (size_t k = 0; k < member.datalist(headers[0]).size(); ++k, ++r)
                {
                    datalist_index[r] = member.id();

                    size_t col = 0;
                    for (const data_key &d : headers)
                        datalist_dbl[col++][r] = member.datalist(d)[k];
                }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(_["id"] = datalist_index);

        size_t col = 0;
        for (const data_key &d : headers)
            output.push_back(datalist_dbl[col++], d);

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting datalist:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : report exports  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export current data configurations to R DataFrame. */
Rcpp::DataFrame MCMC::Impl::export_dataset_report() const
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
        tools::error("While exporting dataset report:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : export house size report  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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
        tools::error("While exporting household size report:\n" +
                     std::string(e.what()));
    }
    return DataFrame();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : export traces  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export MCMC traces of parameter posterior samples. */
Rcpp::NumericMatrix MCMC::Impl::export_traces() const
{
    using namespace Rcpp;

    try
    {
        // UPDATE: I shall allow for export even if _is_power_post
        if (smpl_full_loglike.empty())
            throw std::length_error("No traces to export were found.");

        //----- Collect the column names -----//
        str_vec headers = {"loglikelihood"};
        headers += param_standard_infer_list + param_cache_infer_list;

        // augmented variable of id trace
        if (!smpl_augm_value.empty())
        {
            if (smpl_individual_star.size() == 1)
                headers += augm_list_export();
            else
                for (size_t k = 0; k < smpl_individual_star.size(); ++k)
                    for (const augm_key &a : augm_list_export())
                        headers += {smpl_individual_star[k]->id() + "__" + a};
        }

        // multiple chains
        if (setup.n_chains > 1)
            headers += {"CHAIN__"};

        //----- NumericMatrix with traces -----//
        NumericMatrix output(setup.n_chains * setup.c_sample_size,
                             headers.size());
        size_t col = 0;

        // loglikelihood
        output(_, col++) = NumericVector(wrap(smpl_full_loglike));

        // parameters
        for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
            output(_, col++) = NumericVector(wrap(smpl_param_value.at(p)));

        // augmented variables
        if (!smpl_augm_value.empty())
            for (size_t k = 0; k < smpl_individual_star.size(); ++k)
            {
                const auto &smpl_augm_value_k = smpl_augm_value[k];
                for (const augm_key &a : augm_list_export())
                    output(_, col++) = NumericVector(wrap(smpl_augm_value_k.at(a)));
            }

        // retrieve chain identities
        if (setup.n_chains > 1)
        {
            IntegerVector chain_col(setup.c_sample_size * setup.n_chains);
            for (size_t c = 0; c < setup.n_chains; ++c)
                std::fill(chain_col.begin() + c * setup.c_sample_size,
                          chain_col.begin() + (c + 1) * setup.c_sample_size, c + 1);
            output(_, col++) = chain_col;
        }

        //----- associate column names -----//
        colnames(output) = CharacterVector(wrap(headers));

        //----- associate row names -----//
        IntegerVector sweep(setup.c_sample_size * setup.n_chains);
        for (size_t c = 0; c < setup.n_chains; ++c)
            for (size_t k = 0; k < setup.c_sample_size; k++)
                sweep[k + c * setup.c_sample_size] = setup.n_burnin + (k + 1) * setup.thin;
        rownames(output) = CharacterVector(wrap(sweep));

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting MCMC traces:\n" +
                     std::string(e.what()));
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  MCMC : Impl : export household loglikelihood matrix  //////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export all sampled values of the loglikelihood for all household. */
Rcpp::NumericMatrix MCMC::Impl::export_house_loglike_mtx() const
{
    using namespace Rcpp;

    try
    {
        if (smpl_full_loglike.empty())
            throw std::length_error("No traces to export were found.");

        // create a NumericMatrix with n_house columns
        ulong sample_size = setup.n_chains * setup.c_sample_size;
        NumericMatrix output(sample_size, n_houses());
        CharacterVector headers(n_houses());

        // retrieve restricted index by general index of multiplicity head class
        std::unordered_map<size_t, size_t> restr_idx_map = {};
        restr_idx_map.reserve(restrict_hh_indexes.size());
        for (size_t R = 0; R < restrict_hh_indexes.size(); ++R)
            restr_idx_map.emplace(restrict_hh_indexes[R], R);

        // assign each NumericVector to a column of the matrix
        for (size_t H = 0; H < n_houses(); ++H)
        {
            size_t R = restr_idx_map.at(house_vector[H].multiplicity_class_index());
            output(_, H) = NumericVector(wrap(smpl_house_loglike[R]));
            headers[H] = house_vector[H].house_id();
        }

        // associate column names
        colnames(output) = headers;

        // associate row names
        IntegerVector sweep(setup.c_sample_size * setup.n_chains);
        for (size_t c = 0; c < setup.n_chains; ++c)
            for (size_t k = 0; k < setup.c_sample_size; ++k)
                sweep[k + c * setup.c_sample_size] = setup.n_burnin + (k + 1) * setup.thin;
        rownames(output) = CharacterVector(wrap(sweep));

        return output;
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting household loglikelihood traces:\n" +
                     std::string(e.what()));
    }
    return NumericMatrix();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : infectivity profile export  //////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export a showcase of the selected infectivity profile. */
Rcpp::List MCMC::Impl::showcase_infec_profile(int n_bins) const
{
    using namespace Rcpp;

    try
    {
        // check the call makes sense
        // UPDATE: it would be nice by selecting id_trace...
        if (infec_profile.profile == nullptr)
            throw std::length_error("No infectivity profile to export.");
        else if (infec_profile.is_augmented())
            throw std::length_error(
                "Export of augmented infectivity profile is not supported.");
        else if (!infec_profile.profile->is_able_to_showcase())
            throw std::length_error(
                "Cannot export infectivity profile with undefined arguments.");

        // produce a showcase of the distribution profile
        std::vector<dbl_vec> values = infec_profile.profile->showcase(n_bins);

        // compute mean and variance
        const dbl_vec &x = values[0];
        const dbl_vec &prob = values[1];
        const size_t N = x.size();
        const double Dx = x[1] - x[0];

        double m = x[0] * prob[0] + x[N - 1] * prob[N - 1];
        double m2 = x[0] * x[0] * prob[0] + x[N - 1] * x[N - 1] * prob[N - 1];
        // integration wiht trapezoidal rule for continuos case
        if (!infec_profile.profile->is_discrete())
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

        return List::create(_["name"] = infec_profile.profile->show_name(),
                            _["discr"] = infec_profile.profile->is_discrete(),
                            _["mean"] = m,
                            _["var"] = m2 - m * m,
                            _["gentime"] = wrap(values[0]),
                            _["inst"] = wrap(values[1]),
                            _["cumul"] = wrap(values[2]));
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting infectivity profile:\n" +
                     std::string(e.what()));
    }
    return List();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMC : Impl : infectivity profile export ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Export a list containing all relevant information for a debug. */
Rcpp::List MCMC::Impl::export_for_debug(const std::string &house_id) const
{
    using namespace Rcpp;

    try
    {
        if (!_hhid_to_index.count(house_id))
            throw std::invalid_argument(
                "Index \"" + house_id + "\" does not match any household in the dataset.");

        const Household &house = house_vector[_hhid_to_index.at(house_id)];

        //--------------------------------------------//
        //----- Export dataset, augmented values -----//
        //--------------------------------------------//

        const size_t n_rows_dataset = house.size();
        const str_vec cols_dataset = house.get_member(0).keys_dataset();
        Rcpp::CharacterVector dataset_id = CharacterVector(n_rows_dataset);
        Rcpp::CharacterVector dataset_house_id = CharacterVector(n_rows_dataset);
        Rcpp::IntegerVector dataset_status = IntegerVector(n_rows_dataset);
        std::vector<Rcpp::NumericVector> dataset_dbl;
        dataset_dbl.reserve(cols_dataset.size());
        for (size_t i = 0; i < cols_dataset.size(); ++i)
            dataset_dbl.push_back(Rcpp::NumericVector(n_rows_dataset));

        const str_vec cols_augm = subtr(augm_list_full,
                                        {AUX_IP_INST_PROFILE, AUX_IP_CUMUL_PROFILE});
        std::vector<Rcpp::NumericVector> augm_dbl;
        augm_dbl.reserve(cols_augm.size());
        for (size_t i = 0; i < cols_augm.size(); ++i)
            augm_dbl.push_back(Rcpp::NumericVector(n_rows_dataset));

        size_t r = 0, c = 0;
        for (const Individual &member : house)
        {
            dataset_id[r] = member.id();
            dataset_house_id[r] = house.house_id();
            dataset_status[r] = member.status();

            c = 0;
            for (const data_key &d : cols_dataset)
                dataset_dbl[c++][r] = member.data(d);
            c = 0;
            for (const augm_key &a : cols_augm)
                augm_dbl[c++][r] = member.augm(a);
            ++r;
        }

        DataFrame output_dataset = DataFrame::create(_["id"] = dataset_id,
                                                     _["house_id"] = dataset_house_id,
                                                     _["status"] = dataset_status);
        c = 0;
        for (const data_key &d : cols_dataset)
            output_dataset.push_back(dataset_dbl[c++], d);

        DataFrame output_augm = DataFrame::create();
        if (cols_augm.size() > 0)
        {
            output_augm = DataFrame::create(_["id"] = dataset_id);
            c = 0;
            for (const augm_key &a : cols_augm)
                output_augm.push_back(augm_dbl[c++], a);
        }

        //---------------------------//
        //----- Export datalist -----//
        //---------------------------//

        DataFrame output_datalist;

        if (indiv_datalist.size() == 0)
            output_datalist = DataFrame::create();
        else
        {
            const str_vec cols_datalist = house.get_member(0).keys_dataset();
            size_t n_rows_datalist = 0;

            for (const Individual &member : house)
                n_rows_datalist += member.datalist(cols_datalist[0]).size();

            Rcpp::CharacterVector datalist_id = CharacterVector(n_rows_datalist);
            std::vector<Rcpp::NumericVector> datalist_dbl = {};
            for (size_t s = 0; s < cols_datalist.size(); ++s)
                datalist_dbl.emplace_back(n_rows_datalist);

            r = 0;
            for (const Individual &member : house)
                for (size_t k = 0; k < member.datalist(cols_datalist[0]).size(); ++k, ++r)
                {
                    datalist_id[r] = member.id();

                    c = 0;
                    for (const data_key &d : cols_datalist)
                        datalist_dbl[c++][r] = member.datalist(d)[k];
                }

            // assign vectors to dataframe
            output_datalist = DataFrame::create(_["id"] = datalist_id);

            size_t col = 0;
            for (const data_key &d : cols_datalist)
                output_datalist.push_back(datalist_dbl[col++], d);
        }

        //-----------------------------------//
        //----- Export parameter values -----//
        //-----------------------------------//

        NumericVector output_param;
        if (param_map.empty())
            output_param = NumericVector::create();
        else
        {
            output_param = NumericVector(param_map.size());
            CharacterVector names(param_map.size());

            c = 0;
            for (const auto &this_param : param_map)
            {
                names[c] = this_param.first;
                output_param[c++] = this_param.second.get();
            }
            output_param.names() = names;
        }

        return List::create(_["dataset"] = output_dataset,
                            _["datalist"] = output_datalist,
                            _["param"] = output_param,
                            _["augm"] = output_augm);
    }
    catch (const std::exception &e)
    {
        tools::error("While exporting house \"" + house_id + "\" for debug:\n" +
                     std::string(e.what()));
    }
    return List();
}

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
    tools::read_csv_list(f, setup.deflist, "setup");
    setup.build();
}

/* Read Setup file and produce a definition list. */
void MCMC::Impl::read_infec_profile(const std::string &f)
{
    tools::read_csv_list(f, infec_profile.deflist, "infectivity profile");
    infec_profile.build();
}

/* Read parameter definitions and run assignment. */
void MCMC::Impl::read_param(const std::string &f)
{
    tools::read_csv_table(f, _param_deftable, "parameter");
    build_param();
}

/* Read Model Augmented Variable file. */
void MCMC::Impl::read_augm(const std::string &f)
{
    tools::read_csv_table(f, _augm_deftable, "augmented definitions");
    build_augm();
}

/* Read dataset file and then run assignment. */
void MCMC::Impl::read_dataset(const std::string &f)
{
    dict_str_vec deftable = {};
    tools::read_csv_table(f, deftable, "dataset");
    build_dataset(deftable);
    deftable.clear();
}

/* Read datalist file and then run assignment. */
void MCMC::Impl::read_datalist(const std::string &f)
{
    dict_str_vec deftable = {};
    tools::read_csv_table(f, deftable, "datalist");
    build_datalist(deftable);
    deftable.clear();
}

/* Read multiplicity file and then run assignment. */
void MCMC::Impl::read_multiplicity(const std::string &f)
{
    dict_str_vec deftable = {};
    tools::read_csv_table(f, deftable, "multiplicity");
    build_multiplicity(deftable);
    deftable.clear();
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
const std::string _percent_(int n, int d) { return tools::percentage(static_cast<double>(n) / d, 1); }

/* Auxiliary function to call percentuage. */
const std::string _trail_percent_(int n, int d) { return std::to_string(n) + " (" + _percent_(n, d) + ")"; }

/* Print the report of the database and of the model group rates. */
void MCMC::Impl::show_summary() const
{
    std::ostringstream oss;

    oss << "n. of individuals: " << n_individuals
        << "\n    * n. of cases: "
        << _trail_percent_(get_report_indiv("infected"), n_individuals)
        << "\n    * n. of missing status : "
        << _trail_percent_(get_report_indiv("missing"), n_individuals);

    oss << "\nn. of households: " << n_houses()
        << "\n    * with no cases: "
        << _trail_percent_(get_report_house("no cases"), n_houses())
        << "\n    * with at least 1 case: "
        << _trail_percent_(get_report_house("at least 1 case"), n_houses())
        << "\n    * with at least 2 cases: "
        << _trail_percent_(get_report_house("at least 2 cases"), n_houses());

    oss << "\nmax n. house cases: " << get_report_hh_size("max cases").first
        << " (" << _percent_(get_report_hh_size("max cases").second, n_houses()) << ")"
        << "\nmin n. house members: " << get_report_hh_size("min size").first
        << " (" << _percent_(get_report_hh_size("min size").second, n_houses()) << ")"
        << "\nmax n. house members: " << get_report_hh_size("max size").first
        << " (" << _percent_(get_report_hh_size("max size").second, n_houses()) << ")";

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
std::string MCMC::Impl::Setup::show() const
{
    // common setup variables
    std::ostringstream oss;
    oss << "seed" << C_DELIM << seed
        << "\nn_threads" << C_DELIM << n_threads
        << "\nn_chains" << C_DELIM << n_chains
        << "\nsample_size" << C_DELIM << c_sample_size
        << "\nn_steps" << C_DELIM << n_steps
        << "\nn_adapt" << C_DELIM << n_adapt
        << "\nn_burnin" << C_DELIM << n_burnin
        << "\nthin" << C_DELIM << thin
        << "\nrand_init" << C_DELIM << tools::truth(rand_init_params)
        << "\nprogress" << C_DELIM << tools::truth(progress);
    return oss.str();
}

void MCMC::Impl::show_setup() const
{
    std::string id_trace_txt = HSP_NA_STRING;
    if (!smpl_individual_star.empty())
    {
        id_trace_txt = smpl_individual_star[0]->id();
        for (ushort k = 1; k < smpl_individual_star.size(); ++k)
            id_trace_txt += C_LIST_DELIM + smpl_individual_star[k]->id();
    }

    tools::print_stdout(setup.show() + "\nid_trace" + C_DELIM + id_trace_txt + "\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : show parameters  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Show the user model parameters */
void MCMC::Impl::show_all_params() const
{
    std::ostringstream oss;

    // show the headers
    oss << param_file_headers[0];
    for (size_t k = 1; k < param_file_headers.size(); ++k)
        oss << C_DELIM << param_file_headers[k];
    oss << "\n";

    // show parameter definitions
    std::vector<par_key> all_params = param_standard_infer_list + param_cache_infer_list + param_const_list;
    for (const par_key &p : all_params)
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

    for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
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
            for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
                inferstate_ofstream << C_DELIM << p;

            // augmented variable of id trace
            if (!smpl_augm_value.empty())
            {
                if (smpl_individual_star.size() == 1)
                    for (const augm_key &a : augm_list_export())
                        inferstate_ofstream << C_DELIM << a;
                else
                    for (size_t k = 0; k < smpl_individual_star.size(); ++k)
                        for (const augm_key &a : augm_list_export())
                            inferstate_ofstream << smpl_individual_star[k]->id() + "__" + a;
            }

            // header chain
            if (setup.n_chains > 1)
                inferstate_ofstream << C_DELIM << "CHAIN__";

            inferstate_ofstream << std::endl;

            //----- VALUES -----//

            for (size_t c = 0; c < setup.n_chains; ++c)
                for (size_t k = 0; k < setup.c_sample_size; k++)
                {
                    size_t i = k + c * setup.c_sample_size;

                    // print iteration
                    ulong sweep = setup.n_burnin + (k + 1) * setup.thin;
                    inferstate_ofstream << sweep;

                    // print loglikelihood
                    inferstate_ofstream << C_DELIM << std::setprecision(2)
                                        << std::fixed << smpl_full_loglike[i];

                    // print parameter values
                    for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
                        inferstate_ofstream << C_DELIM << std::setprecision(5)
                                            << std::fixed << smpl_param_value.at(p)[i];

                    // print augmented variable of id trace
                    if (!smpl_augm_value.empty())
                        for (size_t k = 0; k < smpl_individual_star.size(); ++k)
                            for (const augm_key &a : augm_list_export())
                                inferstate_ofstream << C_DELIM << std::setprecision(5)
                                                    << std::fixed << smpl_augm_value[k].at(a)[i];

                    // print which chain
                    if (setup.n_chains > 1)
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
    aux_dataset_ofstream.open(f.c_str());
    if (aux_dataset_ofstream.is_open())
    {
        // print default data headers
        aux_dataset_ofstream << "id" << C_DELIM << "house_id" << C_DELIM << "status";

        // print other data headers
        for (const data_key &d : house_vector[0].get_member(0).keys_dataset())
            aux_dataset_ofstream << C_DELIM << d;

        // print augmented variables headers
        for (const augm_key &a : augm_list_export())
            aux_dataset_ofstream << C_DELIM << a;

        aux_dataset_ofstream << std::endl;
    }
    else
        throw std::runtime_error(
            "Unable to open file: " + f);
}

/* Print all augmented variables of each individual */
void MCMC::Impl::_aux_print_dataset_update_(Household &house)
{
    if (aux_dataset_ofstream.is_open())
        for (const Individual &member : house)
        {
            aux_dataset_ofstream << member.id()
                                 << C_DELIM << house.house_id()
                                 << C_DELIM << member.status();
            for (const data_key &d : house_vector[0].get_member(0).keys_dataset())
                aux_dataset_ofstream << C_DELIM << member.data(d);

            for (const augm_key &a : augm_list_export())
                if (member.status() == INFECTED_STATUS)
                    aux_dataset_ofstream << tools::na_to_str(member.augm(a));
                else
                    aux_dataset_ofstream << C_DELIM << HSP_NA_DOUBLE;

            aux_dataset_ofstream << std::endl;
        }
    else
        throw std::runtime_error(
            "Output file closed unexpectedly.");
}

/* Print all augmented variables of each individual */
void MCMC::Impl::_aux_print_dataset_finalize_()
{
    if (aux_dataset_ofstream.is_open())
        aux_dataset_ofstream.close();
    else
        throw std::runtime_error(
            "Output file closed unexpectedly.");
}

#endif
