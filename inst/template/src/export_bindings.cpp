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
 * export_bindings.cpp
 *
 */

#include <Rcpp.h>

#include "tools/flags.h"
#include "tools/na_utils.h"
#include "tools/error.h"
#include "tools/read.h"
#include "tools/params.h"
#include "function.h"
#include "function_registry.h"
#include "distribution_registry.h"
#include "mcmc_impl.h"
#include "model.h"

RCPP_MODULE(housepi_flags)
{
    Rcpp::class_<CompiledFlagsHousepi>("CompiledFlagsHousepi")
        .constructor<>()
        .method("get_flags", &CompiledFlagsHousepi::get_flags);
}

// macro to define Rcpp wrapper classes to export to R
#define MAKE_HOUSEPI_RCPP_WRAPPER(Wrapper, family, set_eval, pdf_call, cdf_call, ...)          \
    class Wrapper##Housepi                                                                     \
    {                                                                                          \
    public:                                                                                    \
        Wrapper##Housepi(std::string name_, SEXP args_, double min, double max, bool is_discr) \
        {                                                                                      \
            try                                                                                \
            {                                                                                  \
                ptr_ = compose_func_registries(name_, family, {__VA_ARGS__}, false);           \
                if (args_ != R_NilValue)                                                       \
                {                                                                              \
                    Rcpp::List args(args_);                                                    \
                    const Rcpp::CharacterVector arg_names = args.names();                      \
                    for (R_xlen_t i = 0; i < args.size(); ++i)                                 \
                    {                                                                          \
                        const std::string this_arg = Rcpp::as<std::string>(arg_names[i]);      \
                        ptr_->init_attr(this_arg, Attribute(Rcpp::as<double>(args[i])));       \
                    }                                                                          \
                }                                                                              \
                                                                                               \
                if (!hsp::na::is(min))                                                         \
                    ptr_->init_attr(hsp::params::MIN, Attribute(min));                         \
                if (!hsp::na::is(max))                                                         \
                    ptr_->init_attr(hsp::params::MAX, Attribute(max));                         \
                                                                                               \
                if (is_discr)                                                                  \
                    ptr_->init_discretization();                                               \
                ptr_->init_truncation();                                                       \
            }                                                                                  \
            catch (const std::exception &e)                                                    \
            {                                                                                  \
                hsp::msg::error(                                                               \
                    "While initializing " family " '" + name_ + "':\n" + e.what());            \
            }                                                                                  \
        }                                                                                      \
        ~Wrapper##Housepi() = default;                                                         \
        Rcpp::NumericVector pdf_call(Rcpp::NumericVector x) const                              \
        {                                                                                      \
            ptr_->set_eval##_call();                                                           \
            return call_(x);                                                                   \
        }                                                                                      \
        Rcpp::NumericVector cdf_call(Rcpp::NumericVector x) const                              \
        {                                                                                      \
            ptr_->set_eval##_primitive();                                                      \
            return call_(x);                                                                   \
        }                                                                                      \
                                                                                               \
    private:                                                                                   \
        std::unique_ptr<Function> ptr_;                                                        \
                                                                                               \
        Rcpp::NumericVector call_(Rcpp::NumericVector x) const                                 \
        {                                                                                      \
            try                                                                                \
            {                                                                                  \
                const R_xlen_t n = x.size();                                                   \
                Rcpp::NumericVector y(n);                                                      \
                for (R_xlen_t i = 0; i < n; ++i)                                               \
                    y[i] = ptr_->evaluate(x[i]);                                               \
                return y;                                                                      \
            }                                                                                  \
            catch (const std::exception &e)                                                    \
            {                                                                                  \
                hsp::msg::error(                                                               \
                    "While calling '" family "'" + ptr_->name() + "':\n" + e.what());          \
            }                                                                                  \
        }                                                                                      \
    };

//' @name FunctionHousepi
//' @title Function from Housepi
//' @description Call housepi registered function.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)` with function name.
//'   \item Parameter: args A named `numeric()` with arguments.
//'   \item Parameter: min `numeric(1)`\cr Minimum value of the function.
//'   \item Parameter: max `numeric(1)`\cr Maximum value of the function.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discretized.
//' }
//' @field func\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
//' @field prim\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_function)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(Function, "function",
                              set_eval, func, prim,
                              &UserFuncRegistry::inst(),
                              &TranscendentalRegistry::inst(),
                              &ArithmeticRegistry::inst());

    Rcpp::class_<FunctionHousepi>("FunctionHousepi")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("func", &FunctionHousepi::func)
        .method("prim", &FunctionHousepi::prim);
}

//' @name DistributionHousepi
//' @title Distribution from Housepi
//' @description Call housepi registered distribution.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the distribution.
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: min `numeric(1)`\cr Minimum value of the distribution.
//'   \item Parameter: max `numeric(1)`\cr Maximum value of the distribution.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discretized.
//' }
//' @field pdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
//' @field cdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_distribution)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(Distribution, "distribution",
                              set_eval_as_distr, pdf, cdf,
                              &ContinuosDistrRegistry::inst(),
                              &DiscreteDistrRegistry::inst());

    Rcpp::class_<DistributionHousepi>("DistributionHousepi")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("pdf", &DistributionHousepi::pdf)
        .method("cdf", &DistributionHousepi::cdf);
}

//' @name InfecProfileHousepi
//' @title Infectivity Profile from Housepi
//' @description Call housepi registered infectivity profile.
//' @field initialize\cr \itemize{
//'   \item Parameter: name `character(1)`\cr Name of the profile.
//'   \item Parameter: args `numeric()`\cr Custom arguments.
//'   \item Parameter: min `numeric(1)`\cr Minimum value of the profile.
//'   \item Parameter: max `numeric(1)`\cr Maximum value of the profile.
//'   \item Parameter: is_discr `logical(1)`\cr Wheter discretized.
//' }
//' @field pdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
//' @field cdf\cr \itemize{
//'   \item Parameter: x `numeric()`\cr Input vector.
//' }
//' @return `numeric()`
RCPP_MODULE(housepi_infec_profile)
{
    MAKE_HOUSEPI_RCPP_WRAPPER(InfecProfile, "infectivity profile",
                              set_eval_as_distr, pdf, cdf,
                              &InfecProfileRegistry::inst(),
                              &ContinuosDistrRegistry::inst(),
                              &DiscreteDistrRegistry::inst());

    Rcpp::class_<InfecProfileHousepi>("InfecProfileHousepi")
        .constructor<std::string, SEXP, double, double, bool>()
        .method("pdf", &InfecProfileHousepi::pdf)
        .method("cdf", &InfecProfileHousepi::cdf);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : imports from Rcpp /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Auxiliary converter to R NA real values.
inline double hsp_to_rcpp(double x) { return hsp::na::is(x) ? NA_REAL : x; }

/// Convert R list.
u_dict<std::string> rlist_to_strmap(SEXP input)
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

/// Convert R dataframe map of string vectors.
u_dict<str_vec> rdataframe_to_strmap(SEXP input)
{
    using namespace Rcpp;

    if (Rf_isNull(input))
        return {};

    DataFrame df(input);
    CharacterVector colnames = df.names();
    u_dict<str_vec> output;

    for (R_xlen_t i = 0; i < df.size(); ++i)
    {
        CharacterVector col = df[i];
        str_vec col_data(col.begin(), col.end());
        output[as<std::string>(colnames[i])] = std::move(col_data);
    }

    return output;
}

/// Convert R dataframe map of double vectors.
u_dict<std::vector<double>> rdataframe_to_dblmap(SEXP input)
{
    using namespace Rcpp;

    if (Rf_isNull(input))
        return {};

    DataFrame df(input);
    CharacterVector colnames = df.names();
    u_dict<std::vector<double>> output;

    for (R_xlen_t i = 0; i < df.size(); ++i)
    {
        NumericVector col = df[i];
        std::vector<double> col_data(col.begin(), col.end());
        output[as<std::string>(colnames[i])] = std::move(col_data);
    }

    return output;
}

#define MAKE_RCPP_IMPORT_METHOD(module, section_name, r_to_map) \
    void MCMC::import_##module(SEXP input)                      \
    {                                                           \
        try                                                     \
        {                                                       \
            pImpl_->module.build(r_to_map(input));              \
        }                                                       \
        catch (const std::exception &e)                         \
        {                                                       \
            hsp::msg::error(                                    \
                std::string("While importing ") +               \
                section_name + ":\n" + e.what());               \
        }                                                       \
    }

MAKE_RCPP_IMPORT_METHOD(setup, "setup", rlist_to_strmap)
MAKE_RCPP_IMPORT_METHOD(infec_profile, "infectivity profile", rlist_to_strmap)
MAKE_RCPP_IMPORT_METHOD(param, "parameters", rdataframe_to_strmap)
MAKE_RCPP_IMPORT_METHOD(data_wide, "data-wide", rdataframe_to_strmap)
MAKE_RCPP_IMPORT_METHOD(multiplicity, "multiplicity", rdataframe_to_strmap)
MAKE_RCPP_IMPORT_METHOD(data_long, "data-long", rdataframe_to_strmap)
MAKE_RCPP_IMPORT_METHOD(group_rates, "group rates", rdataframe_to_strmap)
MAKE_RCPP_IMPORT_METHOD(augm, "augmented", rdataframe_to_strmap)

void MCMC::import_param_for_expec(SEXP input)
{
    try
    {
        pImpl_->param.build_for_expec(rdataframe_to_dblmap(input));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(
            "While importanting parameter table:\n" +
            std::string(e.what()));
    }
}

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
Rcpp::List MCMC::export_setup() const
{
    using namespace Rcpp;
    const Impl::Setup &setup = pImpl_->setup;

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
    const Impl::InfecProfile &infec_profile = pImpl_->infec_profile;

    if (infec_profile.is_augmented())
        hsp::msg::error(
            "Impossible to export infectivity profile when if definition is not generic.");

    if (infec_profile.generic_profile->is_empty())
        return List::create();

    return List::create(_["name"] = infec_profile.generic_profile->name(),
                        // NOTE: arguments should not be polymorphic when called here
                        _["args"] = rcpp_export_args(*infec_profile.generic_profile),
                        _["min"] = infec_profile.generic_profile->lower(),
                        _["max"] = infec_profile.generic_profile->upper(),
                        _["discr"] = infec_profile.generic_profile->is_discrete());
}

/// Export to R all augmented variables at the current state.
Rcpp::DataFrame MCMC::export_augmented_data_wide() const
{
    using namespace Rcpp;
    const Impl::Augmented &augm = pImpl_->augm;

    const size_t n_cols = pImpl_->data_wide.report_indiv.total;

    int_vec status_col(n_cols);
    std::vector<dbl_vec> augm_cols(augm.list_export.size());
    for (size_t i = 0; i < augm.list_export.size(); ++i)
        augm_cols[i] = dbl_vec(n_cols);

    size_t row = 0;
    for (Household &house : pImpl_->data_wide)
        for (const Individual &member : house)
        {
            status_col[row] = member.status();
            for (size_t i = 0; i < augm.list_export.size(); ++i)
                augm_cols[i][row] = hsp_to_rcpp(member.augm(augm.list_export[i]));
            ++row;
        }

    // create the ouput List
    List output = List::create(_["status"] = IntegerVector(wrap(std::move(status_col))));
    for (size_t i = 0; i < augm.list_export.size(); ++i)
        output[augm.list_export[i]] = NumericVector(wrap(std::move(augm_cols[i])));
    return DataFrame(output);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : DataWide : exports  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Rcpp::DataFrame MCMC::helper_export_data_wide_(const str_vec &columns, bool is_status) const
{
    using namespace Rcpp;
    const Impl::DataWide &data_wide = pImpl_->data_wide;

    try
    {
        if (data_wide.is_empty())
            return DataFrame();

        const size_t n_rows = data_wide.report_indiv.total;
        CharacterVector df_id = CharacterVector(n_rows);
        CharacterVector df_house_id = CharacterVector(n_rows);
        IntegerVector df_status = IntegerVector(n_rows);

        // other non standard data columns
        std::vector<NumericVector> df_double;
        df_double.reserve(columns.size());
        for (size_t c = 0; c < columns.size(); ++c)
            df_double.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : pImpl_->data_wide)
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
    return helper_export_data_wide_(pImpl_->data_wide.headers, true);
}
Rcpp::DataFrame MCMC::data_wide_for_simulation() const
{
    return helper_export_data_wide_(pImpl_->data_wide.headers - pImpl_->augm.list_export, false);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : DataLong : exports  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::export_data_long() const
{
    using namespace Rcpp;
    const Impl::DataLong &data_long = pImpl_->data_long;

    try
    {
        if (data_long.indiv_datalist.size() == 0)
            return DataFrame();

        // retrieve total number of rows
        size_t n_rows = 0;
        const str_vec &headers = data_long.headers;
        for (size_t i = 0; i < data_long.indiv_datalist.size(); ++i)
            n_rows += data_long.indiv_datalist[i].at(headers[0]).size();

        CharacterVector datalist_index = CharacterVector(n_rows);

        // other non standard data-long columns
        std::vector<NumericVector> datalist_dbl = {};
        datalist_dbl.reserve(headers.size());
        for (size_t c = 0; c < headers.size(); ++c)
            datalist_dbl.emplace_back(n_rows);

        size_t r = 0;
        for (const Household &house : pImpl_->data_wide)
            for (const Individual &member : house)
                for (size_t k = 0; k < member.data_long(headers[0]).size(); ++k, ++r)
                {
                    datalist_index[r] = member.id();
                    for (size_t c = 0; c < headers.size(); ++c)
                        datalist_dbl[c][r] = hsp_to_rcpp(member.data_long(headers[c])[k]);
                }

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(_["id"] = datalist_index);
        for (size_t c = 0; c < headers.size(); ++c)
            output.push_back(datalist_dbl[c], headers[c]);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting data (long):\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : report exports  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

inline double format_percentage(double value, double total) { return std::round(10000. * value / total) / 100.; }

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::summarise_cases() const
{
    using namespace Rcpp;
    const Impl::DataWide &data_wide = pImpl_->data_wide;
    const ReportIndividual &report_indiv = data_wide.report_indiv;
    const ReportHousehold &report_house = data_wide.report_house;

    try
    {
        if (data_wide.is_empty())
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

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : export house size report  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export current data configurations to R DataFrame.
Rcpp::DataFrame MCMC::summarise_house_sizes() const
{
    using namespace Rcpp;
    const Impl::DataWide &data_wide = pImpl_->data_wide;
    const ReportIndividual &report_indiv = data_wide.report_indiv;
    const ReportHouseSizes &report_hh_size = data_wide.report_hh_size;

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

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : export draws  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Chain tag of each trace entry.
Rcpp::IntegerVector draws_chain_tags(unsigned long c_sample_size, unsigned long n_chains)
{
    using namespace Rcpp;

    IntegerVector output(c_sample_size * n_chains);
    for (size_t c = 0; c < n_chains; ++c)
        std::fill(output.begin() + c * c_sample_size,
                  output.begin() + (c + 1) * c_sample_size, c + 1);
    return output;
}

/// Iteration name of each trace entry.
Rcpp::CharacterVector draws_iter_names(unsigned long c_sample_size, unsigned long n_chains,
                                       unsigned long n_burnin, unsigned long n_thin)
{
    using namespace Rcpp;

    std::vector<unsigned long> iteration(c_sample_size * n_chains);
    for (size_t c = 0; c < n_chains; ++c)
        for (size_t k = 0; k < c_sample_size; k++)
            iteration[k + c * c_sample_size] = n_burnin + (k + 1) * n_thin;
    return CharacterVector(wrap(std::move(iteration)));
}

/// Export MCMC draws of parameter posterior samples.
Rcpp::NumericMatrix MCMC::export_draws() const
{
    using namespace Rcpp;
    const Impl::Setup &setup = pImpl_->setup;
    const Impl::Parameter &param = pImpl_->param;
    const Impl::Augmented &augm = pImpl_->augm;

    try
    {
        // UPDATE: I shall allow for export even if _is_power_post
        if (pImpl_->smpl_total_log_lik.empty())
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
                for (const std::string &a : pImpl_->augm.list_draws)
                    headers.push_back(augm.smpl_individuals[k]->id() + "_" + a);
        // multiple chains
        headers.push_back(hsp::read::CHAIN_COLUMN);

        // matrix with draws
        NumericMatrix output(setup.c_sample_size * setup.n_chains, headers.size());

        size_t col = 0;
        // log-likelihood
        output(_, col++) = NumericVector(wrap(pImpl_->smpl_total_log_lik));
        // parameter values
        for (const std::string &p : param.list_sample)
            output(_, col++) = NumericVector(wrap(pImpl_->param.smpl_values.at(p)));
        // parameter values
        for (const std::string &p : param.list_sample)
            output(_, col++) = NumericVector(wrap(pImpl_->param.smpl_priors.at(p)));

        // augmented variables
        if (!augm.smpl_values.empty())
            for (size_t k = 0; k < augm.smpl_values.size(); ++k)
                for (const std::string &a : pImpl_->augm.list_draws)
                    output(_, col++) = NumericVector(wrap(augm.smpl_values.at(k).at(a)));

        // asssociate chain tag
        output(_, col++) = draws_chain_tags(setup.c_sample_size, setup.n_chains);
        // associate column names
        colnames(output) = CharacterVector(wrap(headers));
        // associate row names
        rownames(output) = draws_iter_names(setup.c_sample_size, setup.n_chains,
                                            setup.n_burnin, setup.n_thin);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting MCMC draws:\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  MCMC : Impl : export household log-likelihood matrix  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export all sampled values of the log-likelihood for all household.
Rcpp::NumericMatrix MCMC::export_house_log_lik_draws() const
{
    using namespace Rcpp;
    const Impl::Setup &setup = pImpl_->setup;
    const Impl::DataWide &data_wide = pImpl_->data_wide;
    const Impl::Multiplicity &multiplicity = pImpl_->multiplicity;

    try
    {
        if (data_wide.smpl_values.empty())
            throw std::length_error("No draws to export were found.");

        // create a NumericMatrix with n_house columns
        size_t n_houses = data_wide.report_house.total;
        const unsigned long sample_size = setup.n_chains * setup.c_sample_size;
        NumericMatrix output(sample_size, n_houses + 1);
        CharacterVector columns(n_houses + 1);
        columns[n_houses] = hsp::read::CHAIN_COLUMN;

        // retrieve restricted index by general index of multiplicity head class
        std::unordered_map<size_t, size_t> restr_idx_map = {};
        restr_idx_map.reserve(multiplicity.indexes.size());
        for (size_t R = 0; R < multiplicity.indexes.size(); ++R)
            restr_idx_map.emplace(multiplicity.indexes[R], R);

        // assign each NumericVector to a column of the matrix
        for (const Household &house : data_wide)
        {
            size_t R = restr_idx_map.at(house.multiplicity_class_index());
            size_t H = house.house_index();
            output(_, H) = NumericVector(wrap(data_wide.smpl_values[R]));
            columns[H] = house.house_id();
        }

        // asssociate chain tag
        output(_, n_houses) = draws_chain_tags(setup.c_sample_size, setup.n_chains);
        // associate column names
        colnames(output) = columns;
        // associate row names
        rownames(output) = draws_iter_names(setup.c_sample_size, setup.n_chains,
                                            setup.n_burnin, setup.n_thin);

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting household log-likelihood draws:\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : house state snapshot //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Export to R a snapshot of current state of house data.
Rcpp::NumericVector MCMC::helper_param_rcpp_snapshot_() const
{
    using namespace Rcpp;
    const Impl::Parameter &param = pImpl_->param;
    const str_vec list_full = param.list_sample + param.list_const;

    if (list_full.empty())
        return NumericVector::create();

    NumericVector output = NumericVector(list_full.size());
    size_t c = 0;
    for (const std::string &p : list_full)
        output[c++] = param.get(p);
    output.names() = CharacterVector(wrap(list_full));

    return output;
}
/// Export to R a snapshot of current state of house data.
Rcpp::DataFrame MCMC::helper_data_wide_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const Impl::DataWide &data_wide = pImpl_->data_wide;
    const str_vec &headers = data_wide.headers;

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);
    CharacterVector df_house_id = CharacterVector(n_rows);
    IntegerVector df_status = IntegerVector(n_rows);

    std::vector<NumericVector> df_double;
    df_double.reserve(headers.size());
    for (size_t c = 0; c < headers.size(); ++c)
        df_double.emplace_back(NumericVector(n_rows));

    for (const Individual &member : house)
    {
        df_id[member] = member.id();
        df_house_id[member] = house.house_id();
        df_status[member] = member.status();

        for (size_t c = 0; c < headers.size(); ++c)
            df_double[c][member] = hsp_to_rcpp(member.data(headers[c]));
    }

    DataFrame output = DataFrame::create(_["id"] = std::move(df_id),
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
    const Impl::DataLong &data_long = pImpl_->data_long;
    const str_vec &headers = data_long.headers;

    if (data_long.indiv_datalist.empty())
        return DataFrame::create();

    size_t n_rows = 0;

    for (const Individual &member : house)
        n_rows += member.data_long(headers[0]).size();

    CharacterVector datalist_id = CharacterVector(n_rows);

    std::vector<NumericVector> datalist_dbl = {};
    datalist_dbl.reserve(headers.size());
    for (size_t i = 0; i < headers.size(); ++i)
        datalist_dbl.emplace_back(n_rows);

    size_t r = 0;
    for (const Individual &member : house)
        for (size_t i = 0; i < member.data_long(headers[0]).size(); ++i)
        {
            datalist_id[r] = member.id();
            for (size_t c = 0; c < headers.size(); ++c)
                datalist_dbl[c][r] = hsp_to_rcpp(member.data_long(headers[c])[i]);
            ++r;
        }

    // assign vectors to dataframe
    DataFrame output = DataFrame::create(_["id"] = std::move(datalist_id));
    for (size_t c = 0; c < headers.size(); ++c)
        output.push_back(datalist_dbl[c], headers[c]);

    return output;
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::helper_group_rates_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;

    const str_vec &headers = pImpl_->group_rates.headers;

    if (headers.empty())
        return DataFrame::create();

    const size_t n_rows = house.size();
    CharacterVector df_id = CharacterVector(n_rows);

    std::vector<NumericVector> group_rate_dbl;
    group_rate_dbl.reserve(headers.size());
    for (size_t c = 0; c < headers.size(); ++c)
        group_rate_dbl.emplace_back(NumericVector(n_rows));

    for (const Individual &member : house)
    {
        df_id[member] = member.id();

        for (size_t c = 0; c < headers.size(); ++c)
            group_rate_dbl[c][member] = hsp_to_rcpp(member.group_rate(headers[c]));
    }

    DataFrame output = DataFrame::create(_["id"] = df_id);
    for (size_t c = 0; c < headers.size(); ++c)
        output.push_back(group_rate_dbl[c], headers[c]);

    return output;
}

/// Export to R a snapshot of current state of house infectivity profile.
Rcpp::List MCMC::helper_infec_profile_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const Impl::InfecProfile &infec_profile = pImpl_->infec_profile;

    if (!infec_profile.is_augmented() && infec_profile.generic_profile->is_empty())
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
        if (infec_profile.is_augmented())
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
            Rcpp::List generic = export_infec_profile();
            infprof_name[member] = wrap(generic["name"]);
            infprof_discr[member] = generic["discr"];
            infprof_args[member] = generic["args"];
            infprof_min[member] = generic["min"];
            infprof_max[member] = generic["max"];
        }

    return List::create(_["id"] = df_id,
                        _["discr"] = infprof_discr,
                        _["name"] = infprof_name,
                        _["args"] = infprof_args,
                        _["min"] = infprof_min,
                        _["max"] = infprof_max);
}

/// Export to R a snapshot of current state of house augmented variables.
Rcpp::DataFrame MCMC::helper_augm_rcpp_snapshot_(const Household &house) const
{
    using namespace Rcpp;
    const Impl::Augmented &augm = pImpl_->augm;
    const str_vec &list_export = augm.list_export;

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

/// Export a list containing all relevant information for a debug.
Rcpp::List MCMC::export_house_snapshot(const std::string &house_id) const
{
    using namespace Rcpp;

    try
    {
        if (pImpl_->data_wide.is_empty())
            return List::create();

        const Household &house = pImpl_->data_wide(house_id);

        return List::create(_["param"] = helper_param_rcpp_snapshot_(),
                            _["data_wide"] = helper_data_wide_rcpp_snapshot_(house),
                            _["data_long"] = helper_data_long_rcpp_snapshot_(house),
                            _["group_rates"] = helper_group_rates_rcpp_snapshot_(house),
                            _["infec_profile"] = helper_infec_profile_rcpp_snapshot_(house),
                            _["augm"] = helper_augm_rcpp_snapshot_(house));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While exporting snapshot of household '" + house_id + "':\n" +
                        std::string(e.what()));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : Impl : simulate ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Perform simulation and return simulated data.
Rcpp::List MCMC::helper_simulate_(size_t n, bool is_expec_measure)
{
    using namespace Rcpp;
    Impl::Parameter &param = pImpl_->param;

    pImpl_->atomic_error.reset();

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
            if (pImpl_->atomic_error.is_pending())
                continue;
            output[r] = helper_simulate_single_repetition_();
        }
        pImpl_->atomic_error.throw_exception();

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
        const str_vec &augm_list_export = pImpl_->augm.list_export;
        const size_t n_rows = pImpl_->data_wide.report_indiv.total;

        IntegerVector df_status = IntegerVector(n_rows);

        std::vector<NumericVector> df_augm;
        df_augm.reserve(augm_list_export.size());
        for (size_t c = 0; c < augm_list_export.size(); ++c)
            df_augm.emplace_back(n_rows);

        size_t r = 0;
        for (Household house : pImpl_->data_wide) // NOTE: deep copy of households
            try
            {
                pImpl_->simulate_single_household_in_place(house);
                for (const Individual &member : house)
                {
                    df_status[r] = member.status();
                    for (size_t c = 0; c < augm_list_export.size(); ++c)
                        df_augm[c][r] = hsp_to_rcpp(member.augm(augm_list_export[c]));
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
        for (size_t c = 0; c < augm_list_export.size(); ++c)
            output.push_back(df_augm[c], augm_list_export[c]);
    }
    catch (const std::exception &e)
    {
        pImpl_->atomic_error.set_if_first(e);
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
Rcpp::List MCMC::log_evidence(const std::string &scheme, unsigned long n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = pImpl_->power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < pImpl_->setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(_["estimate"] = output,
                            _["scheme"] = scheme,
                            _["n_rungs"] = pImpl_->setup.n_rungs,
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

//' @name ModelHousepi
//' @title Internal housepi model
//' @description Perform inference and simulation from the model.
//' @field infer\cr Perform MCMC inference of the model parameters. \itemize{
//'   \item Parameter: output_file - The file path where the output will be saved
//'   \item Parameter: progbar - Display the progress bar
//' }
//' @field data_wide_for_simulation\cr Dataset without status and augmented entries.
//' @field simulate\cr Simulated status and augmented variables.
//' @field simulate_for_expec\cr Simulation for the computation of expected measure.
//' @field evaluate\cr Evaluate the log-likelihood according to the current state.
//' @field log_evidence\cr Evaluate the logarithm of the marginal likelihood (aka evidence).
//' Uses the method of power posteriors. \itemize{
//'   \item Parameter: scheme
//'   \item Parameter: n_rungs
//' }
//' @field import_setup\cr Import setup from R List.
//' @field import_data_wide\cr Import house indexed data in wide format from R Data.Frame.
//' @field import_data_long\cr Import individual indexed in long format from R Data.Frame.
//' @field import_multiplicity\cr Imoprt multiplicity definitions from R Data.Frame.
//' @field import_group_rates\cr Imoprt group rates definitions from R Data.Frame.
//' @field import_param\cr Import parameter definitions from R Data.Frame.
//' @field import_infec_profile\cr Import infectivity profile definition from R List.
//' @field import_augm\cr Import augmented variable definitions from R Data.Frame.
//' @field import_param_for_expec\cr Import parameter tables for expected measures estimation.
//' @field build_all_rates\cr Assign all contact and group rates.
//' @field export_setup\cr Return a List with current setup.
//' @field export_param_list_const\cr Return a CharacterVector with constant parameters.
//' @field export_param_list_sample_standard\cr Return a CharacterVector with parameters to sample.
//' @field export_param_list_sample_cache\cr Return a CharacterVector with parameters to sample with cached dependecies.
//' @field export_infec_profile\cr Return a List with infectivity profile.
//' @field export_data_wide\cr Return a DataFrame with current database.
//' @field summarise_cases\cr Return a DataFrame with a report of the current database.
//' @field summarise_house_sizes\cr Return a DataFrame with a report of household sizes.
//' @field export_draws\cr Return a NumericMatrix with MCMC draws.
//' @field export_house_log_lik_draws\cr Return a NumericMatrix with household log-likelihood draws.
//' @field export_house_snapshot\cr Return a List with a snapshot containing full household and model details.

RCPP_MODULE(housepi_model_core)
{
    Rcpp::class_<MCMC>("McmcHousepi")
        .constructor()
        .method("init_advanced_setup", &MCMC::init_advanced_setup)
        .method("clear_memory", &MCMC::clear_memory)
        .method("infer", &MCMC::infer)
        .method("data_wide_for_simulation", &MCMC::data_wide_for_simulation)
        .method("simulate", &MCMC::simulate)
        .method("simulate_for_expec", &MCMC::simulate_for_expec)
        .method("evaluate_log_lik", &MCMC::evaluate_log_lik)
        .method("log_evidence", &MCMC::log_evidence)
        .method("import_data_wide", &MCMC::import_data_wide)
        .method("import_data_long", &MCMC::import_data_long)
        .method("import_multiplicity", &MCMC::import_multiplicity)
        .method("import_group_rates", &MCMC::import_group_rates)
        .method("import_setup", &MCMC::import_setup)
        .method("import_param", &MCMC::import_param)
        .method("import_infec_profile", &MCMC::import_infec_profile)
        .method("import_augm", &MCMC::import_augm)
        .method("import_param_for_expec", &MCMC::import_param_for_expec)
        .method("build_all_rates", &MCMC::build_all_rates)
        .method("export_data_wide", &MCMC::export_data_wide)
        .method("export_data_long", &MCMC::export_data_long)
        .method("summarise_cases", &MCMC::summarise_cases)
        .method("summarise_house_sizes", &MCMC::summarise_house_sizes)
        .method("export_setup", &MCMC::export_setup)
        .method("export_param_list_const", &MCMC::export_param_list_const)
        .method("export_param_list_sample_standard", &MCMC::export_param_list_sample_standard)
        .method("export_param_list_sample_cache", &MCMC::export_param_list_sample_cache)
        .method("export_augm_list_draws", &MCMC::export_augm_list_draws)
        .method("export_augm_list_full", &MCMC::export_augm_list_full)
        .method("export_augm_list_export", &MCMC::export_augm_list_export)
        .method("export_augm_list_discrete", &MCMC::export_augm_list_discrete)
        .method("export_augmented_data_wide", &MCMC::export_augmented_data_wide)
        .method("export_infec_profile", &MCMC::export_infec_profile)
        .method("export_draws", &MCMC::export_draws)
        .method("export_house_log_lik_draws", &MCMC::export_house_log_lik_draws)
        .method("export_house_snapshot", &MCMC::export_house_snapshot);

    Rcpp::class_<Model>("ModelHousepi")
        .constructor()
        .derives<MCMC>("McmcHousepi");
}
