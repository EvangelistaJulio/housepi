/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef _STANDALONE
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Rcpp.h>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/operators.h"
#include "tools/string.h"
#include "tools/read.h"
#include "tools/types.h"

#include "mcmc.h"
#include "measures.h"
#include "posterior.h"
#include "registry.h"
#include "augmented.h"
#include "binomial.h"
#include "data_store.h"
#include "data_wide.h"
#include "data_long.h"
#include "distribution.h"
#include "function.h"
#include "group_rate.h"
#include "infec_profile.h"
#include "parameter.h"
#include "variable.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Stdout Message Function  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::msg
{
    /// Throw an error message.
    /// @cite weylandt_suppressing_2018
    [[noreturn]] void error(const std::string &txt)
    {
        throw Rcpp::exception(txt.c_str(), false);
    }

    /// Print an information message to stdout.
    void info(const std::string &txt) noexcept
    {
        Rcpp::Rcout << txt;
        Rcpp::Rcout.flush();
    }

    /// Display a warning message if not already displayed before.
    /// Based on @cite weylandt_suppressing_2018
    void warning(const std::string &txt)
    {
        if (!is_new_warning(txt))
            return;
        Rf_warningcall(R_NilValue, "%s", txt.c_str());
        add_new_warning(txt);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  *** Housepi Rcpp Base : class ***  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using namespace housepi;

    // Base class to wrap internal hosuepi functions
    class ExposeBase
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        ExposeBase(std::string name, const std::string &family,
                   std::vector<std::string> registry_types)
            : family_(family)
        {
            try
            {
                ptr_ = compose_func_registries(name, family_, registry_types, false);
            }
            catch (const std::exception &e)
            {
                msg::error("Exposing " + family_ + " '" + name + "' failed:\n" + e.what());
            }
        }
        virtual ~ExposeBase() = default;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        void init(SEXP args_, double min, double max, bool is_discr)
        {
            try
            {
                if (args_ != R_NilValue)
                {
                    Rcpp::List args(args_);
                    Rcpp::CharacterVector arg_names(args.size());
                    if (args.hasAttribute("names"))
                        arg_names = args.names();

                    for (R_xlen_t i = 0; i < args.size(); ++i)
                        ptr_->init_attr(
                            Rcpp::as<std::string>(arg_names[i]),
                            Attribute(Rcpp::as<double>(args[i])));
                }

                if (!na::is(min))
                    ptr_->init_attr("min", Attribute(min));
                if (!na::is(max))
                    ptr_->init_attr("max", Attribute(max));

                if (is_discr)
                    ptr_->init_discretization();
                ptr_->init_truncation();
                ptr_->assert_is_callable_explicitly();
            }
            catch (const std::exception &e)
            {
                msg::error(
                    "While initializing " + family_ + " '" +
                    std::string(ptr_->name()) + "':\n" + e.what());
            }
        }

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        Rcpp::NumericVector args() const
        {
            Rcpp::NumericVector output(ptr_->n_args());
            Rcpp::CharacterVector names(ptr_->n_args());
            for (size_t i = 0; i < ptr_->n_args(); ++i)
            {
                output[i] = ptr_->arg(i);
                names[i] = std::string(ptr_->arg_name(i));
            }
            output.names() = std::move(names);
            return output;
        }

    protected:
        std::unique_ptr<Function> ptr_;
        std::string family_;

        template <class Fn>
        Rcpp::NumericVector call_(Fn fn, Rcpp::NumericVector x) const
        {
            try
            {
                const R_xlen_t n = x.size();
                Rcpp::NumericVector y(n);
                for (R_xlen_t i = 0; i < n; ++i)
                    y[i] = fn(x[i]);
                return y;
            }
            catch (const std::exception &e)
            {
                msg::error(
                    "While calling " + family_ + " '" + std::string(ptr_->name()) + "':\n" + e.what());
            }
        }

        template <class Fn>
        Rcpp::NumericVector call_ratio_(Fn fn,
                                        Rcpp::NumericVector xo,
                                        Rcpp::NumericVector xn) const
        {
            try
            {
                const R_xlen_t n = xo.size();
                if (xn.size() != n)
                    throw std::invalid_argument(
                        "Size of 'xo' must be identical to size of 'xn'.");
                Rcpp::NumericVector y(n);
                for (R_xlen_t i = 0; i < n; ++i)
                    y[i] = fn(xo[i], xn[i]);
                return y;
            }
            catch (const std::exception &e)
            {
                msg::error(
                    "While calling " + family_ + " '" + std::string(ptr_->name()) + "':\n" + e.what());
            }
        }
    };

#define MAKE_EXPOSE_METHOD_CALL(name, expr)                 \
    Rcpp::NumericVector name(Rcpp::NumericVector x) const   \
    {                                                       \
        return call_([this](double v) { return expr; }, x); \
    }
#define MAKE_EXPOSE_METHOD_CALL_RATIO(name, expr)                                  \
    Rcpp::NumericVector name(Rcpp::NumericVector xo, Rcpp::NumericVector xn) const \
    {                                                                              \
        return call_ratio_([this](double vo, double vn) { return expr; }, xo, xn); \
    }

    class ExposeFunction : public ExposeBase
    {
    public:
        ExposeFunction(std::string name_)
            : ExposeBase(name_, "function",
                         {"custom", "transcendental", "arithmetic"}) {}

        MAKE_EXPOSE_METHOD_CALL(func, ptr_->func(v))
        MAKE_EXPOSE_METHOD_CALL(log_func, ptr_->log_func(v))
        MAKE_EXPOSE_METHOD_CALL_RATIO(log_ratio, ptr_->log_ratio_safe(vo, vn))
        MAKE_EXPOSE_METHOD_CALL(prim, ptr_->prim(v))
    };

    class ExposeDistribution : public ExposeBase
    {
    public:
        ExposeDistribution(std::string name_)
            : ExposeBase(name_, "distribution",
                         {"continuos", "discrete"}) {}

        MAKE_EXPOSE_METHOD_CALL(pdf, ptr_->pdf(v))
        MAKE_EXPOSE_METHOD_CALL(log_pdf, distr_()->log_pdf_safe(v))
        MAKE_EXPOSE_METHOD_CALL_RATIO(log_ratio, distr_()->log_ratio_safe(vo, vn))
        MAKE_EXPOSE_METHOD_CALL(cdf, ptr_->cdf(v))
        MAKE_EXPOSE_METHOD_CALL(quantile, distr_()->quantile(v))

        /// Sample n random draws. Seed is drawn from R's RNG so that
        /// set.seed() in R controls the result.
        Rcpp::NumericVector rand(int n) const
        {
            try
            {
                if (n < 0)
                    throw std::invalid_argument("'n' must be non-negative.");

                // generate seed from R
                Rcpp::RNGScope scope;
                std::uint32_t seed = static_cast<std::uint32_t>(
                    R::unif_rand() * std::numeric_limits<std::uint32_t>::max());

                rnd_gen prng(seed);
                Rcpp::NumericVector y(n);
                for (int i = 0; i < n; ++i)
                    y[i] = distr_()->rand(prng);
                return y;
            }
            catch (const std::exception &e)
            {
                msg::error(
                    "While sampling " + family_ + " '" +
                    std::string(ptr_->name()) + "':\n" + e.what());
            }
        }

    private:
        Distribution *distr_() const { return static_cast<Distribution *>(ptr_.get()); }
    };

    class ExposeInfecProfile : public ExposeBase
    {
    public:
        ExposeInfecProfile(std::string name_)
            : ExposeBase(name_, "infectivity profile",
                         {"infecprofile", "continuos", "discrete"}) {}

        MAKE_EXPOSE_METHOD_CALL(inst, ptr_->pdf(v))
        MAKE_EXPOSE_METHOD_CALL(cumul, ptr_->cdf(v))
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMC : imports from Rcpp ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Auxiliary converter to R NA real values.
    double hsp_to_rcpp(double x) { return na::is(x) ? NA_REAL : x; }
    double hsp_to_rcpp(const Attribute &attr)
    {
        if (attr.is_empty())
            return NA_REAL;
        if (attr.is_function() && !attr.function()->is_callable_implicitly())
            return NA_REAL;
        return hsp_to_rcpp(attr.value());
    }

    /// Convert R dataframe map of string vectors.
    str_table rdataframe_to_strmap(SEXP input)
    {
        if (Rf_isNull(input))
            return {};

        Rcpp::DataFrame df(input);
        Rcpp::CharacterVector colnames = df.names();
        str_table output;

        for (R_xlen_t i = 0; i < df.size(); ++i)
        {
            Rcpp::CharacterVector col = df[i];
            std::vector<std::string> col_data(col.begin(), col.end());
            output[Rcpp::as<std::string>(colnames[i])] = std::move(col_data);
        }

        return output;
    }

    /// Convert R dataframe map of double vectors.
    dbl_table rdataframe_to_dblmap(SEXP input)
    {
        if (Rf_isNull(input))
            return {};

        Rcpp::DataFrame df(input);
        Rcpp::CharacterVector colnames = df.names();
        dbl_table output;

        for (R_xlen_t i = 0; i < df.size(); ++i)
        {
            Rcpp::NumericVector col = df[i];
            std::vector<double> col_data(col.begin(), col.end());
            output[Rcpp::as<std::string>(colnames[i])] = std::move(col_data);
        }

        return output;
    }

    /// Export Function argument list to R.
    Rcpp::NumericVector rcpp_export_args(const Function *func)
    {
        size_t n_args = func->n_args();
        Rcpp::NumericVector vals(n_args);
        Rcpp::CharacterVector keys(n_args);
        for (size_t k = 0; k < n_args; ++k)
        {
            vals[k] = func->arg(k);
            keys[k] = std::string(func->arg_name(k));
        }
        vals.names() = std::move(keys);
        return vals;
    }

    /// Standaridized percentage formatting.
    double format_percentage(size_t value, size_t total)
    {
        return std::round(10000. * static_cast<double>(value) / static_cast<double>(total)) / 100.;
    }

#define MAKE_RCPP_IMPORT_METHOD(module, section_name)            \
    void import_##module(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)  \
    {                                                            \
        try                                                      \
        {                                                        \
            mcmc_ptr->module.build(rdataframe_to_strmap(input)); \
        }                                                        \
        catch (const std::exception &e)                          \
        {                                                        \
            msg::error(std::string("While importing ") +         \
                       section_name + ":\n" + e.what());         \
        }                                                        \
    }
    MAKE_RCPP_IMPORT_METHOD(setup, "setup")

    void import_param(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->param.build(rdataframe_to_strmap(input), mcmc_ptr->translator);
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing parameters:\n") + e.what());
        }
    }

    void import_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->detach_data_store();
            mcmc_ptr->data_wide().build(rdataframe_to_strmap(input));
            mcmc_ptr->init_equivalence_class_compression();
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing data-wide:\n") + e.what());
        }
    }

    void import_data_long(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            /// FIXME: copy-on-write should be also on long? Not clear...
            mcmc_ptr->data_long().build(mcmc_ptr->data_wide(), rdataframe_to_strmap(input));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing data-long:\n") + e.what());
        }
    }

    void import_group_rate(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->group_rate.build(mcmc_ptr->data_wide(), mcmc_ptr->param, rdataframe_to_strmap(input));
            mcmc_ptr->init_equivalence_class_compression();
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing group rates:\n") + e.what());
        }
    }

    void import_contact_rate(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->contact_rate.build(mcmc_ptr->data_wide(), mcmc_ptr->param, rdataframe_to_strmap(input));
            mcmc_ptr->init_equivalence_class_compression();
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing contact rates:\n") + e.what());
        }
    }

    void import_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->infec_profile.build(mcmc_ptr->translator, rdataframe_to_strmap(input));
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing infectivity profile:\n") + e.what());
        }
    }

    void import_augm(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->augm.build(mcmc_ptr->infec_profile, rdataframe_to_strmap(input));
            mcmc_ptr->init_equivalence_class_compression();
        }
        catch (const std::exception &e)
        {
            msg::error(std::string("While importing augmented:\n") + e.what());
        }
    }

    void import_param_table(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->param.build_for_predic(rdataframe_to_dblmap(input));
        }
        catch (const std::exception &e)
        {
            msg::error("While importanting parameter table:\n" +
                       std::string(e.what()));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMC : exports to Rcpp //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export to R definitions of infectivity profile.
    Rcpp::List helper_export_infec_profile(const InfecProfile &infec_profile)
    {
        if (infec_profile.is_augmented())
            /// FIXME: I need to reverse parser here + cast to numeric
            return Rcpp::List::create(
                Rcpp::_["func"] = Rcpp::CharacterVector(Rcpp::wrap(infec_profile.augm_definitions.at("func"))),
                Rcpp::_["func_args"] = Rcpp::CharacterVector(Rcpp::wrap(infec_profile.augm_definitions.at("func"))),
                Rcpp::_["min"] = Rcpp::CharacterVector(Rcpp::wrap(infec_profile.augm_definitions.at("min"))),
                Rcpp::_["max"] = Rcpp::CharacterVector(Rcpp::wrap(infec_profile.augm_definitions.at("max"))),
                Rcpp::_["discr"] = Rcpp::CharacterVector(Rcpp::wrap(infec_profile.augm_definitions.at("discr"))));

        if (!infec_profile.generic)
            return Rcpp::List::create();

        // NOTE: arguments should not be polymorphic when called here
        return Rcpp::List::create(
            Rcpp::_["func"] = std::string(infec_profile.generic->name()),
            Rcpp::_["func_args"] = rcpp_export_args(infec_profile.generic.get()),
            Rcpp::_["min"] = hsp_to_rcpp(infec_profile.generic->lower()),
            Rcpp::_["max"] = hsp_to_rcpp(infec_profile.generic->upper()),
            Rcpp::_["discr"] = infec_profile.generic->is_discrete());
    }

    /// Export current setup configurations to R List.
    Rcpp::List export_setup(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const Setup &setup = mcmc_ptr->setup;

        return Rcpp::List::create(
            Rcpp::_["steps"] = setup.n_steps,
            Rcpp::_["burn"] = setup.n_burn,
            Rcpp::_["thin"] = setup.n_thin,
            Rcpp::_["adapt"] = setup.n_adapt,
            Rcpp::_["chains"] = setup.n_chains,
            Rcpp::_["threads"] = setup.n_threads,
            Rcpp::_["id_draws"] = Rcpp::CharacterVector(Rcpp::wrap(setup.id_draws)),
            Rcpp::_["seed"] = setup.seed,
            Rcpp::_["progress"] = setup.progress);
    }

    Rcpp::List export_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return helper_export_infec_profile(mcmc_ptr->infec_profile);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  MCMC : DataWide : exports  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Rcpp::List helper_export_data_wide(const DataWide &data_wide,
                                       const std::vector<std::string> &columns,
                                       bool is_status)
    {
        size_t n_rows = data_wide.report_indiv.total;
        Rcpp::CharacterVector id_col = Rcpp::CharacterVector(n_rows);
        Rcpp::CharacterVector house_id_col = Rcpp::CharacterVector(n_rows);
        Rcpp::IntegerVector status_col = Rcpp::IntegerVector(n_rows);

        // other non standard data columns
        std::vector<size_t> list_export = data_wide.khold(columns);
        std::vector<Rcpp::NumericVector> data_col(columns.size());
        for (size_t d = 0; d < columns.size(); ++d)
            data_col[d] = Rcpp::NumericVector(n_rows);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &member : house)
            {
                id_col[row] = member.id();
                house_id_col[row] = house.house_id();
                status_col[row] = member.status();

                for (size_t d = 0; d < columns.size(); ++d)
                    data_col[d][row] = hsp_to_rcpp(member.data(list_export[d]));
                ++row;
            }

        Rcpp::List output;
        output["id"] = id_col;
        output["house_id"] = house_id_col;
        if (is_status)
            output["status"] = status_col;
        for (size_t d = 0; d < columns.size(); ++d)
            output[columns[d]] = data_col[d];

        return output;
    }

    Rcpp::DataFrame export_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        return Rcpp::DataFrame(helper_export_data_wide(data_wide, data_wide.khold().all(), true));
    }

    /// Export to R all augmented variables at the current state.
    Rcpp::DataFrame export_augmented_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const Augmented &augm = mcmc_ptr->augm;

        size_t n_cols = augm.list_export.size();
        size_t n_rows = data_wide.report_indiv.total;

        std::vector<Rcpp::NumericVector> augm_cols(n_cols);
        for (size_t i = 0; i < n_cols; ++i)
            augm_cols[i] = Rcpp::NumericVector(n_rows);

        const std::vector<std::string> data_names = data_wide.khold().all() -
                                                    augm.khold(augm.list_export);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &member : house)
            {
                size_t c = 0;
                for (size_t a : augm.list_export)
                    augm_cols[c++][row] = hsp_to_rcpp(mcmc_ptr->get_augm_varia(a, member));
                ++row;
            }

        Rcpp::List output = helper_export_data_wide(data_wide, data_names, true);
        size_t i = 0;
        for (size_t a : augm.list_export)
            output[augm.khold(a)] = augm_cols[i++];
        return Rcpp::DataFrame(output);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  MCMCBase : DataLong : exports  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export current data configurations to R DataFrame.
    Rcpp::DataFrame export_data_long(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const DataLong &data_long = mcmc_ptr->data_long();

        if (data_long.indiv_datalist.size() == 0)
            return Rcpp::DataFrame();

        // retrieve total number of rows
        size_t n_rows = 0;
        for (const auto &tbl : data_long.indiv_datalist)
            if (!tbl.vec.empty())
                n_rows += tbl.vec.begin()->size();

        Rcpp::CharacterVector datalist_index(n_rows);

        // other non standard data-long columns
        std::vector<Rcpp::NumericVector> datalist_dbl(data_long.khold().size());
        for (size_t d = 0; d < datalist_dbl.size(); ++d)
            datalist_dbl[d] = Rcpp::NumericVector(n_rows);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &member : house)
            {
                const std::vector<std::vector<double>> *ptr = member.get_data_long_ptr();
                if (ptr && !ptr->empty())
                    for (size_t k = 0; k < ptr->begin()->size(); ++k)
                    {
                        datalist_index[row] = member.id();
                        for (size_t d = 0; d < datalist_dbl.size(); ++d)
                            datalist_dbl[d][row] = hsp_to_rcpp(ptr->at(d)[k]);
                        ++row;
                    }
            }

        // assign vectors to dataframe
        Rcpp::List output;
        output["id"] = datalist_index;
        for (size_t d = 0; d < datalist_dbl.size(); ++d)
            output.push_back(datalist_dbl[d], data_long.khold(d));

        return Rcpp::DataFrame(output);
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  MCMCBase : report exports  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export current data configurations to R DataFrame.
    Rcpp::DataFrame summarise_house_cases(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const ReportIndividual &report_indiv = data_wide.report_indiv;
        const ReportHousehold &report_house = data_wide.report_house;

        Rcpp::CharacterVector level = Rcpp::wrap(std::vector<std::string>(4, "individual") +
                                                 std::vector<std::string>(5, "household"));
        Rcpp::CharacterVector metric = {"total", "missing", "not infected", "infected",
                                        "total", "distinct", "no cases", "at least 1 case", "at least 2 cases"};
        Rcpp::IntegerVector count(level.size());
        Rcpp::NumericVector percent(level.size());

        // individual total
        count[0] = report_indiv.total;
        percent[0] = 100.;
        // individual missing
        count[1] = report_indiv.missing;
        percent[1] = format_percentage(report_indiv.missing, report_indiv.total);
        // individual not infected
        count[2] = report_indiv.susceptible;
        percent[2] = format_percentage(report_indiv.susceptible, report_indiv.total);
        // individual infected
        count[3] = report_indiv.infectious;
        percent[3] = format_percentage(report_indiv.infectious, report_indiv.total);
        // household total
        count[4] = report_house.total;
        percent[4] = 100.;
        // household equivalence classes
        count[5] = report_house.distinct;
        percent[5] = 100.;
        // household no cases
        count[6] = report_house.no_cases;
        percent[6] = format_percentage(report_house.no_cases, report_house.total);
        // household at least 1 case
        count[7] = report_house.at_least_1_case;
        percent[7] = format_percentage(report_house.at_least_1_case, report_house.total);
        // household at least 2 cases
        count[8] = report_house.at_least_2_cases;
        percent[8] = format_percentage(report_house.at_least_2_cases, report_house.total);

        return Rcpp::DataFrame::create(
            Rcpp::_["level"] = std::move(level),
            Rcpp::_["metric"] = std::move(metric),
            Rcpp::_["count"] = std::move(count),
            Rcpp::_["percent"] = std::move(percent));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  MCMCBase : export house size report  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export current data configurations to R DataFrame.
    Rcpp::DataFrame summarise_house_sizes(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const ReportIndividual &report_indiv = data_wide.report_indiv;
        const ReportHouseSizes &report_hh_size = data_wide.report_hh_size;

        Rcpp::CharacterVector metric({"most cases", "most missing", "smallest size", "largest size"});
        R_xlen_t n_rows = metric.size();
        Rcpp::IntegerVector value(n_rows), count(n_rows);
        Rcpp::NumericVector percent(n_rows);

        value[0] = report_hh_size.max_cases.first;
        count[0] = report_hh_size.max_cases.second;
        percent[0] = format_percentage(count[0], report_indiv.total);

        value[1] = report_hh_size.max_miss.first;
        count[1] = report_hh_size.max_miss.second;
        percent[1] = format_percentage(count[1], report_indiv.total);

        value[2] = report_hh_size.min_size.first;
        count[2] = report_hh_size.min_size.second;
        percent[2] = format_percentage(count[2], report_indiv.total);

        value[3] = report_hh_size.max_size.first;
        count[3] = report_hh_size.max_size.second;
        percent[3] = format_percentage(count[3], report_indiv.total);

        return Rcpp::DataFrame::create(
            Rcpp::_["metric"] = std::move(metric),
            Rcpp::_["value"] = std::move(value),
            Rcpp::_["count"] = std::move(count),
            Rcpp::_["percent"] = std::move(percent));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Posterior : Rcpp export methods  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Move a vector<double> into a fresh NumericVector and release the source.
    /// One allocation, one copy, no enclosing matrix.
    Rcpp::NumericVector flush_vector_to_numvec(std::vector<double> &src)
    {
        Rcpp::NumericVector col(src.begin(), src.end());
        std::vector<double>().swap(src);
        return col;
    }

    /// Attach data.frame class + row.names to a List of equal-length columns.
    void helper_set_dataframe_attrs(Rcpp::List &cols, R_xlen_t n_rows)
    {
        cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(n_rows));
        cols.attr("class") = "data.frame";
    }

    /// Flush each vector<double> as a NumericVector column of a data.frame,
    /// release all C++ storage, set column names. Skips the contiguous
    /// NumericMatrix intermediate so R can use the columns in place.
    Rcpp::List flush_2d_draws_to_df(
        std::vector<std::vector<double>> &draws,
        const std::vector<std::string> &col_names,
        size_t n_rows,
        const char *err_type)
    {
        if (draws.empty())
            return Rcpp::List();

        size_t n_cols = draws.size();
        Rcpp::List cols(n_cols);

        for (size_t c = 0; c < n_cols; ++c)
            cols[c] = flush_vector_to_numvec(draws[c]);
        std::vector<std::vector<double>>().swap(draws);

        cols.attr("names") = Rcpp::wrap(col_names);
        helper_set_dataframe_attrs(cols, static_cast<R_xlen_t>(n_rows));
        return cols;
    }

    /// Flush total transmission log-likelihood draws as a data.frame.
    Rcpp::List flush_tot_log_lik_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        return flush_2d_draws_to_df(post.tot_transm_log_lik,
                                    post.tot_transm_log_lik_names,
                                    post.count(), "total log-likelihood");
    }

    /// Flush total augmented-variable log-prior draws as a data.frame.
    /// Empty when the model has no augmented prior.
    Rcpp::List flush_tot_augm_log_prior_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        if (post.tot_augm_log_prior.empty())
            return Rcpp::List();
        return flush_2d_draws_to_df(post.tot_augm_log_prior,
                                    post.tot_augm_log_prior_names,
                                    post.count(), "total augmented log-prior");
    }

    /// Rcpp wrapper: parameter value draws as a data.frame.
    Rcpp::List flush_param_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        return flush_2d_draws_to_df(post.param_values, post.param_names,
                                    post.count(), "parameter");
    }

    /// Rcpp wrapper: parameter log-prior draws as a data.frame.
    Rcpp::List flush_param_log_prior_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        return flush_2d_draws_to_df(post.param_priors, post.param_prior_names,
                                    post.count(), "parameter log-prior");
    }

    /// Rcpp wrapper: augmented variable draws as a data.frame.
    Rcpp::List flush_augm_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;

        if (post.augm_values.empty())
            return Rcpp::List();

        size_t n_cols = post.augm_col_names.size();
        Rcpp::List cols(n_cols);

        size_t col = 0;
        for (auto &member_matrix : post.augm_values)
        {
            for (auto &var_vec : member_matrix)
                cols[col++] = flush_vector_to_numvec(var_vec);
            std::vector<std::vector<double>>().swap(member_matrix);
        }
        std::vector<std::vector<std::vector<double>>>().swap(post.augm_values);

        cols.attr("names") = Rcpp::wrap(post.augm_col_names);
        helper_set_dataframe_attrs(cols, static_cast<R_xlen_t>(post.count()));
        return cols;
    }

    /// Rcpp wrapper: household log-likelihood draws as a data.frame.
    /// Households of same equivalence class share NumericVector SEXP.
    /// R's reference counting keeps storage at n_classes * n_draws
    /// until a downstream consumer (e.g. as_draws_matrix) materialises it.
    Rcpp::List flush_house_log_lik_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        const DataWide &data_wide = mcmc_ptr->data_wide();

        size_t n_classes = post.house_transm_log_lik.size();
        size_t n_houses = data_wide.report_house.total;

        // One NumericVector per class; detached copy, C++ source freed.
        std::vector<Rcpp::NumericVector> class_cols(n_classes);
        for (size_t R = 0; R < n_classes; ++R)
            class_cols[R] = flush_vector_to_numvec(post.house_transm_log_lik[R]);
        std::vector<std::vector<double>>().swap(post.house_transm_log_lik);

        // Identity case (n_classes == n_houses) is naturally a no-op.
        Rcpp::List cols(n_houses);
        for (size_t H = 0; H < n_houses; ++H)
            cols[H] = class_cols[data_wide.get_eq_class(H)];

        cols.attr("names") = Rcpp::wrap(post.house_names);
        helper_set_dataframe_attrs(cols, static_cast<R_xlen_t>(post.count()));
        return cols;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMC : house state snapshot /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export to R a snapshot of current state of house data.
    /// NOTE: there are not private parameters so far.
    Rcpp::NumericVector helper_param_rcpp_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const Parameter &param = mcmc_ptr->param;

        if (param.khold().empty())
            return Rcpp::NumericVector::create();

        size_t n = param.khold().size();
        Rcpp::NumericVector output(n);
        Rcpp::CharacterVector headers(n);
        for (size_t p = 0; p < n; ++p)
        {
            output[p] = param.at(p).value();
            headers[p] = std::string(param.at(p).name());
        }
        output.names() = Rcpp::CharacterVector(headers);

        return output;
    }

    /// Export to R a snapshot of current state of house data.
    Rcpp::DataFrame helper_data_wide_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const std::vector<std::string> &headers = data_wide.khold().all();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);
        Rcpp::CharacterVector house_id_col(n_rows);
        Rcpp::IntegerVector status_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(headers.size());
        for (size_t d = 0; d < headers.size(); ++d)
            df_double[d] = Rcpp::NumericVector(n_rows);

        for (const Individual &member : house)
        {
            size_t i = member.hhindex();
            id_col[i] = member.id();
            house_id_col[i] = house.house_id();
            status_col[i] = member.status();

            for (size_t d = 0; d < headers.size(); ++d)
                df_double[d][i] = hsp_to_rcpp(member.data(d));
        }

        auto output = Rcpp::DataFrame::create(
            Rcpp::_["id"] = std::move(id_col),
            Rcpp::_["house_id"] = std::move(house_id_col),
            Rcpp::_["status"] = std::move(status_col));
        for (size_t d = 0; d < headers.size(); ++d)
            output.push_back(df_double[d], headers[d]);

        return output;
    }

    /// Export to R a snapshot of current state of house data-long.
    Rcpp::DataFrame helper_data_long_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
         const DataLong &data_long = mcmc_ptr->data_long();

        if (data_long.indiv_datalist.empty())
            return Rcpp::DataFrame::create();

        size_t n_rows = 0;
        for (const Individual &member : house)
        {
            const std::vector<std::vector<double>> *ptr = member.get_data_long_ptr();
            if (ptr && !ptr->empty())
                n_rows += ptr->begin()->size();
        }

        Rcpp::CharacterVector datalist_id(n_rows);
        std::vector<Rcpp::NumericVector> df_double(data_long.khold().size());
        for (size_t d = 0; d < df_double.size(); ++d)
            df_double[d] = Rcpp::NumericVector(n_rows);

        size_t r = 0;
        for (const Individual &member : house)
        {
            const std::vector<std::vector<double>> *ptr = member.get_data_long_ptr();
            if (ptr && !ptr->empty())
                for (size_t k = 0; k < ptr->begin()->size(); ++k)
                {
                    for (size_t d = 0; d < df_double.size(); ++d)
                        df_double[d][r] = hsp_to_rcpp(ptr->at(d)[k]);
                    datalist_id[r++] = member.id();
                }
        }

        // assign vectors to dataframe
        Rcpp::List output;
        output["id"] = datalist_id;
        for (size_t d = 0; d < df_double.size(); ++d)
            output[data_long.khold(d)] = df_double[d];

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_group_rates_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const GroupRate &group_rate = mcmc_ptr->group_rate;
        const std::vector<std::string> &headers = group_rate.khold().all();

        if (headers.empty())
            return Rcpp::DataFrame::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(headers.size());
        for (size_t c = 0; c < headers.size(); ++c)
            df_double[c] = Rcpp::NumericVector(n_rows);

        for (const Individual &member : house)
        {
            size_t i = member.hhindex();
            id_col[i] = member.id();

            for (size_t c = 0; c < headers.size(); ++c)
                df_double[c][i] = hsp_to_rcpp(member.get_group_rate_attribute(c));
        }

        Rcpp::List output;
        output["id"] = id_col;
        for (size_t c = 0; c < headers.size(); ++c)
            output.push_back(df_double[c], headers[c]);

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current contact rate attribute values
    /// for a household as a (from x to) numeric matrix with member id dimnames.
    /// Returns an empty 0x0 matrix when no contact rate was defined.
    Rcpp::NumericMatrix helper_contact_rates_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const ContactRate &contact_rate = mcmc_ptr->contact_rate;

        if (contact_rate.is_empty())
            return Rcpp::NumericMatrix(0, 0);

        const size_t n = house.size();
        Rcpp::NumericMatrix mtx(n, n);
        Rcpp::CharacterVector ids(n);

        for (size_t i = 0; i < n; ++i)
        {
            ids[i] = house.member(i).id();
            for (size_t j = 0; j < n; ++j)
                mtx(i, j) = hsp_to_rcpp(house.contact_rate(i, j));
        }

        mtx.attr("dimnames") = Rcpp::List::create(ids, ids);
        return mtx;
    }

    /// Export to R a snapshot of current state of house infectivity profile.
    Rcpp::List helper_infec_profile_rcpp_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                  const Household &house)
    {
        const InfecProfile &infec_profile = mcmc_ptr->infec_profile;

        if (!infec_profile.generic && !infec_profile.is_augmented())
            return Rcpp::List::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        Rcpp::CharacterVector infprof_name(n_rows, NA_STRING);
        Rcpp::LogicalVector infprof_discr(n_rows, NA_LOGICAL);
        Rcpp::List infprof_args(n_rows, Rcpp::NumericVector(0));
        Rcpp::NumericVector infprof_min(n_rows, NA_REAL);
        Rcpp::NumericVector infprof_max(n_rows, NA_REAL);

        for (size_t i = 0; i < house.size(); ++i)
            id_col[i] = house.member(i).id();

        for (size_t i : house.infected_hhindexes())
            if (infec_profile.is_augmented())
            {
                const auto &f = mcmc_ptr->get_augm_varia(infec_profile.augm_key, house.member(i)).function();

                infprof_name[i] = std::string(f->name());
                infprof_discr[i] = f->is_discrete();
                infprof_args[i] = rcpp_export_args(f);
                infprof_min[i] = hsp_to_rcpp(f->lower());
                infprof_max[i] = hsp_to_rcpp(f->upper());
            }
            else
            {
                Rcpp::List generic = helper_export_infec_profile(infec_profile);
                infprof_name[i] = Rcpp::as<std::string>(generic["func"]);
                infprof_discr[i] = Rcpp::as<bool>(generic["discr"]);
                infprof_args[i] = generic["func_args"];
                infprof_min[i] = Rcpp::as<double>(generic["min"]);
                infprof_max[i] = Rcpp::as<double>(generic["max"]);
            }

        return Rcpp::List::create(
            Rcpp::_["id"] = id_col,
            Rcpp::_["discr"] = infprof_discr,
            Rcpp::_["func"] = infprof_name,
            Rcpp::_["func_args"] = infprof_args,
            Rcpp::_["min"] = infprof_min,
            Rcpp::_["max"] = infprof_max);
    }

    Rcpp::DataFrame helper_augm_log_prior_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const Augmented &augm = mcmc_ptr->augm;

        size_t n_cols = augm.list_export.size();

        if (n_cols == 0)
            return Rcpp::DataFrame::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(n_cols);
        for (size_t c = 0; c < n_cols; ++c)
            df_double[c] = Rcpp::NumericVector(n_rows, NA_REAL);

        for (size_t i = 0; i < house.size(); ++i)
        {
            const Individual &member_i = house.member(i);
            id_col[i] = member_i.id();

            size_t c = 0;
            for (size_t a : augm.list_export)
            {
                const Variable &avar = mcmc_ptr->get_augm_varia(a, member_i);
                if (avar.has_prior())
                    df_double[c++][i] = hsp_to_rcpp(avar.get_log_prior());
            }
        }

        Rcpp::List output;
        output["id"] = id_col;
        size_t c = 0;
        for (size_t a : augm.list_export)
            output.push_back(df_double[c++], augm.khold(a));

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_augm_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const Augmented &augm = mcmc_ptr->augm;
        size_t n_cols = augm.list_export.size();

        if (n_cols == 0)
            return Rcpp::DataFrame::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(n_cols);
        for (size_t c = 0; c < n_cols; ++c)
            df_double[c] = Rcpp::NumericVector(n_rows);

        for (size_t i = 0; i < house.size(); ++i)
        {
            const Individual &member_i = house.member(i);
            id_col[i] = member_i.id();
            size_t c = 0;
            for (size_t a : augm.list_export)
                df_double[c++][i] = hsp_to_rcpp(mcmc_ptr->get_augm_varia(a, member_i));
        }

        Rcpp::List output;
        output["id"] = id_col;
        size_t i = 0;
        for (size_t a : augm.list_export)
            output.push_back(df_double[i++], augm.khold(a));

        return Rcpp::DataFrame(output);
    }

    /// Export a list containing all relevant information for a debug.
    Rcpp::List export_house_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr, const std::string &house_id)
    {
        try
        {
            if (mcmc_ptr->data_wide().is_empty())
                return Rcpp::List::create();

            const Household &house = mcmc_ptr->data_wide().at(house_id);

            return Rcpp::List::create(
                Rcpp::_["transm_log_lik"] = house.transm_log_lik(),
                Rcpp::_["augm_log_prior"] = helper_augm_log_prior_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["param"] = helper_param_rcpp_snapshot(mcmc_ptr),
                Rcpp::_["data_wide"] = helper_data_wide_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["data_long"] = helper_data_long_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["group_rate"] = helper_group_rates_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["contact_rate"] = helper_contact_rates_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["infec_profile"] = helper_infec_profile_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["augm"] = helper_augm_rcpp_snapshot(mcmc_ptr, house));
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting snapshot of household '" + house_id + "':\n" +
                       std::string(e.what()));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////  MCMC : simulate ///////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace simulate
    {
        /// Perform a single repetition.
        void single_repetition(MCMC &engine,
                               Rcpp::IntegerVector &status_col,
                               std::vector<Rcpp::NumericVector> &augm_cols,
                               const std::vector<size_t> &row_offsets)
        {
            DataWide &data_wide = engine.data_wide();
            const Augmented &augm = engine.augm;

            housepi::parallel::for_each(
                data_wide.report_house.total,
                [&](size_t H)
                {
                    try
                    {
                        Household &house = data_wide.at(H);
                        engine.model->simulate_single_house(house);

                        size_t row = row_offsets[H]; // thread-local: no shared counter
                        for (const Individual &member : house)
                        {
                            status_col[row] = member.status();
                            size_t c = 0;
                            for (size_t a : augm.list_export)
                                augm_cols[c++][row] = hsp_to_rcpp(engine.get_augm_varia(a, member));
                            ++row;
                        }
                    }
                    catch (const std::exception &e)
                    {
                        throw std::runtime_error(
                            "While simulating " + data_wide.at(H).err_name() + ":\n" +
                            std::string(e.what()));
                    }
                });
        }

        // pre-compute per-household row offsets (fixes OpenMP race on row counter)
        Rcpp::DataFrame build_simulated(const Rcpp::List &shared_df,
                                        const Rcpp::IntegerVector &status_col,
                                        const std::vector<Rcpp::NumericVector> &augm_cols,
                                        const Rcpp::CharacterVector &sim_names)
        {
            // build simulated list
            const size_t n_shared(shared_df.size());
            const size_t n_augm_cols(augm_cols.size());
            Rcpp::List sim_list(n_shared + 1 + n_augm_cols);
            for (size_t i = 0; i < n_shared; ++i)
                sim_list[i] = shared_df[i];
            sim_list[n_shared] = status_col;
            for (size_t c = 0; c < n_augm_cols; ++c)
                sim_list[n_shared + 1 + static_cast<int>(c)] = augm_cols[c];
            sim_list.attr("names") = sim_names;
            return Rcpp::DataFrame(sim_list);
        }
    }

    /// Perform simulation and return simulated data with separate wide and long post-processing.
    Rcpp::List run_simulation(Rcpp::XPtr<MCMC> mcmc_ptr,
                              Rcpp::Nullable<Rcpp::Function> obs_process,
                              Rcpp::Nullable<Rcpp::Function> apply_wide,
                              Rcpp::Nullable<Rcpp::Function> obs_process_and_apply_long)
    {
        MCMC &engine = *mcmc_ptr;
        Setup &setup = engine.setup;
        const Augmented &augm = engine.augm;
        const DataWide &data_wide = engine.data_wide();

        try
        {
            engine.model->init_model_for_simulation(); // model specific initialization

            //////////////////////////////////////////////
            //////  Observation Process & Measures  //////
            //////////////////////////////////////////////

            // UPDATE: only work-around I found for NULL functions...
            std::unique_ptr<Rcpp::Function> f_obs_process = nullptr, f_wide = nullptr, f_long = nullptr;

            if (!obs_process.isNull())
                f_obs_process = std::make_unique<Rcpp::Function>(obs_process);
            if (!apply_wide.isNull())
                f_wide = std::make_unique<Rcpp::Function>(apply_wide);
            if (!obs_process_and_apply_long.isNull())
                f_long = std::make_unique<Rcpp::Function>(obs_process_and_apply_long);

            //////////////////////////////
            //////  Initialization  //////
            //////////////////////////////

            Rcpp::List out_wide(setup.n_steps); // allocate memory for wide results
            Rcpp::List out_long;                // long results are optional
            if (f_long)
                out_long = Rcpp::List(setup.n_steps);

            const std::vector<std::string> static_data_names = engine.get_static_data_names();
            const Rcpp::List shared_df = helper_export_data_wide(data_wide, static_data_names, false);

            // simulated data column name vector
            Rcpp::CharacterVector sim_names = Rcpp::wrap(engine.get_col_names());
            // row offsets for house parallelization
            const std::vector<size_t> row_offsets = engine.get_row_offsets();

            // pre-allocate mutable output columns once (reused across iterations)
            Rcpp::IntegerVector status_col(data_wide.report_indiv.total);
            std::vector<Rcpp::NumericVector> augm_cols(augm.list_export.size());
            for (size_t c = 0; c < augm_cols.size(); ++c)
                augm_cols[c] = Rcpp::NumericVector(data_wide.report_indiv.total);

            ///////////////////////////
            //////  Simulations  //////
            ///////////////////////////

            setup.prog_open_for_simulation();
            for (size_t r = 0; r < setup.n_steps; ++r, setup.prog_update())
            {
                engine.augm.clear_all();
                engine.data_wide().recover_all_status();
                engine.param.load_sample_values(r);

                // run single repetition
                simulate::single_repetition(
                    engine, status_col, augm_cols, row_offsets);

                // build simulated list
                Rcpp::DataFrame simulated(simulate::build_simulated(
                    shared_df, status_col, augm_cols, sim_names));

                // observation process
                if (f_obs_process)
                    simulated = (*f_obs_process)(simulated);

                // apply measure to data-wide or store it
                if (f_wide)
                    out_wide[r] = (*f_wide)(simulated);
                else
                    out_wide[r] = Rcpp::clone(simulated); // clone when storing
                // long observation process and apply measure
                if (f_long)
                    out_long[r] = (*f_long)(simulated);
            }
            setup.prog_close();

            // restore observed status so post-simulation reads of data_wide
            // return the original data (e.g. running observed measures again).
            engine.augm.clear_all();
            engine.data_wide().recover_all_status();

            return Rcpp::List::create(Rcpp::_["wide"] = out_wide,
                                      Rcpp::_["long"] = out_long);
        }
        catch (std::exception &e)
        {
            setup.crash(e);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  MCMC : log-evidence /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /* --FUTURE-LOG-EVIDENCE--
    /// Estimate log_evidence according to a given scheme and number of rungs.
    /// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
    /// Based on @cite friel_improving_2014
    Rcpp::List log_evidence(Rcpp::XPtr<MCMC> mcmc_ptr, const std::string &scheme, unsigned long n_rungs)
    {
        MCMC &engine = *mcmc_ptr;
        const auto pp = engine.power_posteriors(scheme, n_rungs);
        const std::vector<double> &inv_tem = pp[0];
        const std::vector<double> &avg_vec = pp[1];
        const std::vector<double> &var_vec = pp[2];

        double output = 0;
        for (unsigned long k = 0; k < engine.setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return Rcpp::List::create(
            Rcpp::_["estimate"] = output,
            Rcpp::_["scheme"] = scheme,
            Rcpp::_["n_rungs"] = engine->setup.c_sample_size,
            Rcpp::_["inv_temper"] = Rcpp::wrap(inv_tem),
            Rcpp::_["mean_log_lik"] = Rcpp::wrap(avg_vec),
            Rcpp::_["var_log_lik"] = Rcpp::wrap(var_vec));
    }
    */

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  *** Rcpp exports ***  ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_EXPORT_CHARACTER_VECTOR(method)                                \
    Rcpp::CharacterVector method(Rcpp::XPtr<housepi::ModelBase> model_xptr) \
    {                                                                       \
        return Rcpp::wrap(model_xptr->method());                            \
    }

    MAKE_EXPORT_CHARACTER_VECTOR(get_DataKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_LongKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_ParamKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_GroupRateKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_AugmKey)

    // UPDATE: I am just using this to check flags
    Rcpp::XPtr<ModelBase> new_model(std::string module_name)
    {
        try
        {
            auto model = housepi::ModelRegistry::instance().create(module_name);
            return Rcpp::XPtr<ModelBase>(model.release(), true);
        }
        catch (const std::exception &e)
        {
            msg::error("While creating new model: " + std::string(e.what()));
        }
    }

    Rcpp::XPtr<MCMC> create_mcmc_runner(std::string module_name)
    {
        Rcpp::XPtr<MCMC> mcmc(new MCMC(module_name), true);

        // options to be passed to progress setup
        const unsigned width = Rcpp::as<unsigned>(Rcpp::Function("getOption")("width"));
        const bool is_interactive = Rcpp::as<bool>(Rcpp::Function("interactive")());
        mcmc->setup.init_R_context(
            std::min(width, 78u) / 2,
            is_interactive,
            []
            { Rcpp::checkUserInterrupt(); });

        return mcmc;
    }

    /// Share the DataStore between two MCMC objects (COW copy).
    void share_data_store(Rcpp::XPtr<MCMC> dest, Rcpp::XPtr<MCMC> src)
    {
        dest->set_data_store(src->pDataStore);
    }

    bool has_contact_rates(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return mcmc_ptr->model->has_contact_rates();
    }

    Rcpp::NumericVector evaluate_log_posterior(Rcpp::XPtr<MCMC> mcmc_ptr, std::string house_id)
    {
        try
        {
            // ordered terms by convention; R assigns the names
            return Rcpp::wrap(mcmc_ptr->evaluate_log_posterior(house_id));
        }
        catch (const std::exception &e)
        {
            mcmc_ptr->setup.crash(e);
        }
        return Rcpp::NumericVector();
    }
    void run_inference(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        try
        {
            mcmc_ptr->run_inference();
        }
        catch (std::exception &e)
        {
            mcmc_ptr->setup.crash(e);
        }
    }

    /// Aggregate individuals by house_id. NA status becomes negative (NA_INTEGER is INT_MIN).
    measure::cases::Report summarise_houses_from_df(Rcpp::DataFrame data_wide)
    {
        Rcpp::CharacterVector col_names = data_wide.names();

        std::vector<std::string> house_id = {};
        if (std::find(col_names.begin(), col_names.end(), "house_id") != col_names.end())
            house_id = Rcpp::as<std::vector<std::string>>(data_wide["house_id"]);

        std::vector<int> status = {};
        if (std::find(col_names.begin(), col_names.end(), "status") != col_names.end())
            status = Rcpp::as<std::vector<int>>(data_wide["status"]);

        return measure::summarise_houses(house_id, status);
    }

    Rcpp::DataFrame measure_sar_by_house_size(Rcpp::DataFrame data_wide)
    {
        measure::SecAttackRates res = measure::sar_by_house_size(summarise_houses_from_df(data_wide));
        unsigned grid_len = static_cast<unsigned>(res.sar_mean.size());

        Rcpp::IntegerVector out_size(grid_len);
        Rcpp::NumericVector out_sar(grid_len);
        for (unsigned i = 0; i < grid_len; ++i)
        {
            out_size[i] = res.grid_min + i;
            out_sar[i] = res.count[i] > 0 ? res.sar_mean[i] : NA_REAL;
        }
        return Rcpp::DataFrame::create(
            Rcpp::_["house_size"] = out_size,
            Rcpp::_["sar"] = out_sar);
    }

    Rcpp::DataFrame measure_cases_by_house_size(Rcpp::DataFrame data_wide)
    {
        measure::cases::Table res = measure::cases_by_house_size(summarise_houses_from_df(data_wide));
        unsigned n_rows = static_cast<unsigned>(res.counts.size());

        Rcpp::IntegerVector out_size(n_rows);
        for (unsigned i = 0; i < n_rows; ++i)
            out_size[i] = res.min_size + i;

        Rcpp::List result;
        result.push_back(out_size, "house_size");

        unsigned n_case_cols = res.max_size + 1;
        for (unsigned c = 0; c < n_case_cols; ++c)
        {
            Rcpp::IntegerVector col(n_rows);
            for (unsigned r = 0; r < n_rows; ++r)
                col[r] = (c > res.min_size + r) ? NA_INTEGER : res.counts[r][c];

            result.push_back(col, std::to_string(c));
        }
        Rcpp::DataFrame df(result);
        df.attr("names") = result.names();
        return df;
    }

    Rcpp::DataFrame count_non_zero_cases(Rcpp::DataFrame data_wide)
    {
        measure::cases::NonZero res = measure::non_zero_cases(summarise_houses_from_df(data_wide));
        unsigned grid_len = static_cast<unsigned>(res.count.size());

        Rcpp::IntegerVector out_size(grid_len);
        Rcpp::IntegerVector out_count(grid_len);
        for (unsigned i = 0; i < grid_len; ++i)
        {
            out_size[i] = 2 + i;
            out_count[i] = res.count[i];
        }
        return Rcpp::DataFrame::create(
            Rcpp::_["house_size"] = out_size,
            Rcpp::_["count"] = out_count);
    }
}

RCPP_MODULE(housepi)
{
    Rcpp::class_<ExposeBase>("ExposeBase")
        .method("init", &ExposeBase::init)
        .method("args", &ExposeBase::args);

    Rcpp::class_<ExposeFunction>("ExposeFunction")
        .constructor<std::string>()
        .derives<ExposeBase>("ExposeBase")
        .method("func", &ExposeFunction::func)
        .method("log_func", &ExposeFunction::log_func)
        .method("log_ratio", &ExposeFunction::log_ratio)
        .method("prim", &ExposeFunction::prim);

    Rcpp::class_<ExposeDistribution>("ExposeDistribution")
        .constructor<std::string>()
        .derives<ExposeBase>("ExposeBase")
        .method("pdf", &ExposeDistribution::pdf)
        .method("log_pdf", &ExposeDistribution::log_pdf)
        .method("log_ratio", &ExposeDistribution::log_ratio)
        .method("cdf", &ExposeDistribution::cdf)
        .method("quantile", &ExposeDistribution::quantile)
        .method("rand", &ExposeDistribution::rand);

    Rcpp::class_<ExposeInfecProfile>("ExposeInfecProfile")
        .constructor<std::string>()
        .derives<ExposeBase>("ExposeBase")
        .method("inst", &ExposeInfecProfile::inst)
        .method("cumul", &ExposeInfecProfile::cumul);

    Rcpp::function("get_DataKey", &get_DataKey);
    Rcpp::function("get_LongKey", &get_LongKey);
    Rcpp::function("get_ParamKey", &get_ParamKey);
    Rcpp::function("get_GroupRateKey", &get_GroupRateKey);
    Rcpp::function("get_AugmKey", &get_AugmKey);

    Rcpp::function("measure_sar_by_house_size", &measure_sar_by_house_size);
    Rcpp::function("measure_cases_by_house_size", &measure_cases_by_house_size);
    Rcpp::function("count_non_zero_cases", &count_non_zero_cases);

    Rcpp::function("auxiliary_chain_binomial_debugger",
                   &housepi::chbi::auxiliary_chain_binomial_debugger);

    Rcpp::function("new_model", &new_model);
    Rcpp::function("create_mcmc_runner", &create_mcmc_runner);
    Rcpp::function("run_inference", &run_inference);
    Rcpp::function("run_simulation", &run_simulation);
    Rcpp::function("evaluate_log_posterior", &evaluate_log_posterior);
    Rcpp::function("import_data_wide", &import_data_wide);
    Rcpp::function("import_data_long", &import_data_long);
    Rcpp::function("share_data_store", &share_data_store);
    Rcpp::function("import_group_rate", &import_group_rate);
    Rcpp::function("import_contact_rate", &import_contact_rate);
    Rcpp::function("import_setup", &import_setup);
    Rcpp::function("import_param", &import_param);
    Rcpp::function("import_infec_profile", &import_infec_profile);
    Rcpp::function("import_augm", &import_augm);
    Rcpp::function("import_param_table", &import_param_table);
    Rcpp::function("has_contact_rates", &has_contact_rates);
    Rcpp::function("export_data_wide", &export_data_wide);
    Rcpp::function("export_data_long", &export_data_long);
    Rcpp::function("summarise_house_cases", &summarise_house_cases);
    Rcpp::function("summarise_house_sizes", &summarise_house_sizes);
    Rcpp::function("export_setup", &export_setup);
    Rcpp::function("export_augmented_data_wide", &export_augmented_data_wide);
    Rcpp::function("export_infec_profile", &export_infec_profile);
    Rcpp::function("flush_tot_log_lik_draws", &flush_tot_log_lik_draws);
    Rcpp::function("flush_tot_augm_log_prior_draws", &flush_tot_augm_log_prior_draws);
    Rcpp::function("flush_param_draws", &flush_param_draws);
    Rcpp::function("flush_param_log_prior_draws", &flush_param_log_prior_draws);
    Rcpp::function("flush_augm_draws", &flush_augm_draws);
    Rcpp::function("flush_house_log_lik_draws", &flush_house_log_lik_draws);
    Rcpp::function("export_house_snapshot", &export_house_snapshot);
}
#endif
