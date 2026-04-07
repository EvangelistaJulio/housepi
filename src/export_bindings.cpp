/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef _STANDALONE
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "tools/atomic.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/types.h"

#include "mcmc.h"
#include "registry.h"
#include "augmented.h"
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

namespace hsp::msg
{
    /// Throw an error message.
    /// Ref:
    /// -----
    /// Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
    /// Michael Weylandt — written Mar 27, 2018
    /// https://gallery.rcpp.org/articles/quiet-stop-and-warning/
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
    /// Ref:
    /// -----
    /// Suppressing Call Stack Info in Rcpp-Generated Errors and Warnings
    /// Michael Weylandt — written Mar 27, 2018
    /// https://gallery.rcpp.org/articles/quiet-stop-and-warning/
    void warning(const std::string &txt)
    {
        if (!is_new_housepi_warning(txt))
            return;
        Rf_warningcall(R_NilValue, "%s", txt.c_str());
        add_new_housepi_warning(txt);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  *** Housepi Rcpp Base : class ***  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using namespace hsp;

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

        Rcpp::NumericVector call_(Rcpp::NumericVector x) const
        {
            try
            {
                const R_xlen_t n = x.size();
                Rcpp::NumericVector y(n);
                for (R_xlen_t i = 0; i < n; ++i)
                    y[i] = ptr_->evaluate(x[i]);
                return y;
            }
            catch (const std::exception &e)
            {
                msg::error(
                    "While calling " + family_ + " '" + std::string(ptr_->name()) + "':\n" + e.what());
            }
        }
    };

    class ExposeFunction : public ExposeBase
    {
    public:
        ExposeFunction(std::string name_)
            : ExposeBase(name_, "function",
                         {"custom", "transcendental", "arithmetic"}) {}

        Rcpp::NumericVector func(Rcpp::NumericVector x) const
        {
            ptr_->set_eval_call();
            return call_(x);
        }
        Rcpp::NumericVector prim(Rcpp::NumericVector x) const
        {
            ptr_->set_eval_primitive();
            return call_(x);
        }
    };

    class ExposeDistribution : public ExposeBase
    {
    public:
        ExposeDistribution(std::string name_)
            : ExposeBase(name_, "distribution",
                         {"continuos", "discrete"}) {}

        Rcpp::NumericVector pdf(Rcpp::NumericVector x) const
        {
            ptr_->set_eval_as_distr_call();
            return call_(x);
        }
        Rcpp::NumericVector cdf(Rcpp::NumericVector x) const
        {
            ptr_->set_eval_as_distr_primitive();
            return call_(x);
        }
    };

    class ExposeInfecProfile : public ExposeBase
    {
    public:
        ExposeInfecProfile(std::string name_)
            : ExposeBase(name_, "infectivity profile",
                         {"infecprofile", "continuos", "discrete"}) {}

        Rcpp::NumericVector inst(Rcpp::NumericVector x) const
        {
            ptr_->set_eval_as_distr_call();
            return call_(x);
        }
        Rcpp::NumericVector cumul(Rcpp::NumericVector x) const
        {
            ptr_->set_eval_as_distr_primitive();
            return call_(x);
        }
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

    // Copies src into dest. If destroy_src is true, releases src's storage after copy.
    template <typename DestNumeric>
    void helper_vector_flush(DestNumeric dest, std::vector<double> &src)
    {
        std::copy(src.begin(), src.end(), dest.begin());
        std::vector<double>().swap(src);
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
    MAKE_RCPP_IMPORT_METHOD(param, "parameters")

    void import_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        try
        {
            mcmc_ptr->detach_data_store();
            mcmc_ptr->data_wide().build(rdataframe_to_strmap(input));
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
            mcmc_ptr->infec_profile.build(mcmc_ptr->param, rdataframe_to_strmap(input));
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
        try
        {
            if (data_wide.is_empty())
                return Rcpp::List();

            size_t n_rows = data_wide.report_indiv.total;
            Rcpp::CharacterVector id_col = Rcpp::CharacterVector(n_rows);
            Rcpp::CharacterVector house_id_col = Rcpp::CharacterVector(n_rows);
            Rcpp::IntegerVector status_col = Rcpp::IntegerVector(n_rows);

            // other non standard data columns
            std::vector<size_t> list_export = data_wide.pKeys->from_string(columns);
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
        catch (const std::exception &e)
        {
            msg::error("While exporting data (wide):\n" + std::string(e.what()));
        }
    }

    Rcpp::DataFrame export_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        return Rcpp::DataFrame(helper_export_data_wide(data_wide, data_wide.pKeys->get_all(), true));
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

        const std::vector<std::string> data_names = data_wide.pKeys->get_all() -
                                                    augm.pKeys->as_string(augm.list_export);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &member : house)
            {
                size_t c = 0;
                for (size_t a : augm.list_export)
                    augm_cols[c++][row] = hsp_to_rcpp(member.get_augm_variable(a));
                ++row;
            }

        Rcpp::List output = helper_export_data_wide(data_wide, data_names, true);
        size_t i = 0;
        for (size_t a : augm.list_export)
            output[augm.pKeys->as_string(a)] = augm_cols[i++];
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

        try
        {
            if (data_long.indiv_datalist.size() == 0)
                return Rcpp::DataFrame();

            // retrieve total number of rows
            size_t n_rows = 0;
            for (const auto &tbl : data_long.indiv_datalist)
                if (!tbl.vec.empty())
                    n_rows += tbl.vec.begin()->size();

            Rcpp::CharacterVector datalist_index(n_rows);

            // other non standard data-long columns
            std::vector<Rcpp::NumericVector> datalist_dbl(data_long.pKeys->size());
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
                output.push_back(datalist_dbl[d], data_long.pKeys->as_string(d));

            return Rcpp::DataFrame(output);
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting data (long):\n" + std::string(e.what()));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  MCMCBase : report exports  ///////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export current data configurations to R DataFrame.
    Rcpp::DataFrame summarise_cases(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const ReportIndividual &report_indiv = data_wide.report_indiv;
        const ReportHousehold &report_house = data_wide.report_house;

        try
        {
            if (data_wide.is_empty())
                msg::warning("Called 'summarise_cases' on empty data.");

            Rcpp::CharacterVector type = Rcpp::wrap(std::vector<std::string>(4, "individual") +
                                                    std::vector<std::string>(4, "household"));
            Rcpp::CharacterVector comment = {"total", "missing", "not infected", "infected",
                                             "total", "no cases", "at least 1 case", "at least 2 cases"};
            Rcpp::IntegerVector counts(type.size());
            Rcpp::NumericVector pct(type.size());

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

            return Rcpp::DataFrame::create(
                Rcpp::_["type"] = std::move(type),
                Rcpp::_["comment"] = std::move(comment),
                Rcpp::_["counts"] = std::move(counts),
                Rcpp::_["pct"] = std::move(pct));
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting data report:\n" +
                       std::string(e.what()));
        }
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

        try
        {
            Rcpp::CharacterVector comment({"max cases", "max miss", "min size", "max size"});
            R_xlen_t n_rows = comment.size();
            Rcpp::IntegerVector value(n_rows), counts(n_rows);
            Rcpp::NumericVector pct(n_rows);

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

            return Rcpp::DataFrame::create(
                Rcpp::_["comment"] = std::move(comment),
                Rcpp::_["value"] = std::move(value),
                Rcpp::_["counts"] = std::move(counts),
                Rcpp::_["pct"] = std::move(pct));
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting household size report:\n" +
                       std::string(e.what()));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMCBase : export draws  ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export MCMC draws of total log-likelihood posterior samples.
    Rcpp::NumericMatrix flush_total_log_lik_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        TotalLogLike &tot_log_lik = mcmc_ptr->tot_log_lik;

        try
        {
            Rcpp::NumericMatrix output(tot_log_lik.count(), 1);
            helper_vector_flush(output(Rcpp::_, 0), tot_log_lik.smpl_values);
            colnames(output) = Rcpp::CharacterVector({"log_lik"});
            return output;
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting total log-likelihood draws:\n" +
                       std::string(e.what()));
        }
    }

    /// Export MCMC draws of parameter values posterior samples.
    Rcpp::NumericMatrix flush_param_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Parameter &param = mcmc_ptr->param;
        const TotalLogLike &tot_log_lik = mcmc_ptr->tot_log_lik;

        if (param.list_sample.empty())
            return Rcpp::NumericMatrix();

        try
        {
            size_t n_cols = param.list_sample.size();
            Rcpp::CharacterVector headers = Rcpp::wrap(param.pKeys->as_string(param.list_sample));
            Rcpp::NumericMatrix output(tot_log_lik.count(), n_cols);

            for (size_t c = 0; c < n_cols; ++c)
                helper_vector_flush(output(Rcpp::_, c), param.smpl_values[c]);

            colnames(output) = std::move(headers);
            param.smpl_values.clear();
            return output;
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting parameters draws:\n" +
                       std::string(e.what()));
        }
    }

    /// Export MCMC draws of parameter log-priors posterior samples.
    Rcpp::NumericMatrix flush_param_lprior_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Parameter &param = mcmc_ptr->param;
        const TotalLogLike &tot_log_lik = mcmc_ptr->tot_log_lik;

        if (param.list_sample.empty())
            return Rcpp::NumericMatrix();

        try
        {
            size_t n_cols = param.list_sample.size();
            Rcpp::CharacterVector headers(n_cols);
            Rcpp::NumericMatrix output(tot_log_lik.count(), n_cols);

            size_t c = 0;
            for (size_t p : param.list_sample)
            {
                headers[c] = "lprior_" + param.pKeys->as_string(p);
                helper_vector_flush(output(Rcpp::_, c), param.smpl_priors[c]);
                ++c;
            }
            colnames(output) = std::move(headers);
            param.smpl_priors.clear();
            return output;
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting parameter log-priors draws:\n" +
                       std::string(e.what()));
        }
    }

    /// Export MCMC draws of augmented posterior samples.
    Rcpp::NumericMatrix flush_augm_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Augmented &augm = mcmc_ptr->augm;
        const TotalLogLike &tot_log_lik = mcmc_ptr->tot_log_lik;

        if (augm.smpl_values.empty())
            return Rcpp::NumericMatrix();

        try
        {
            size_t n_cols = 0;
            for (const Individual *member : augm.smpl_individuals)
                n_cols += augm.list_draws(*member).size();

            Rcpp::CharacterVector headers(n_cols);
            Rcpp::NumericMatrix output(tot_log_lik.count(), n_cols);

            size_t col = 0;
            for (size_t i = 0; i < augm.smpl_values.size(); ++i)
            {
                const Individual &member = *augm.smpl_individuals[i];
                auto &member_matrix = augm.smpl_values[i];

                size_t j = 0;
                for (size_t a : augm.list_draws(member))
                {
                    headers[col] = member.id() + "_" + augm.pKeys->as_string(a);
                    helper_vector_flush(output(Rcpp::_, col++), member_matrix[j++]);
                }
            }

            colnames(output) = headers;
            augm.smpl_values.clear();
            return output;
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting augmented draws:\n" +
                       std::string(e.what()));
        }
    }

    /// Export all sampled values of the log-likelihood for all household.
    Rcpp::NumericMatrix flush_house_log_lik_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        DataWide &data_wide = mcmc_ptr->data_wide();
        const TotalLogLike &tot_log_lik = mcmc_ptr->tot_log_lik;

        try
        {
            if (data_wide.smpl_values.empty())
                return Rcpp::NumericMatrix();

            // create a NumericMatrix with n_house columns
            size_t n_houses = data_wide.report_house.total;
            Rcpp::NumericMatrix output(tot_log_lik.count(), n_houses);
            Rcpp::CharacterVector headers(n_houses);

            // retrieve restricted index by general index of multiplicity head class
            size_t n_classes = data_wide.restricted_size();
            std::vector<std::vector<size_t>> house_classes(n_classes);

            for (size_t H = 0; H < n_houses; ++H)
            {
                size_t R = data_wide.get_restricted_index(H);
                house_classes[R].push_back(H);
                headers[H] = data_wide.at(H).house_id();
            }

            for (size_t R = 0; R < n_classes; ++R)
            {
                std::vector<double> &src = data_wide.smpl_values[R];

                for (size_t H : house_classes[R])
                    std::copy(src.begin(), src.end(), output(Rcpp::_, H).begin());
                std::vector<double>().swap(src);
            }
            colnames(output) = std::move(headers);

            return output;
        }
        catch (const std::exception &e)
        {
            msg::error("While exporting household log-likelihood draws:\n" +
                       std::string(e.what()));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMC : house state snapshot /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Export to R a snapshot of current state of house data.
    /// NOTE: there are not private parameters so far.
    Rcpp::NumericVector helper_param_rcpp_snapshot(const Parameter &param)
    {
        if (param.pKeys->empty())
            return Rcpp::NumericVector::create();

        size_t n = param.pKeys->size();
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
    Rcpp::DataFrame helper_data_wide_rcpp_snapshot(const DataWide &data_wide,
                                                   const Household &house)
    {
        const std::vector<std::string> &headers = data_wide.pKeys->get_all();

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
    Rcpp::DataFrame helper_data_long_rcpp_snapshot(const DataLong &data_long,
                                                   const Household &house)
    {
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
        std::vector<Rcpp::NumericVector> df_double(data_long.pKeys->size());
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
            output[data_long.pKeys->as_string(d)] = df_double[d];

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_group_rates_rcpp_snapshot(const GroupRate &group_rate,
                                                     const Household &house)
    {
        const std::vector<std::string> &headers = group_rate.pKeys->get_all();

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

    /// Export to R a snapshot of current state of house infectivity profile.
    Rcpp::List helper_infec_profile_rcpp_snapshot(const InfecProfile &infec_profile,
                                                  const Household &house)
    {
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
                const auto &f = house.member(i).get_augm_variable(infec_profile.augm_key).function();

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

    Rcpp::DataFrame helper_augm_log_lik_rcpp_snapshot(const Augmented &augm, const Household &house)
    {
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
                const Variable &xaugm = member_i.get_augm_variable(a);
                if (xaugm.has_model())
                    df_double[c++][i] = hsp_to_rcpp(xaugm.get_distr_log_pdf());
            }
        }

        Rcpp::List output;
        output["id"] = id_col;
        size_t c = 0;
        for (size_t a : augm.list_export)
            output.push_back(df_double[c++], augm.pKeys->as_string(a));

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_augm_rcpp_snapshot(const Augmented &augm, const Household &house)
    {
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
                df_double[c++][i] = hsp_to_rcpp(member_i.get_augm_variable(a));
        }

        Rcpp::List output;
        output["id"] = id_col;
        size_t i = 0;
        for (size_t a : augm.list_export)
            output.push_back(df_double[i++], augm.pKeys->as_string(a));

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
                Rcpp::_["obs_log_lik"] = house.obs_log_lik(),
                Rcpp::_["augm_log_lik"] = helper_augm_log_lik_rcpp_snapshot(mcmc_ptr->augm, house),
                Rcpp::_["param"] = helper_param_rcpp_snapshot(mcmc_ptr->param),
                Rcpp::_["data_wide"] = helper_data_wide_rcpp_snapshot(mcmc_ptr->data_wide(), house),
                Rcpp::_["data_long"] = helper_data_long_rcpp_snapshot(mcmc_ptr->data_long(), house),
                Rcpp::_["group_rate"] = helper_group_rates_rcpp_snapshot(mcmc_ptr->group_rate, house),
                Rcpp::_["infec_profile"] = helper_infec_profile_rcpp_snapshot(mcmc_ptr->infec_profile, house),
                Rcpp::_["augm"] = helper_augm_rcpp_snapshot(mcmc_ptr->augm, house));
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
        // pre-compute per-household row offsets (fixes OpenMP race on row counter)
        std::vector<size_t> get_row_offsets(const DataWide &data_wide)
        {
            std::vector<size_t> output(data_wide.report_house.total + 1, 0);
            for (size_t H = 0; H < data_wide.report_house.total; ++H)
                output[H + 1] = output[H] + data_wide.at(H).size();
            return output;
        }

        // pre-compute per-household row offsets (fixes OpenMP race on row counter)
        std::vector<std::string> get_col_names(const Augmented &augm,
                                               const std::vector<std::string> &static_data_names)
        {
            // pre-build the column name vector for the per-iteration shallow list
            const size_t n_shared = static_data_names.size();
            const size_t n_augm_cols = augm.list_export.size();

            std::vector<std::string> output(2 + n_shared + 1 + n_augm_cols);

            output[0] = "id";
            output[1] = "house_id";
            for (size_t i = 0; i < n_shared; ++i)
                output[2 + i] = static_data_names[i];

            output[2 + n_shared] = "status";
            for (size_t c = 0; c < n_augm_cols; ++c)
                output[3 + n_shared + c] = augm.pKeys->as_string(augm.list_export[c]);

            return output;
        }

        /// Perform a single repetition.
        void single_repetition(MCMC &engine,
                               Rcpp::IntegerVector &status_col,
                               std::vector<Rcpp::NumericVector> &augm_cols,
                               const std::vector<size_t> &row_offsets)
        {
            DataWide &data_wide = engine.data_wide();
            const Augmented &augm = engine.augm;

#pragma omp parallel for
            for (size_t H = 0; H < data_wide.report_house.total; ++H)
                try
                {
                    if (engine.atomic_error.is_pending())
                        continue;

                    Household &house = data_wide.at(H);
                    engine.model->simulate_single_house(house);

                    size_t row = row_offsets[H]; // thread-local: no shared counter
                    for (const Individual &member : house)
                    {
                        status_col[row] = member.status();
                        size_t c = 0;
                        for (size_t a : augm.list_export)
                            augm_cols[c++][row] = hsp_to_rcpp(member.get_augm_variable(a));
                        ++row;
                    }
                }
                catch (const std::exception &e)
                {
                    engine.atomic_error.set(
                        "While simulating " + data_wide.at(H).err_name() + ":\n" +
                        std::string(e.what()));
                }
            engine.atomic_error.throw_exception();
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
            //////////////////////////
            //////  Assertions  //////
            //////////////////////////

            data_wide.assert_non_empty();              // check if data is empty
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

            // column names of data shared across simulations
            const std::vector<std::string> static_data_names = data_wide.pKeys->get_all() -
                                                               augm.pKeys->as_string(augm.list_export);
            const Rcpp::List shared_df = helper_export_data_wide(data_wide, static_data_names, false);

            // simulated data column name vector
            Rcpp::CharacterVector sim_names = Rcpp::wrap(simulate::get_col_names(augm, static_data_names));

            // row offsets for house parallelization
            const std::vector<size_t> row_offsets = simulate::get_row_offsets(data_wide);

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
    ///
    /// Ref.
    /// -----
    /// Ref: Friel, N., Hurn, M. & Wyse, J.
    /// Improving power posterior estimation of statistical evidence.
    /// Stat Comput 24, 709–723 (2014).
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
    //////////////////////////////////  *** housepi - Rcpp exports ***  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_EXPORT_CHARACTER_VECTOR(method)                            \
    Rcpp::CharacterVector method(Rcpp::XPtr<hsp::ModelBase> model_xptr) \
    {                                                                   \
        return Rcpp::wrap(model_xptr->method());                        \
    }

    MAKE_EXPORT_CHARACTER_VECTOR(compile_flags)
    MAKE_EXPORT_CHARACTER_VECTOR(get_DataKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_LongKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_ParamKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_GroupRateKey)
    MAKE_EXPORT_CHARACTER_VECTOR(get_AugmKey)

    // UPDATE: I am just using this to check flags
    Rcpp::XPtr<ModelBase> new_model(std::string module_name)
    {
        auto model = hsp::ModelRegistry::instance().create(module_name);
        return Rcpp::XPtr<ModelBase>(model.release(), true);
    }

    Rcpp::XPtr<MCMC> create_mcmc_runner(std::string module_name)
    {
        auto mcmc = Rcpp::XPtr<MCMC>(new MCMC(module_name), true);

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

    double evaluate_log_lik(Rcpp::XPtr<MCMC> mcmc_ptr, std::string house_id)
    {
        try
        {
            return mcmc_ptr->evaluate_log_lik(house_id);
        }
        catch (const std::exception &e)
        {
            mcmc_ptr->setup.crash(e);
        }
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
        .method("prim", &ExposeFunction::prim);

    Rcpp::class_<ExposeDistribution>("ExposeDistribution")
        .constructor<std::string>()
        .derives<ExposeBase>("ExposeBase")
        .method("pdf", &ExposeDistribution::pdf)
        .method("cdf", &ExposeDistribution::cdf);

    Rcpp::class_<ExposeInfecProfile>("ExposeInfecProfile")
        .constructor<std::string>()
        .derives<ExposeBase>("ExposeBase")
        .method("inst", &ExposeInfecProfile::inst)
        .method("cumul", &ExposeInfecProfile::cumul);

    Rcpp::function("compile_flags", &compile_flags);
    Rcpp::function("get_DataKey", &get_DataKey);
    Rcpp::function("get_LongKey", &get_LongKey);
    Rcpp::function("get_ParamKey", &get_ParamKey);
    Rcpp::function("get_GroupRateKey", &get_GroupRateKey);
    Rcpp::function("get_AugmKey", &get_AugmKey);

    Rcpp::function("new_model", &new_model);
    Rcpp::function("create_mcmc_runner", &create_mcmc_runner);
    Rcpp::function("run_inference", &run_inference);
    Rcpp::function("run_simulation", &run_simulation);
    Rcpp::function("evaluate_log_lik", &evaluate_log_lik);
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
    Rcpp::function("summarise_cases", &summarise_cases);
    Rcpp::function("summarise_house_sizes", &summarise_house_sizes);
    Rcpp::function("export_setup", &export_setup);
    Rcpp::function("export_augmented_data_wide", &export_augmented_data_wide);
    Rcpp::function("export_infec_profile", &export_infec_profile);
    Rcpp::function("flush_total_log_lik_draws", &flush_total_log_lik_draws);
    Rcpp::function("flush_param_draws", &flush_param_draws);
    Rcpp::function("flush_param_lprior_draws", &flush_param_lprior_draws);
    Rcpp::function("flush_augm_draws", &flush_augm_draws);
    Rcpp::function("flush_house_log_lik_draws", &flush_house_log_lik_draws);
    Rcpp::function("export_house_snapshot", &export_house_snapshot);
}
#endif
