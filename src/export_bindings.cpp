/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef _STANDALONE
#include <functional>
#include <memory>
#include <exception>
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
#include "group_risk.h"
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
    [[noreturn]] void stop(const std::string &txt)
    {
        throw housepi::CppError(txt);
    }

    /// Print an information message to stdout.
    void info(const std::string &txt) noexcept
    {
        Rcpp::Rcout << txt;
        Rcpp::Rcout.flush();
    }

    /// Display a warning message.
    /// Raised as a classed condition to identify housepi-authored warnings.
    void warn(const std::string &txt)
    {
        Rcpp::List cond = Rcpp::List::create(Rcpp::_["message"] = txt,
                                             Rcpp::_["call"] = R_NilValue);
        cond.attr("class") = Rcpp::CharacterVector::create(
            "housepi_cpp_warning", "warning", "condition");
        Rcpp::Function("warning")(cond);
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

        ExposeBase(
            std::string name,
            ErrType family,
            std::vector<Registry> registry_types)
            : ptr_(compose_func_registries(name, family, registry_types, false)),
              family_(family) {}
        virtual ~ExposeBase() = default;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        void init(SEXP args_, double min, double max, bool is_discr)
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
        ErrType family_;

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
                msg::stop(e.what());
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
                        "Size of " + msg::arg("xo") + " must be identical to size of " + msg::arg("xn") +
                        " for " + err_type_name(family_) + " " + msg::val(ptr_->name()) + ".");
                Rcpp::NumericVector y(n);
                for (R_xlen_t i = 0; i < n; ++i)
                    y[i] = fn(xo[i], xn[i]);
                return y;
            }
            catch (const std::exception &e)
            {
                msg::stop(e.what());
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
        static ErrType family() { return ErrType::Function; }
        static std::vector<Registry> registries() { return {Registry::Custom, Registry::Transcendental, Registry::Arithmetic}; }

        ExposeFunction(std::string name_)
            : ExposeBase(name_, family(), registries()) {}

        MAKE_EXPOSE_METHOD_CALL(func, ptr_->func(v))
        MAKE_EXPOSE_METHOD_CALL(log_func, ptr_->log_func(v))
        MAKE_EXPOSE_METHOD_CALL_RATIO(log_ratio, ptr_->guarded_log_ratio(vo, vn))
        MAKE_EXPOSE_METHOD_CALL(prim, ptr_->prim(v))
    };

    class ExposeDistribution : public ExposeBase
    {
    public:
        static ErrType family() { return ErrType::Distribution; }
        static std::vector<Registry> registries() { return {Registry::Continuos, Registry::Discrete}; }

        ExposeDistribution(std::string name_)
            : ExposeBase(name_, family(), registries()) {}

        MAKE_EXPOSE_METHOD_CALL(pdf, ptr_->pdf(v))
        MAKE_EXPOSE_METHOD_CALL(log_pdf, distr_()->guarded_log_pdf(v))
        MAKE_EXPOSE_METHOD_CALL_RATIO(log_ratio, distr_()->guarded_log_ratio(vo, vn))
        MAKE_EXPOSE_METHOD_CALL(cdf, ptr_->cdf(v))
        MAKE_EXPOSE_METHOD_CALL(quantile, distr_()->quantile(v))

        /// Sample n random draws. Seed is drawn from R's RNG so that
        /// set.seed() in R controls the result.
        Rcpp::NumericVector rand(int n) const
        {
            try
            {
                if (n < 0)
                    throw std::invalid_argument(
                        msg::arg("n") + " must be non-negative for " +
                        err_type_name(family_) + " " + msg::val(ptr_->name()) + ".");

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
                msg::stop(e.what());
            }
        }

    private:
        Distribution *distr_() const { return static_cast<Distribution *>(ptr_.get()); }
    };

    class ExposeInfecProfile : public ExposeBase
    {
    public:
        static ErrType family() { return ErrType::InfecProfile; }
        static std::vector<Registry> registries() { return {Registry::InfecProfile, Registry::Continuos, Registry::Discrete}; }

        ExposeInfecProfile(std::string name_)
            : ExposeBase(name_, family(), registries())
        {
            // core errors from the underlying object should say "infectivity
            // profile X", not the Distribution default "distribution X".
            if (ptr_)
                ptr_->set_err_family(ErrType::InfecProfile);
        }

        MAKE_EXPOSE_METHOD_CALL(inst, ptr_->pdf(v))
        MAKE_EXPOSE_METHOD_CALL(cumul, ptr_->cdf(v))
    };

    // Factory free functions: construct + initialise in one call, so the whole
    // build routes through the `.engine` proxy interceptor like every other
    // module free function. Ownership of the returned object passes to R.
    template <class Expose>
    Expose *expose_(std::string name, SEXP args, double min, double max, bool discr)
    {
        auto obj = std::make_unique<Expose>(std::move(name));
        obj->init(args, min, max, discr);
        return obj.release();
    }

    // List every name the family accepts from the inherent registries.
    template <class Expose>
    std::vector<std::string> list_options_()
    {
        return registry_options(Expose::registries(), Expose::family());
    }

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
    /// Auxiliary converter of an individual status, missing status to R NA integer.
    int hsp_to_rcpp(const Individual &indiv)
    {
        return indiv.is_missing_status() ? NA_INTEGER : static_cast<int>(indiv.status());
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

    // Import helpers run serially at the Rcpp boundary; their build() paths
    // raise housepi::CppError (via msg::stop) directly, which Rcpp forwards to R
    // with the right class, so no per-import try/catch conversion is needed.
    void import_setup(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->setup.build(rdataframe_to_strmap(input));
    }

    void import_param(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->param.build(rdataframe_to_strmap(input), mcmc_ptr->translator);
    }

    void import_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->detach_data_store();
        mcmc_ptr->data_wide().build(rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_data_long(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        /// FIXME: copy-on-write should be also on long? Not clear...
        mcmc_ptr->data_long().build(mcmc_ptr->data_wide(), rdataframe_to_strmap(input));
    }

    void import_group_risk(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->group_risk.build(mcmc_ptr->data_wide(), mcmc_ptr->param, rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_contact_risk(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->contact_risk.build(mcmc_ptr->data_wide(), mcmc_ptr->param, rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->infec_profile.build(mcmc_ptr->translator, rdataframe_to_strmap(input));
    }

    void import_augm(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->augm.build(mcmc_ptr->infec_profile, rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_param_table(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->param.build_for_predic(rdataframe_to_dblmap(input));
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
            Rcpp::_["func"] = infec_profile.generic->name(),
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
            Rcpp::_["progress"] = setup.progress,
            Rcpp::_["inv_temp"] = setup.inv_temp,
            Rcpp::_["time_step"] = setup.time_step);
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
        std::vector<size_t> data_keys = data_wide.khold(columns);
        std::vector<Rcpp::NumericVector> data_col(columns.size());
        for (size_t d = 0; d < columns.size(); ++d)
            data_col[d] = Rcpp::NumericVector(n_rows);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &indiv : house)
            {
                id_col[row] = data_wide.id(indiv);
                house_id_col[row] = data_wide.house_id(house);
                status_col[row] = hsp_to_rcpp(indiv);

                for (size_t d = 0; d < columns.size(); ++d)
                    data_col[d][row] = hsp_to_rcpp(data_wide.data(data_keys[d], indiv));
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
        MCMC &engine = *mcmc_ptr;

        // a deterministic observation stands in for its data column: build it here so the
        // exported frame reports the value the model reads, run or no run. Finalizing
        // collapses each wired cell to a constant, so no read-cell outlives the Variable
        // it points at (the next definition import clears them) -- which is also why a
        // failed build rewinds the read-surface before the error travels on.
        if (!engine.augm.list_surrogate.empty())
        {
            engine.data_wide().reset_readcells(engine.augm.get_produced_cols());
            try
            {
                engine.augm.materialize_surrogate_columns(true);
                engine.data_wide().finalize();
            }
            catch (const std::exception &)
            {
                engine.data_wide().reset_readcells();
                throw;
            }
        }

        const DataWide &data_wide = engine.data_wide();
        return Rcpp::DataFrame(helper_export_data_wide(data_wide, data_wide.imported_cols(), true));
    }

    /// Export to R all augmented variables at the current state.
    Rcpp::DataFrame export_augmented_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const Augmented &augm = mcmc_ptr->augm;

        size_t n_cols = augm.list_augm_export.size();
        size_t n_rows = data_wide.report_indiv.total;

        std::vector<Rcpp::NumericVector> augm_cols(n_cols);
        for (size_t i = 0; i < n_cols; ++i)
            augm_cols[i] = Rcpp::NumericVector(n_rows);

        const std::vector<std::string> data_names = data_wide.khold().all() -
                                                    augm.khold(augm.list_augm_export);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &indiv : house)
            {
                size_t c = 0;
                for (size_t a : augm.list_augm_export)
                    augm_cols[c++][row] = hsp_to_rcpp(mcmc_ptr->get_augm_varia(a, indiv));
                ++row;
            }

        Rcpp::List output = helper_export_data_wide(data_wide, data_names, true);
        size_t i = 0;
        for (size_t a : augm.list_augm_export)
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
            for (const Individual &indiv : house)
            {
                const std::vector<std::vector<double>> *ptr = data_long.get_data_long_ptr(indiv);
                if (!ptr->empty())
                    for (size_t k = 0; k < ptr->begin()->size(); ++k)
                    {
                        datalist_index[row] = data_wide.id(indiv);
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
        count[2] = report_indiv.escaped;
        percent[2] = format_percentage(report_indiv.escaped, report_indiv.total);
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
        for (auto &indiv_matrix : post.augm_values)
        {
            for (auto &var_vec : indiv_matrix)
                cols[col++] = flush_vector_to_numvec(var_vec);
            std::vector<std::vector<double>>().swap(indiv_matrix);
        }
        std::vector<std::vector<std::vector<double>>>().swap(post.augm_values);

        cols.attr("names") = Rcpp::wrap(post.augm_col_names);
        helper_set_dataframe_attrs(cols, static_cast<R_xlen_t>(post.count()));
        return cols;
    }

    /// Rcpp wrapper: pointwise per-household log-likelihood draws as a data.frame.
    /// Each column is the complete-data household log-likelihood = transmission
    /// (+ observation, when an observation node is defined). Observation nodes force
    /// unique equivalence classes, so the obs store aligns 1:1 with the transmission
    /// store. Households of same equivalence class share NumericVector SEXP.
    /// R's reference counting keeps storage at n_classes * n_draws
    /// until a downstream consumer (e.g. as_draws_matrix) materialises it.
    Rcpp::List flush_house_log_lik_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        const DataWide &data_wide = mcmc_ptr->data_wide();

        size_t n_classes = post.house_log_lik.size();
        size_t n_houses = data_wide.report_house.total;

        // One NumericVector per class = transmission (+ observation when present),
        // already summed at store time.
        std::vector<Rcpp::NumericVector> class_cols(n_classes);
        for (size_t R = 0; R < n_classes; ++R)
            class_cols[R] = flush_vector_to_numvec(post.house_log_lik[R]);
        std::vector<std::vector<double>>().swap(post.house_log_lik);

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

        for (const Individual &indiv : house)
        {
            size_t i = Index::member(indiv);
            id_col[i] = data_wide.id(indiv);
            house_id_col[i] = data_wide.house_id(house);
            status_col[i] = hsp_to_rcpp(indiv);

            for (size_t d = 0; d < headers.size(); ++d)
                df_double[d][i] = hsp_to_rcpp(data_wide.data(d, indiv));
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
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const DataLong &data_long = mcmc_ptr->data_long();

        if (data_long.indiv_datalist.empty())
            return Rcpp::DataFrame::create();

        size_t n_rows = 0;
        for (const Individual &indiv : house)
        {
            const std::vector<std::vector<double>> *ptr = data_long.get_data_long_ptr(indiv);
            if (!ptr->empty())
                n_rows += ptr->begin()->size();
        }

        Rcpp::CharacterVector datalist_id(n_rows);
        std::vector<Rcpp::NumericVector> df_double(data_long.khold().size());
        for (size_t d = 0; d < df_double.size(); ++d)
            df_double[d] = Rcpp::NumericVector(n_rows);

        size_t r = 0;
        for (const Individual &indiv : house)
        {
            const std::vector<std::vector<double>> *ptr = data_long.get_data_long_ptr(indiv);
            if (!ptr->empty())
                for (size_t k = 0; k < ptr->begin()->size(); ++k)
                {
                    for (size_t d = 0; d < df_double.size(); ++d)
                        df_double[d][r] = hsp_to_rcpp(ptr->at(d)[k]);
                    datalist_id[r++] = data_wide.id(indiv);
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
    Rcpp::DataFrame helper_group_risks_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const GroupRisk &group_risk = mcmc_ptr->group_risk;
        const std::vector<std::string> &headers = group_risk.khold().all();

        if (headers.empty())
            return Rcpp::DataFrame::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(headers.size());
        for (size_t c = 0; c < headers.size(); ++c)
            df_double[c] = Rcpp::NumericVector(n_rows);

        for (const Individual &indiv : house)
        {
            size_t i = Index::member(indiv);
            id_col[i] = data_wide.id(indiv);

            for (size_t c = 0; c < headers.size(); ++c)
                df_double[c][i] = hsp_to_rcpp(group_risk.at(c, indiv));
        }

        Rcpp::List output;
        output["id"] = id_col;
        for (size_t c = 0; c < headers.size(); ++c)
            output.push_back(df_double[c], headers[c]);

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current contact risk attribute values
    /// for a household as a (from x to) numeric matrix with member id dimnames.
    /// Returns an empty 0x0 matrix when no contact risk was defined.
    Rcpp::NumericMatrix helper_contact_risks_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const ContactRisk &contact_risk = mcmc_ptr->contact_risk;

        if (contact_risk.is_empty())
            return Rcpp::NumericMatrix(0, 0);

        const size_t n = house.size();
        Rcpp::NumericMatrix mtx(n, n);
        Rcpp::CharacterVector ids(n);

        for (size_t i = 0; i < n; ++i)
        {
            ids[i] = data_wide.id(house.at(i));
            for (size_t j = 0; j < n; ++j)
                mtx(i, j) = hsp_to_rcpp(contact_risk.value(house, i, j));
        }

        mtx.attr("dimnames") = Rcpp::List::create(ids, ids);
        return mtx;
    }

    /// Export to R a snapshot of current state of house infectivity profile.
    Rcpp::List helper_infec_profile_rcpp_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                  const Household &house)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
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
            id_col[i] = data_wide.id(house.at(i));

        for (MemberIndex i : house.infected_members())
            if (infec_profile.is_augmented())
            {
                const auto &f = mcmc_ptr->get_augm_varia(infec_profile.augm_key, house.at(i)).function();

                infprof_name[i] = f->name();
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

    Rcpp::DataFrame helper_augm_log_lik_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const Augmented &augm = mcmc_ptr->augm;

        size_t n_cols = augm.list_augm_export.size();

        if (n_cols == 0)
            return Rcpp::DataFrame::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(n_cols);
        for (size_t c = 0; c < n_cols; ++c)
            df_double[c] = Rcpp::NumericVector(n_rows, NA_REAL);

        for (size_t i = 0; i < house.size(); ++i)
        {
            const Individual &indiv_i = house.at(i);
            id_col[i] = data_wide.id(indiv_i);

            size_t c = 0;
            for (size_t a : augm.list_augm_export)
            {
                const Variable &avar = mcmc_ptr->get_augm_varia(a, indiv_i);
                if (avar.has_distr())
                    df_double[c++][i] = hsp_to_rcpp(avar.log_pdf());
            }
        }

        Rcpp::List output;
        output["id"] = id_col;
        size_t c = 0;
        for (size_t a : augm.list_augm_export)
            output.push_back(df_double[c++], augm.khold(a));

        return Rcpp::DataFrame(output);
    }

    /// Export to R a snapshot of current state of house augmented variables.
    Rcpp::DataFrame helper_augm_rcpp_snapshot(
        Rcpp::XPtr<MCMC> mcmc_ptr, const Household &house)
    {
        const DataWide &data_wide = mcmc_ptr->data_wide();
        const Augmented &augm = mcmc_ptr->augm;
        size_t n_cols = augm.list_augm_export.size();

        if (n_cols == 0)
            return Rcpp::DataFrame::create();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(n_cols);
        for (size_t c = 0; c < n_cols; ++c)
            df_double[c] = Rcpp::NumericVector(n_rows);

        for (size_t i = 0; i < house.size(); ++i)
        {
            const Individual &indiv_i = house.at(i);
            id_col[i] = data_wide.id(indiv_i);
            size_t c = 0;
            for (size_t a : augm.list_augm_export)
                df_double[c++][i] = hsp_to_rcpp(mcmc_ptr->get_augm_varia(a, indiv_i));
        }

        Rcpp::List output;
        output["id"] = id_col;
        size_t i = 0;
        for (size_t a : augm.list_augm_export)
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
                Rcpp::_["param"] = helper_param_rcpp_snapshot(mcmc_ptr),
                Rcpp::_["data_wide"] = helper_data_wide_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["data_long"] = helper_data_long_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["group_risk"] = helper_group_risks_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["contact_risk"] = helper_contact_risks_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["infec_profile"] = helper_infec_profile_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["augm"] = helper_augm_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["augm_log_lik"] = helper_augm_log_lik_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["transm_log_lik"] = mcmc_ptr->data_wide().transm_log_lik(house).value());
        }
        catch (const std::exception &e)
        {
            msg::stop(std::string(e.what()) +
                      "\nFor household " + msg::val(house_id) + ".");
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
                        for (const Individual &indiv : house)
                        {
                            status_col[row] = hsp_to_rcpp(indiv);
                            size_t c = 0;
                            for (size_t a : augm.list_augm_export)
                                augm_cols[c++][row] = hsp_to_rcpp(engine.get_augm_varia(a, indiv));
                            // observation nodes: generated for infected members, NA otherwise
                            for (size_t a : augm.list_obs_member)
                                augm_cols[c++][row] = indiv.is_infected_status()
                                                          ? hsp_to_rcpp(engine.get_augm_varia(a, indiv))
                                                          : NA_REAL;
                            ++row;
                        }

                        // surrogate columns: a household reduction is one value shared by
                        // every member row, a member column is read on each member --
                        // including those never infected, which the model reads too
                        size_t c = augm.list_augm_export.size() + augm.list_obs_member.size();
                        for (const Augmented::Surrogate &surr : augm.list_surrogate)
                        {
                            if (!Augmented::is_own_column(surr))
                                continue;
                            if (surr.is_house)
                            {
                                const double v = hsp_to_rcpp(augm.house_varia[H][surr.slot].value());
                                for (size_t r = row_offsets[H]; r < row_offsets[H + 1]; ++r)
                                    augm_cols[c][r] = v;
                            }
                            else
                            {
                                const size_t a = augm.khold(augm.surrogate_name(surr));
                                size_t r = row_offsets[H];
                                for (const Individual &indiv : house)
                                    augm_cols[c][r++] = hsp_to_rcpp(engine.get_augm_varia(a, indiv));
                            }
                            ++c;
                        }
                    }
                    catch (const std::exception &e)
                    {
                        throw std::runtime_error(
                            std::string(e.what()) +
                            "\nFor " + data_wide.err_name(data_wide.at(H)) + ".");
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
            engine.augm.init_house_for_simulation();   // build household-level observations

            // FIXME: augmented variables are not built yet here
            // -> impossible to assert profile whe Mode::Augmented
            // engine.assert_valid_infectivity_profile(); // assert profile is valid

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

            // pre-allocate mutable output columns once (reused across iterations):
            // exported augmented, observation and surrogate columns
            Rcpp::IntegerVector status_col(data_wide.report_indiv.total);
            std::vector<Rcpp::NumericVector> augm_cols(
                augm.list_augm_export.size() +
                augm.list_obs_member.size() +
                augm.n_own_columns());
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

            // restore observed status and point every read-cell back at its imported
            // value, so post-simulation reads of data_wide return the input: the
            // generated values were handed back in the simulated frames, and the cells
            // wired to produce them refer to Variables this run is about to clear.
            engine.augm.clear_all();
            engine.data_wide().recover_all_status();
            engine.data_wide().reset_readcells();

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

    /// Return power-posterior thermodynamic ladder.
    /// Based on @cite friel_improving_2014
    Rcpp::List run_power_posteriors(Rcpp::XPtr<MCMC> mcmc_ptr, const std::string &scheme, unsigned long n_rungs)
    {
        MCMC &engine = *mcmc_ptr;

        try
        {
            const PowerLadder ladder = engine.run_power_posteriors(scheme, n_rungs);

            const size_t n_rung = ladder.draws.size();
            const size_t n_draw = n_rung ? ladder.draws[0].size() : 0;
            Rcpp::List draws(n_rung);
            Rcpp::CharacterVector col_names(n_rung);
            for (size_t r = 0; r < n_rung; ++r)
            {
                Rcpp::NumericVector col(n_draw);
                for (size_t c = 0; c < n_draw; ++c)
                    col[c] = ladder.draws[r][c];
                draws[r] = col;
                col_names[r] = "log_lik[" + std::to_string(r + 1) + "]";
            }
            // set names via attribute so the non-syntactic `log_lik[r]` survives
            // the data.frame conversion (which would otherwise sanitise them)
            Rcpp::DataFrame draws_df(draws);
            draws_df.attr("names") = col_names;

            return Rcpp::List::create(
                Rcpp::_["scheme"] = ladder.scheme,
                Rcpp::_["inv_temp"] = Rcpp::wrap(ladder.inv_temp),
                Rcpp::_["mean_log_lik"] = Rcpp::wrap(ladder.mean_log_lik),
                Rcpp::_["var_log_lik"] = Rcpp::wrap(ladder.var_log_lik),
                Rcpp::_["log_lik_draws"] = draws_df);
        }
        catch (std::exception &e)
        {
            engine.setup.crash(e);
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////  *** Rcpp exports ***  ////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

/// Keys the model declares through HOUSEPI_DECLARE_REGISTER_KEYS.
#define MAKE_EXPORT_CHARACTER_VECTOR(method)                         \
    Rcpp::CharacterVector method(Rcpp::XPtr<housepi::MCMC> mcmc_ptr) \
    {                                                                \
        return Rcpp::wrap(mcmc_ptr->model->method());                \
    }

    MAKE_EXPORT_CHARACTER_VECTOR(key_Data)
    MAKE_EXPORT_CHARACTER_VECTOR(key_Long)
    MAKE_EXPORT_CHARACTER_VECTOR(key_Param)
    MAKE_EXPORT_CHARACTER_VECTOR(key_GroupRisk)
    MAKE_EXPORT_CHARACTER_VECTOR(key_Augm)

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

    bool has_contact_risks(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return mcmc_ptr->model->has_contact_risks();
    }

    bool has_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return mcmc_ptr->model->has_infec_profile();
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

    /// Aggregate individuals by house_id. A missing status reads NA (NA_INTEGER is INT_MIN).
    measure::cases::Report summarise_houses_from_df(Rcpp::DataFrame data_wide)
    {
        Rcpp::CharacterVector col_names = data_wide.names();

        std::vector<std::string> house_id = {};
        if (contains(col_names, "house_id"))
            house_id = Rcpp::as<std::vector<std::string>>(data_wide["house_id"]);

        std::vector<int> status = {};
        if (contains(col_names, "status"))
            status = Rcpp::as<std::vector<int>>(data_wide["status"]);

        return measure::summarise_houses(house_id, status);
    }

    Rcpp::DataFrame count_houses_by_hhsize(Rcpp::DataFrame data_wide)
    {
        measure::cases::HouseCounts res = measure::house_counts(summarise_houses_from_df(data_wide));
        unsigned grid_len = static_cast<unsigned>(res.total.size());

        Rcpp::IntegerVector out_size(grid_len);
        Rcpp::IntegerVector out_total(grid_len);
        Rcpp::IntegerVector out_affected(grid_len);
        for (unsigned i = 0; i < grid_len; ++i)
        {
            out_size[i] = res.grid_min + i;
            out_total[i] = res.total[i];
            out_affected[i] = res.at_least_1_case[i];
        }
        return Rcpp::DataFrame::create(
            Rcpp::_["hhsize"] = out_size,
            Rcpp::_["total"] = out_total,
            Rcpp::_["at_least_1_case"] = out_affected);
    }

    Rcpp::DataFrame measure_ar_by_hhsize(Rcpp::DataFrame data_wide)
    {
        measure::AttackRates res = measure::attack_rate_by_hhsize(summarise_houses_from_df(data_wide));
        unsigned grid_len = static_cast<unsigned>(res.ar_mean.size());

        Rcpp::IntegerVector out_size(grid_len);
        Rcpp::NumericVector out_ar(grid_len);
        for (unsigned i = 0; i < grid_len; ++i)
        {
            out_size[i] = res.grid_min + i;
            out_ar[i] = res.count[i] > 0 ? res.ar_mean[i] : NA_REAL;
        }
        return Rcpp::DataFrame::create(
            Rcpp::_["hhsize"] = out_size,
            Rcpp::_["ar"] = out_ar);
    }

    Rcpp::DataFrame measure_sar_by_hhsize(Rcpp::DataFrame data_wide)
    {
        measure::SecAttackRates res = measure::sar_by_hhsize(summarise_houses_from_df(data_wide));
        unsigned grid_len = static_cast<unsigned>(res.sar_mean.size());

        Rcpp::IntegerVector out_size(grid_len);
        Rcpp::NumericVector out_sar(grid_len);
        for (unsigned i = 0; i < grid_len; ++i)
        {
            out_size[i] = res.grid_min + i;
            out_sar[i] = res.count[i] > 0 ? res.sar_mean[i] : NA_REAL;
        }
        return Rcpp::DataFrame::create(
            Rcpp::_["hhsize"] = out_size,
            Rcpp::_["sar"] = out_sar);
    }

    Rcpp::DataFrame measure_cases_by_hhsize(Rcpp::DataFrame data_wide)
    {
        measure::cases::Table res = measure::cases_by_hhsize(summarise_houses_from_df(data_wide));
        unsigned n_rows = static_cast<unsigned>(res.counts.size());

        Rcpp::IntegerVector out_size(n_rows);
        for (unsigned i = 0; i < n_rows; ++i)
            out_size[i] = res.min_size + i;

        Rcpp::List result;
        result.push_back(out_size, "hhsize");

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
}

RCPP_MODULE(housepi)
{
    Rcpp::class_<ExposeBase>("ExposeBase")
        .method("args", &ExposeBase::args);

    Rcpp::class_<ExposeFunction>("ExposeFunction")
        .derives<ExposeBase>("ExposeBase")
        .method("func", &ExposeFunction::func)
        .method("log_func", &ExposeFunction::log_func)
        .method("log_ratio", &ExposeFunction::log_ratio)
        .method("prim", &ExposeFunction::prim);

    Rcpp::class_<ExposeDistribution>("ExposeDistribution")
        .derives<ExposeBase>("ExposeBase")
        .method("pdf", &ExposeDistribution::pdf)
        .method("log_pdf", &ExposeDistribution::log_pdf)
        .method("log_ratio", &ExposeDistribution::log_ratio)
        .method("cdf", &ExposeDistribution::cdf)
        .method("quantile", &ExposeDistribution::quantile)
        .method("rand", &ExposeDistribution::rand);

    Rcpp::class_<ExposeInfecProfile>("ExposeInfecProfile")
        .derives<ExposeBase>("ExposeBase")
        .method("inst", &ExposeInfecProfile::inst)
        .method("cumul", &ExposeInfecProfile::cumul);

    Rcpp::function("expose_function", &expose_<ExposeFunction>);
    Rcpp::function("expose_distribution", &expose_<ExposeDistribution>);
    Rcpp::function("expose_infec_profile", &expose_<ExposeInfecProfile>);

    Rcpp::function("list_functions", &list_options_<ExposeFunction>);
    Rcpp::function("list_distributions", &list_options_<ExposeDistribution>);
    Rcpp::function("list_infec_profiles", &list_options_<ExposeInfecProfile>);

    Rcpp::function("key_Data", &key_Data);
    Rcpp::function("key_Long", &key_Long);
    Rcpp::function("key_Param", &key_Param);
    Rcpp::function("key_GroupRisk", &key_GroupRisk);
    Rcpp::function("key_Augm", &key_Augm);

    Rcpp::function("measure_sar_by_hhsize", &measure_sar_by_hhsize);
    Rcpp::function("measure_ar_by_hhsize", &measure_ar_by_hhsize);
    Rcpp::function("measure_cases_by_hhsize", &measure_cases_by_hhsize);
    Rcpp::function("count_houses_by_hhsize", &count_houses_by_hhsize);

    Rcpp::function("aux_chain_binomial_debugger", &housepi::chain_binom::aux_debugger);

    Rcpp::function("create_mcmc_runner", &create_mcmc_runner);
    Rcpp::function("run_inference", &run_inference);
    Rcpp::function("run_simulation", &run_simulation);
    Rcpp::function("run_power_posteriors", &run_power_posteriors);
    Rcpp::function("evaluate_log_lik", &evaluate_log_lik);
    Rcpp::function("import_data_wide", &import_data_wide);
    Rcpp::function("import_data_long", &import_data_long);
    Rcpp::function("share_data_store", &share_data_store);
    Rcpp::function("import_group_risk", &import_group_risk);
    Rcpp::function("import_contact_risk", &import_contact_risk);
    Rcpp::function("import_setup", &import_setup);
    Rcpp::function("import_param", &import_param);
    Rcpp::function("import_infec_profile", &import_infec_profile);
    Rcpp::function("import_augm", &import_augm);
    Rcpp::function("import_param_table", &import_param_table);
    Rcpp::function("has_contact_risks", &has_contact_risks);
    Rcpp::function("has_infec_profile", &has_infec_profile);
    Rcpp::function("export_data_wide", &export_data_wide);
    Rcpp::function("export_data_long", &export_data_long);
    Rcpp::function("summarise_house_cases", &summarise_house_cases);
    Rcpp::function("summarise_house_sizes", &summarise_house_sizes);
    Rcpp::function("export_setup", &export_setup);
    Rcpp::function("export_augmented_data_wide", &export_augmented_data_wide);
    Rcpp::function("flush_param_draws", &flush_param_draws);
    Rcpp::function("flush_param_log_prior_draws", &flush_param_log_prior_draws);
    Rcpp::function("flush_augm_draws", &flush_augm_draws);
    Rcpp::function("flush_house_log_lik_draws", &flush_house_log_lik_draws);
    Rcpp::function("export_house_snapshot", &export_house_snapshot);
}
#endif
