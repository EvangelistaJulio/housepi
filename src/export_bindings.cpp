/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef _STANDALONE
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <exception>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Rcpp.h>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "mcmc.h"
#include "measures.h"
#include "mcmc/posterior.h"
#include "registry.h"
#include "augmented.h"
#include "chain_binom.h"
#include "data.h"
#include "data/wide.h"
#include "data/long.h"
#include "distribution.h"
#include "function.h"
#include "risks.h"
#include "infec_profile.h"
#include "parameter.h"
#include "variable.h"

//////////////////////////////  msg  //////////////////////////////

namespace housepi
{
    /// Exception raised by housepi code.
    /// @note Keep its qualified name in sync with the R side.
    struct CppError : std::exception
    {
        std::string m;
        explicit CppError(std::string s) : m(std::move(s)) {}
        const char *what() const noexcept override { return m.c_str(); }
    };
}

namespace housepi::msg
{
    /// Throw `txt` as a CppError.
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

    /// Raise warning `txt` as a housepi_cpp_warning condition.
    void warn(const std::string &txt)
    {
        Rcpp::List cond = Rcpp::List::create(Rcpp::_["message"] = txt,
                                             Rcpp::_["call"] = R_NilValue);
        cond.attr("class") = Rcpp::CharacterVector::create(
            "housepi_cpp_warning", "warning", "condition");
        Rcpp::Function("warning")(cond);
    }
}

//////////////////////////////  Rcpp helpers  //////////////////////////////

namespace
{
    using namespace housepi;

    /// Convert a value to an R double, NA_REAL when missing.
    double hsp_to_rcpp(double x) { return na::is(x) ? NA_REAL : x; }
    double hsp_to_rcpp(const Attribute &attr)
    {
        if (attr.is_empty())
            return NA_REAL;
        if (attr.is_function() && !attr.function()->is_callable_implicitly())
            return NA_REAL;
        return hsp_to_rcpp(attr.value());
    }
    /// Convert the status of `indiv` to an R integer, NA_INTEGER when missing.
    int hsp_to_rcpp(const Individual &indiv)
    {
        return indiv.is_missing_status() ? NA_INTEGER : static_cast<int>(indiv.status());
    }

    /// Convert `input` to an IntegerVector when `discr` is requested.
    /// Values are assumed whole.
    /// - value beyond the int range: `input` returned unchanged.
    SEXP rcpp_typed_column(const Rcpp::NumericVector &input, bool discr)
    {
        static int machine_max_int = static_cast<double>(std::numeric_limits<int>::max());

        if (!discr)
            return input;
        Rcpp::IntegerVector out(input.size());
        for (R_xlen_t i = 0; i < input.size(); ++i)
        {
            if (ISNAN(input[i]))
                out[i] = NA_INTEGER;
            else if (std::abs(input[i]) > machine_max_int)
                return input;
            else
                out[i] = static_cast<int>(input[i]);
        }
        return out;
    }

    /// Convert an interval to an R list.
    /// - unbounded side: -Inf or Inf.
    Rcpp::List rcpp_interval(double min, double max,
                             Interval below, Interval above)
    {
        return Rcpp::List::create(
            Rcpp::_["min"] = below == Interval::Unbounded ? -inf::DOUBLE : hsp_to_rcpp(min),
            Rcpp::_["max"] = above == Interval::Unbounded ? inf::DOUBLE : hsp_to_rcpp(max),
            Rcpp::_["closed"] = Rcpp::LogicalVector::create(
                Rcpp::_["min"] = below == Interval::Closed,
                Rcpp::_["max"] = above == Interval::Closed));
    }

    /// Stamp a list of equally-sized columns as an R data.frame.
    /// Column names pass verbatim: Rcpp's data.frame conversion would sanitise them.
    /// They are data keys and definition names.
    Rcpp::DataFrame rcpp_data_frame(Rcpp::List output)
    {
        R_xlen_t n_rows = output.size() ? Rf_xlength(output[0]) : 0;
        // compact row names, the form R stores for the 1..n sequence
        output.attr("row.names") = Rcpp::IntegerVector::create(
            NA_INTEGER, -static_cast<int>(n_rows));
        output.attr("class") = "data.frame";
        return output;
    }

    /// Base wrapper exposing a housepi Function to R.
    class ExposeBase
    {
    public:
        //////  Constructor  //////

        ExposeBase(std::string name,
                   ErrType family,
                   std::vector<registry::Kind> reg_kinds)
            : ptr_(registry::create_function(name, family, reg_kinds, false)),
              family_(family) {}
        virtual ~ExposeBase() = default;

        //////  Setter  //////

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

            ptr_->reshape(is_discr);
        }

        //////  Getter  //////

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

        /// Whether the definition takes discrete values.
        /// - plain function: the values it returns;
        /// - distribution: the values of its support, so drawing families override this.
        virtual bool is_discrete() const { return ptr_->codomain().on_integers(); }
        /// Whether the definition is provably strictly positive.
        virtual bool is_positive() const { return ptr_->is_positive(); }
        /// Values returned.
        Rcpp::List codomain() const
        {
            const Codomain c = ptr_->codomain();
            return rcpp_interval(c.inf, c.sup, c.below, c.above);
        }
        /// Values it may be evaluated on, in the frame relative to the location.
        Rcpp::List domain() const
        {
            const Domain &d = ptr_->domain();
            return rcpp_interval(d.min.value(), d.max.value(), d.below, d.above);
        }
        /// Whether each bound reads an argument rather than holding its own value.
        Rcpp::LogicalVector bound_aliases() const
        {
            return Rcpp::LogicalVector::create(
                Rcpp::_["min"] = ptr_->min().is_pointer(),
                Rcpp::_["max"] = ptr_->max().is_pointer());
        }
        /// Location shift, NA without a location argument.
        double loc() const { return ptr_->loc().value(); }

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
                        " for " + str::format(family_) + " " + msg::val(ptr_->name()) + ".");
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
        static std::vector<registry::Kind> reg_kinds()
        {
            return {registry::Kind::Custom,
                    registry::Kind::Transcendental,
                    registry::Kind::Arithmetic};
        }

        ExposeFunction(std::string name_, bool /*discr*/)
            : ExposeBase(name_, family(), reg_kinds()) {}

        MAKE_EXPOSE_METHOD_CALL(func, ptr_->func(v))
        MAKE_EXPOSE_METHOD_CALL(log_func, ptr_->log_func(v))
        MAKE_EXPOSE_METHOD_CALL_RATIO(log_ratio, ptr_->guarded_log_ratio(vo, vn))
        MAKE_EXPOSE_METHOD_CALL(prim, ptr_->prim(v))
    };

    class ExposeDistribution : public ExposeBase
    {
    public:
        static ErrType family() { return ErrType::Distribution; }
        static std::vector<registry::Kind> reg_kinds(bool discr = false)
        {
            return registry::distr_kinds(discr);
        }

        ExposeDistribution(std::string name_, bool discr)
            : ExposeBase(name_, family(), reg_kinds(discr)) {}

        bool is_discrete() const override final { return distr_()->support().on_integers(); }
        bool is_positive() const override final { return distr_()->support().is_positive(); }

        MAKE_EXPOSE_METHOD_CALL(pdf, ptr_->pdf(v))
        MAKE_EXPOSE_METHOD_CALL(log_pdf, distr_()->guarded_log_pdf(v))
        MAKE_EXPOSE_METHOD_CALL_RATIO(log_ratio, distr_()->guarded_log_ratio(vo, vn))
        MAKE_EXPOSE_METHOD_CALL(cdf, ptr_->cdf(v))
        MAKE_EXPOSE_METHOD_CALL(quantile, distr_()->quantile(v))

        /// Draw `n` random values, integer-typed on a discrete distribution.
        /// The seed is drawn from R, so set.seed() controls the result.
        SEXP rand(int n) const
        {
            try
            {
                if (n < 0)
                    throw std::invalid_argument(
                        msg::arg("n") + " must be non-negative for " +
                        str::format(family_) + " " + msg::val(ptr_->name()) + ".");

                // generate seed from R
                Rcpp::RNGScope scope;
                std::uint32_t seed = static_cast<std::uint32_t>(
                    R::unif_rand() * std::numeric_limits<std::uint32_t>::max());

                rnd_gen prng(seed);
                Rcpp::NumericVector y(n);
                for (int i = 0; i < n; ++i)
                    y[i] = distr_()->rand(prng);
                return rcpp_typed_column(y, distr_()->support().on_integers());
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
        static std::vector<registry::Kind> reg_kinds(bool discr = false)
        {
            std::vector<registry::Kind> out = {registry::Kind::InfecProfile};
            vect::append(out, registry::distr_kinds(discr));
            return out;
        }

        ExposeInfecProfile(std::string name_, bool discr)
            : ExposeBase(name_, family(), reg_kinds(discr))
        {
            // core errors say "infectivity profile X", not the Distribution default "distribution X"
            if (ptr_)
                ptr_->set_err_family(ErrType::InfecProfile);
        }

        /// A profile is evaluated, never drawn: its discreteness is the lattice it is called on.
        /// reshape() puts that lattice on the integers.
        bool is_discrete() const override final { return ptr_->is_discretized(); }

        MAKE_EXPOSE_METHOD_CALL(inst, ptr_->pdf(v))
        MAKE_EXPOSE_METHOD_CALL(cumul, ptr_->cdf(v))
    };

    /// Construct and initialize an Expose object in one call.
    /// The whole build routes through the `.engine` proxy interceptor, like every module free function.
    /// Ownership of the returned object passes to R.
    template <class Expose>
    Expose *expose_(std::string name, SEXP args, double min, double max, bool discr)
    {
        auto obj = std::make_unique<Expose>(std::move(name), discr);
        obj->init(args, min, max, discr);
        return obj.release();
    }

    /// List every name the family accepts, from the registry kinds it looks up.
    template <class Expose>
    std::vector<std::string> list_options_()
    {
        return registry::list_names(Expose::reg_kinds());
    }

    //////////////////////////////  MCMC : imports from Rcpp  //////////////////////////////

    /// Convert an R data frame to a map of string vectors.
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

    /// Convert an R data frame to a map of double vectors.
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

    /// Export the argument list of `func` to R.
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

    /// Format `value` as a percentage of `total`.
    double format_percentage(size_t value, size_t total)
    {
        return std::round(10000. * static_cast<double>(value) / static_cast<double>(total)) / 100.;
    }

    // import helpers run serially at the Rcpp boundary
    // their build() paths raise housepi::CppError through msg::stop
    // Rcpp forwards it to R with its class
    // no import needs its own try/catch conversion
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
        mcmc_ptr->data.build_wide(rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_data_long(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        // @todo Decide whether long data should be copy-on-write too.
        mcmc_ptr->data.build_long(rdataframe_to_strmap(input));
    }

    void import_group_risk(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->risks.build_group(mcmc_ptr->data.wide(), mcmc_ptr->param, rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_contact_risk(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->risks.build_contact(mcmc_ptr->data.wide(), mcmc_ptr->param, rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->infec_profile.build(mcmc_ptr->translator, rdataframe_to_strmap(input));
    }

    void import_latent(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->augm.build(mcmc_ptr->infec_profile, rdataframe_to_strmap(input));
        mcmc_ptr->init_equivalence_class_compression();
    }

    void import_param_table(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
    {
        mcmc_ptr->param.build_for_predic(rdataframe_to_dblmap(input));
    }

    //////////////////////////////  MCMC : exports to Rcpp  //////////////////////////////

    /// Export the current setup to an R list.
    Rcpp::List export_setup(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const Setup &setup = mcmc_ptr->setup;

        return Rcpp::List::create(
            Rcpp::_["steps"] = setup.n_steps,
            Rcpp::_["burn"] = setup.n_burn,
            Rcpp::_["thin"] = setup.n_thin,
            Rcpp::_["adapt"] = setup.n_adapt,
            Rcpp::_["adapt_window"] = setup.n_adapt_window,
            Rcpp::_["chains"] = setup.n_chains,
            Rcpp::_["threads"] = setup.n_threads,
            Rcpp::_["id_draws"] = Rcpp::CharacterVector(Rcpp::wrap(setup.id_draws)),
            Rcpp::_["seed"] = setup.seed,
            Rcpp::_["progress"] = setup.progress,
            Rcpp::_["time_step"] = setup.time_step);
    }

    //////  MCMC : data::Wide : exports  //////

    Rcpp::List helper_export_data_wide(const Data &data,
                                       const std::vector<std::string> &columns,
                                       bool is_status)
    {
        const data::Wide &data_wide = data.wide();
        size_t n_rows = data_wide.report_indiv.total;
        Rcpp::CharacterVector id_col = Rcpp::CharacterVector(n_rows);
        Rcpp::CharacterVector hhid_col = Rcpp::CharacterVector(n_rows);
        Rcpp::IntegerVector status_col = Rcpp::IntegerVector(n_rows);

        // other non-standard data columns
        std::vector<size_t> data_keys = data_wide.khold(columns);
        std::vector<Rcpp::NumericVector> data_col(columns.size());
        for (size_t d = 0; d < columns.size(); ++d)
            data_col[d] = Rcpp::NumericVector(n_rows);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &indiv : house)
            {
                id_col[row] = data_wide.id(indiv);
                hhid_col[row] = data_wide.hhid(house);
                status_col[row] = hsp_to_rcpp(indiv);

                for (size_t d = 0; d < columns.size(); ++d)
                    data_col[d][row] = hsp_to_rcpp(data.at(data_keys[d], indiv).value());
                ++row;
            }

        Rcpp::List output;
        output["id"] = id_col;
        output["hhid"] = hhid_col;
        if (is_status)
            output["status"] = status_col;
        for (size_t d = 0; d < columns.size(); ++d)
            output[columns[d]] = data_col[d];

        return output;
    }

    /// Export the imported data columns to an R data frame.
    Rcpp::DataFrame export_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        MCMC &engine = *mcmc_ptr;

        // a deterministic observation stands in for its data column: build it here
        // the exported frame then reports the value the model reads, with or without a run
        // finalizing collapses each wired cell to a constant: no read-cell outlives its Variable
        // the next definition import clears the Variables
        // a failed build rewinds the read-surface before rethrowing
        if (!engine.augm.scheme().list_surrogate.empty())
        {
            engine.data.reset_readcells(engine.augm.get_produced_cols());
            try
            {
                engine.augm.materialize_surrogate_columns(true);
                engine.data.finalize();
            }
            catch (const std::exception &)
            {
                engine.data.reset_readcells();
                throw;
            }
        }

        return rcpp_data_frame(helper_export_data_wide(engine.data, engine.data.wide().imported_cols(), true));
    }

    /// Export data and latents at the current state to an R data frame.
    Rcpp::DataFrame export_augmented_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const Augmented &augm = mcmc_ptr->augm;

        size_t n_cols = augm.scheme().list_latent_export.size();
        size_t n_rows = data_wide.report_indiv.total;

        std::vector<Rcpp::NumericVector> latent_cols(n_cols);
        for (size_t i = 0; i < n_cols; ++i)
            latent_cols[i] = Rcpp::NumericVector(n_rows);

        const std::vector<std::string> data_names =
            vect::difference(data_wide.khold().all(), augm.khold(augm.scheme().list_latent_export));

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &indiv : house)
            {
                size_t c = 0;
                for (size_t a : augm.scheme().list_latent_export)
                    latent_cols[c++][row] = hsp_to_rcpp(mcmc_ptr->augm.at(a, indiv));
                ++row;
            }

        // discrete_export_cols starts with the list_latent_export columns
        const std::vector<bool> discr_cols = augm.discrete_export_cols();

        Rcpp::List output = helper_export_data_wide(mcmc_ptr->data, data_names, true);
        size_t i = 0;
        for (size_t a : augm.scheme().list_latent_export)
        {
            output.push_back(rcpp_typed_column(latent_cols[i], discr_cols[i]), augm.khold(a));
            ++i;
        }
        return rcpp_data_frame(output);
    }

    //////  MCMCBase : data::Long : exports  //////

    /// Export the long data to an R data frame.
    Rcpp::DataFrame export_data_long(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const data::Long &data_long = mcmc_ptr->data.long_format();

        if (data_long.indiv_datalist.size() == 0)
            return Rcpp::DataFrame();

        // retrieve total number of rows
        size_t n_rows = 0;
        for (const auto &tbl : data_long.indiv_datalist)
            if (!tbl.empty())
                n_rows += tbl.begin()->size();

        Rcpp::CharacterVector datalist_index(n_rows);

        // other non-standard data-long columns
        std::vector<Rcpp::NumericVector> datalist_dbl(data_long.khold().size());
        for (size_t d = 0; d < datalist_dbl.size(); ++d)
            datalist_dbl[d] = Rcpp::NumericVector(n_rows);

        size_t row = 0;
        for (const Household &house : data_wide)
            for (const Individual &indiv : house)
            {
                const data::Long::IndivLong *ptr = data_long.get_data_long_ptr(indiv);
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

        return rcpp_data_frame(output);
    }

    //////  MCMCBase : report exports  //////

    /// Summarise household cases to an R data frame.
    Rcpp::DataFrame summarise_house_cases(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const data::ReportIndividual &report_indiv = data_wide.report_indiv;
        const data::ReportHousehold &report_house = data_wide.report_house;

        Rcpp::CharacterVector level = Rcpp::wrap(
            vect::concatenate(std::vector<std::string>(4, "individual"),
                              std::vector<std::string>(5, "household")));
        Rcpp::CharacterVector metric = {"total", "missing", "not infected", "infected",
                                        "total", "no cases", "at least 1 case", "at least 2 cases", "classes"};
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
        // household no cases
        count[5] = report_house.no_cases;
        percent[5] = format_percentage(report_house.no_cases, report_house.total);
        // household at least 1 case
        count[6] = report_house.at_least_1_case;
        percent[6] = format_percentage(report_house.at_least_1_case, report_house.total);
        // household at least 2 cases
        count[7] = report_house.at_least_2_cases;
        percent[7] = format_percentage(report_house.at_least_2_cases, report_house.total);
        // household equivalence classes
        count[8] = report_house.distinct;
        percent[8] = NA_REAL;

        return Rcpp::DataFrame::create(
            Rcpp::_["level"] = std::move(level),
            Rcpp::_["metric"] = std::move(metric),
            Rcpp::_["count"] = std::move(count),
            Rcpp::_["percent"] = std::move(percent));
    }

    //////  Posterior : Rcpp export methods  //////

    /// Move `src` into a fresh NumericVector and release it.
    /// One allocation, one copy, no enclosing matrix.
    Rcpp::NumericVector flush_vector_to_numvec(std::vector<double> &src)
    {
        Rcpp::NumericVector col(src.begin(), src.end());
        std::vector<double>().swap(src);
        return col;
    }

    /// Attach data.frame class and row names to `cols`, a list of equal-length columns.
    void helper_set_dataframe_attrs(Rcpp::List &cols,
                                    R_xlen_t n_rows)
    {
        cols.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(n_rows));
        cols.attr("class") = "data.frame";
    }

    /// Flush each column of `draws` into a data.frame column, releasing C++ storage.
    /// Skipping the NumericMatrix intermediate lets R use the columns in place.
    Rcpp::List flush_2d_draws_to_df(std::vector<std::vector<double>> &draws,
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

    /// Flush parameter value draws as a data.frame.
    Rcpp::List flush_param_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        return flush_2d_draws_to_df(post.param_values, post.param_names,
                                    post.sample_size(), "parameter");
    }

    /// Flush parameter log-prior draws as a data.frame.
    Rcpp::List flush_param_log_prior_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        return flush_2d_draws_to_df(post.param_priors, post.param_prior_names,
                                    post.sample_size(), "parameter log-prior");
    }

    /// Flush latent draws as a data.frame.
    Rcpp::List flush_latent_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;

        if (post.latent_values.empty())
            return Rcpp::List();

        size_t n_cols = post.latent_col_names.size();
        Rcpp::List cols(n_cols);

        size_t col = 0;
        for (auto &indiv_matrix : post.latent_values)
        {
            for (auto &var_vec : indiv_matrix)
                cols[col++] = flush_vector_to_numvec(var_vec);
            std::vector<std::vector<double>>().swap(indiv_matrix);
        }
        std::vector<std::vector<std::vector<double>>>().swap(post.latent_values);

        cols.attr("names") = Rcpp::wrap(post.latent_col_names);
        helper_set_dataframe_attrs(cols, static_cast<R_xlen_t>(post.sample_size()));
        return cols;
    }

    /// Flush pointwise household log-likelihood draws as a data.frame.
    ///
    /// Each column is the complete-data household log-likelihood: transmission, observation, latent.
    /// Observation nodes force unique equivalence classes: every store aligns with transmission.
    /// Households of one equivalence class share a NumericVector SEXP.
    /// R's reference counting keeps storage at n_classes * n_draws until a consumer materialises it.
    Rcpp::List flush_house_log_lik_draws(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        Posterior &post = mcmc_ptr->posterior;
        const data::Wide &data_wide = mcmc_ptr->data.wide();

        size_t n_classes = post.house_log_lik.size();
        size_t n_houses = data_wide.report_house.total;

        // one NumericVector per class, its terms summed at store time
        std::vector<Rcpp::NumericVector> class_cols(n_classes);
        for (size_t R = 0; R < n_classes; ++R)
            class_cols[R] = flush_vector_to_numvec(post.house_log_lik[R]);
        std::vector<std::vector<double>>().swap(post.house_log_lik);

        // identity case (n_classes == n_houses) is a no-op
        Rcpp::List cols(n_houses);
        for (size_t H = 0; H < n_houses; ++H)
            cols[H] = class_cols[data_wide.get_eq_class(H)];

        cols.attr("names") = Rcpp::wrap(post.house_names);
        helper_set_dataframe_attrs(cols, static_cast<R_xlen_t>(post.sample_size()));
        return cols;
    }

    //////  MCMC : house state snapshot  //////

    /// Snapshot the current parameter values, shared by every household.
    Rcpp::NumericVector rcpp_snapshot_param(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const Parameter &param = mcmc_ptr->param;

        if (!param.is_allocated())
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

    /// Snapshot the data of `house`.
    Rcpp::DataFrame rcpp_snapshot_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr,
                                            const Household &house)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const std::vector<std::string> &headers = data_wide.khold().all();

        size_t n_rows = house.size();
        Rcpp::CharacterVector id_col(n_rows);
        Rcpp::CharacterVector hhid_col(n_rows);
        Rcpp::IntegerVector status_col(n_rows);

        std::vector<Rcpp::NumericVector> df_double(headers.size());
        for (size_t d = 0; d < headers.size(); ++d)
            df_double[d] = Rcpp::NumericVector(n_rows);

        for (const Individual &indiv : house)
        {
            size_t i = Index::member(indiv);
            id_col[i] = data_wide.id(indiv);
            hhid_col[i] = data_wide.hhid(house);
            status_col[i] = hsp_to_rcpp(indiv);

            for (size_t d = 0; d < headers.size(); ++d)
                df_double[d][i] = hsp_to_rcpp(mcmc_ptr->data.at(d, indiv).value());
        }

        Rcpp::List output;
        output["id"] = std::move(id_col);
        output["hhid"] = std::move(hhid_col);
        output["status"] = std::move(status_col);
        for (size_t d = 0; d < headers.size(); ++d)
            output.push_back(df_double[d], headers[d]);

        return rcpp_data_frame(output);
    }

    /// Snapshot the long data of `house`.
    Rcpp::DataFrame rcpp_snapshot_data_long(Rcpp::XPtr<MCMC> mcmc_ptr,
                                            const Household &house)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const data::Long &data_long = mcmc_ptr->data.long_format();

        if (data_long.indiv_datalist.empty())
            return Rcpp::DataFrame::create();

        size_t n_rows = 0;
        for (const Individual &indiv : house)
        {
            const data::Long::IndivLong *ptr = data_long.get_data_long_ptr(indiv);
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
            const data::Long::IndivLong *ptr = data_long.get_data_long_ptr(indiv);
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

        return rcpp_data_frame(output);
    }

    /// Snapshot the group risks of the members of `house`.
    Rcpp::DataFrame rcpp_snapshot_group_risks(Rcpp::XPtr<MCMC> mcmc_ptr,
                                              const Household &house)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const risks::Group &group_risk = mcmc_ptr->risks.group();
        const std::vector<std::string> &headers = group_risk.khold().all();

        if (headers.empty() || group_risk.is_empty())
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

        return rcpp_data_frame(output);
    }

    /// Snapshot the contact risks of `house` as a from-by-to matrix, with member ids as dimnames.
    /// - no contact risk defined: empty 0x0 matrix.
    Rcpp::NumericMatrix rcpp_snapshot_contact_risks(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                    const Household &house)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const risks::Contact &contact_risk = mcmc_ptr->risks.contact();

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

    /// Snapshot the infectivity profile of `house`, as resolved for each infected member.
    Rcpp::List rcpp_snapshot_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr,
                                           const Household &house)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const InfecProfile &infec_profile = mcmc_ptr->infec_profile;
        const bool is_latent = infec_profile.is_latent();

        if (!infec_profile.generic && !is_latent)
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

        // a latent profile reads NA until its storage is allocated
        if (!is_latent || mcmc_ptr->augm.is_allocated())
            for (MemberIdx i : house.infected_members())
            {
                const Function *f =
                    is_latent ? mcmc_ptr->augm.at(infec_profile.latent_key, house.at(i)).function()
                              : infec_profile.generic.get();
                if (!f)
                    continue;
                infprof_name[i] = f->name();
                infprof_discr[i] = f->is_discretized();
                infprof_args[i] = rcpp_export_args(f);
                infprof_min[i] = hsp_to_rcpp(f->min());
                infprof_max[i] = hsp_to_rcpp(f->max());
            }

        return Rcpp::List::create(
            Rcpp::_["id"] = id_col,
            Rcpp::_["discr"] = infprof_discr,
            Rcpp::_["func"] = infprof_name,
            Rcpp::_["func_args"] = infprof_args,
            Rcpp::_["min"] = infprof_min,
            Rcpp::_["max"] = infprof_max);
    }

    /// Snapshot the log-density of each key of `list`, per member of `house`.
    /// - key with no distribution for a member: NA.
    Rcpp::DataFrame helper_log_lik_rcpp_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                 const Household &house,
                                                 const std::vector<size_t> &list)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const Augmented &augm = mcmc_ptr->augm;

        size_t n_cols = list.size();

        if (n_cols == 0 || !augm.is_allocated())
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

            for (size_t c = 0; c < n_cols; ++c)
            {
                const Variable &avar = mcmc_ptr->augm.at(list[c], indiv_i);
                if (avar.has_distr())
                    df_double[c][i] = hsp_to_rcpp(avar.guarded_log_pdf());
            }
        }

        Rcpp::List output;
        output["id"] = id_col;
        for (size_t c = 0; c < n_cols; ++c)
            output.push_back(df_double[c], augm.khold(list[c]));

        return rcpp_data_frame(output);
    }

    /// Snapshot the latents of `house`.
    Rcpp::DataFrame helper_latent_rcpp_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                const Household &house)
    {
        const data::Wide &data_wide = mcmc_ptr->data.wide();
        const Augmented &augm = mcmc_ptr->augm;
        size_t n_cols = augm.scheme().list_latent_export.size();

        if (n_cols == 0 || !augm.is_allocated())
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
            for (size_t a : augm.scheme().list_latent_export)
                df_double[c++][i] = hsp_to_rcpp(mcmc_ptr->augm.at(a, indiv_i));
        }

        // discrete_export_cols starts with the list_latent_export columns
        const std::vector<bool> discr_cols = augm.discrete_export_cols();

        Rcpp::List output;
        output.push_back(id_col, "id");
        size_t i = 0;
        for (size_t a : augm.scheme().list_latent_export)
        {
            output.push_back(rcpp_typed_column(df_double[i], discr_cols[i]), augm.khold(a));
            ++i;
        }

        return rcpp_data_frame(output);
    }

    /// Export every snapshot of household `hhid`, for debugging.
    Rcpp::List export_house_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr, const std::string &hhid)
    {
        try
        {
            if (mcmc_ptr->data.wide().is_empty())
                return Rcpp::List::create();

            const Household &house = mcmc_ptr->data.wide().at(hhid);
            const Augmented &augm = mcmc_ptr->augm;

            return Rcpp::List::create(
                Rcpp::_["data_wide"] = rcpp_snapshot_data_wide(mcmc_ptr, house),
                Rcpp::_["data_long"] = rcpp_snapshot_data_long(mcmc_ptr, house),
                Rcpp::_["parameters"] = rcpp_snapshot_param(mcmc_ptr),
                Rcpp::_["group_risks"] = rcpp_snapshot_group_risks(mcmc_ptr, house),
                Rcpp::_["contact_risks"] = rcpp_snapshot_contact_risks(mcmc_ptr, house),
                Rcpp::_["latent"] = helper_latent_rcpp_snapshot(mcmc_ptr, house),
                Rcpp::_["infectivity_profiles"] = rcpp_snapshot_infec_profile(mcmc_ptr, house),
                Rcpp::_["latent_log_lik"] = helper_log_lik_rcpp_snapshot(mcmc_ptr, house, augm.scheme().list_latent_export),
                Rcpp::_["observation_log_lik"] = helper_log_lik_rcpp_snapshot(mcmc_ptr, house, augm.scheme().list_obs_member),
                Rcpp::_["transmission_log_lik"] = mcmc_ptr->data.state.transm_log_lik(house).value());
        }
        catch (const std::exception &e)
        {
            msg::stop(msg::join(e.what(), "for household " + msg::val(hhid)));
        }
    }

    //////  MCMC : simulate  //////

    namespace simulate
    {
        /// Simulate one repetition over every household.
        void single_repetition(MCMC &engine,
                               Rcpp::IntegerVector &status_col,
                               std::vector<Rcpp::NumericVector> &latent_cols,
                               const std::vector<size_t> &row_offsets)
        {
            data::Wide &data_wide = engine.data.wide();
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
                            for (size_t a : augm.scheme().list_latent_export)
                                latent_cols[c++][row] = hsp_to_rcpp(engine.augm.at(a, indiv));
                            // observation nodes: generated for infected members, NA otherwise
                            for (size_t a : augm.scheme().list_obs_member)
                                latent_cols[c++][row] = indiv.is_infected_status()
                                                            ? hsp_to_rcpp(engine.augm.at(a, indiv))
                                                            : NA_REAL;
                            ++row;
                        }

                        // surrogate columns:
                        // - household reduction: one value shared by every member row;
                        // - member column: read on each member, never-infected ones included, as the model reads them.
                        size_t c = augm.scheme().list_latent_export.size() + augm.scheme().list_obs_member.size();
                        for (const augm::Surrogate &surr : augm.scheme().list_surrogate)
                        {
                            if (!surr.is_own_column())
                                continue;
                            if (surr.is_house)
                            {
                                const double v = hsp_to_rcpp(augm.state.house_varia[H][surr.slot].value());
                                for (size_t r = row_offsets[H]; r < row_offsets[H + 1]; ++r)
                                    latent_cols[c][r] = v;
                            }
                            else
                            {
                                const size_t a = augm.khold(augm.def()[surr.row].name);
                                size_t r = row_offsets[H];
                                for (const Individual &indiv : house)
                                    latent_cols[c][r++] = hsp_to_rcpp(engine.augm.at(a, indiv));
                            }
                            ++c;
                        }
                    }
                    catch (const std::exception &e)
                    {
                        throw std::runtime_error(
                            msg::join(e.what(), "for " + data_wide.str(data_wide.at(H))));
                    }
                },
                engine.setup.budget.inner);
        }

        /// Build the simulated data frame: shared columns, status, then exported columns.
        Rcpp::DataFrame build_simulated(const Rcpp::List &shared_df,
                                        const Rcpp::IntegerVector &status_col,
                                        const std::vector<Rcpp::NumericVector> &latent_cols,
                                        const std::vector<bool> &discr_cols,
                                        const Rcpp::CharacterVector &sim_names)
        {
            // build simulated list
            const size_t n_shared(shared_df.size());
            const size_t n_latent_cols(latent_cols.size());
            Rcpp::List sim_list(n_shared + 1 + n_latent_cols);
            for (size_t i = 0; i < n_shared; ++i)
                sim_list[i] = shared_df[i];
            sim_list[n_shared] = status_col;
            for (size_t c = 0; c < n_latent_cols; ++c)
                sim_list[n_shared + 1 + static_cast<int>(c)] =
                    rcpp_typed_column(latent_cols[c], discr_cols[c]);
            sim_list.attr("names") = sim_names;
            return rcpp_data_frame(sim_list);
        }
    }

    /// Run the simulation, returning wide and long results after post-processing.
    Rcpp::List run_simulation(Rcpp::XPtr<MCMC> mcmc_ptr,
                              Rcpp::Nullable<Rcpp::Function> obs_process,
                              Rcpp::Nullable<Rcpp::Function> apply_wide,
                              Rcpp::Nullable<Rcpp::Function> obs_process_and_apply_long)
    {
        MCMC &engine = *mcmc_ptr;
        Setup &setup = engine.setup;
        const Augmented &augm = engine.augm;
        const data::Wide &data_wide = engine.data.wide();

        try
        {
            engine.model->init_model_for_simulation(); // model-specific initialization
            engine.augm.init_house_for_simulation();   // build household-level observations
            data_wide.assert_required_everyone(engine.data.state, augm.get_live_cols());
            engine.risks.assert_identifiability(engine.param);

            // @todo Assert the infectivity profile here.
            // latents are not built yet, so a latent profile cannot be checked
            // engine.assert_valid_infectivity_profile(); // assert profile is valid

            //////  Observation Process & Measures  //////

            // @todo Hold nullable R functions without a pointer work-around.
            std::unique_ptr<Rcpp::Function> f_obs_process = nullptr, f_wide = nullptr, f_long = nullptr;

            if (!obs_process.isNull())
                f_obs_process = std::make_unique<Rcpp::Function>(obs_process);
            if (!apply_wide.isNull())
                f_wide = std::make_unique<Rcpp::Function>(apply_wide);
            if (!obs_process_and_apply_long.isNull())
                f_long = std::make_unique<Rcpp::Function>(obs_process_and_apply_long);

            //////  Initialization  //////

            Rcpp::List out_wide(setup.n_steps); // allocate memory for wide results
            Rcpp::List out_long;                // long results are optional
            if (f_long)
                out_long = Rcpp::List(setup.n_steps);

            const std::vector<std::string> static_data_names = engine.get_static_data_names();
            const Rcpp::List shared_df = helper_export_data_wide(engine.data, static_data_names, false);

            // simulated data column names
            Rcpp::CharacterVector sim_names = Rcpp::wrap(engine.get_col_names());
            // row offsets for house parallelization
            const std::vector<size_t> row_offsets = engine.get_row_offsets();

            // pre-allocate output columns once, reused across repetitions
            // they hold exported latent, observation and surrogate columns
            Rcpp::IntegerVector status_col(data_wide.report_indiv.total);

            std::vector<Rcpp::NumericVector> latent_cols(augm.n_export_cols());
            for (size_t c = 0; c < latent_cols.size(); ++c)
                latent_cols[c] = Rcpp::NumericVector(data_wide.report_indiv.total);

            // flags for discrete columns, set once variables are bound after the first repetition
            std::vector<bool> discr_cols(latent_cols.size(), false);

            //////  Simulations  //////

            setup.prog_open_for_simulation();
            for (size_t r = 0; r < setup.n_steps; ++r, setup.prog_update())
            {
                engine.augm.clear_all();
                engine.data.wide().recover_all_status();
                engine.param.load_sample_values(r);
                engine.param.draw_sample(setup.prng()); // draw parameters with no pinned value

                // run one repetition
                simulate::single_repetition(
                    engine, status_col, latent_cols, row_offsets);

                if (r == 0)
                    discr_cols = augm.discrete_export_cols();

                // build simulated list
                Rcpp::DataFrame simulated(simulate::build_simulated(
                    shared_df, status_col, latent_cols, discr_cols, sim_names));

                // observation process
                if (f_obs_process)
                    simulated = (*f_obs_process)(simulated);

                // apply the measure to data-wide, or store it
                if (f_wide)
                    out_wide[r] = (*f_wide)(simulated);
                else
                    out_wide[r] = Rcpp::clone(simulated); // clone when storing
                // long observation process and apply measure
                if (f_long)
                    out_long[r] = (*f_long)(simulated);
            }
            setup.prog_close();

            // restore observed status and point every read-cell back at its imported value
            // post-simulation reads of data_wide then return the input
            // generated values leave through the simulated frames
            // the wired cells refer to Variables this run is about to clear
            engine.augm.clear_all();
            engine.data.wide().recover_all_status();
            engine.data.reset_readcells();

            return Rcpp::List::create(Rcpp::_["wide"] = out_wide,
                                      Rcpp::_["long"] = out_long);
        }
        catch (std::exception &e)
        {
            setup.crash(e);
        }
    }

    //////////////////////////////  Rcpp exports  //////////////////////////////

/// Keys the model declares through HOUSEPI_DECLARE_REGISTER_KEYS.
#define MAKE_EXPORT_CHARACTER_VECTOR(method)                \
    Rcpp::CharacterVector method(Rcpp::XPtr<MCMC> mcmc_ptr) \
    {                                                       \
        return Rcpp::wrap(mcmc_ptr->model->method());       \
    }

    MAKE_EXPORT_CHARACTER_VECTOR(key_Data)
    MAKE_EXPORT_CHARACTER_VECTOR(key_Long)
    MAKE_EXPORT_CHARACTER_VECTOR(key_Param)
    MAKE_EXPORT_CHARACTER_VECTOR(key_GroupRisk)
    MAKE_EXPORT_CHARACTER_VECTOR(key_Latent)

/// Keys the model requests as discrete.
#define MAKE_EXPORT_DISCRETE_KEYS(method, module)            \
    Rcpp::CharacterVector method(Rcpp::XPtr<MCMC> mcmc_ptr)  \
    {                                                        \
        return Rcpp::wrap(mcmc_ptr->module.discrete_keys()); \
    }

    MAKE_EXPORT_DISCRETE_KEYS(key_Param_discrete, param)
    MAKE_EXPORT_DISCRETE_KEYS(key_Latent_discrete, augm)

    /// Data columns the model requests as discrete.
    Rcpp::CharacterVector key_Data_discrete(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return Rcpp::wrap(mcmc_ptr->data.wide().discrete_keys());
    }

    /// Baselines the model anchors to risks, keyed by parameter name.
    Rcpp::List risk_anchors(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const std::vector<risks::Anchor> &anchors = mcmc_ptr->risks.anchors();

        Rcpp::List out(anchors.size());
        Rcpp::CharacterVector names(anchors.size());
        for (size_t a = 0; a < anchors.size(); ++a)
        {
            names[a] = mcmc_ptr->param.khold(anchors[a].baseline);
            out[a] = Rcpp::List::create(
                Rcpp::_["risks"] = Rcpp::wrap(mcmc_ptr->risks.group().khold(anchors[a].risks)),
                Rcpp::_["contacts"] = anchors[a].contacts);
        }
        out.names() = names;
        return out;
    }

    /// Group risks the model reads on members with infected status only.
    Rcpp::CharacterVector risk_status_stratified(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        const risks::Group &group = mcmc_ptr->risks.group();
        return Rcpp::wrap(group.khold(group.status_stratified()));
    }

    Rcpp::XPtr<MCMC> create_mcmc_runner(std::string module_name)
    {
        Rcpp::XPtr<MCMC> mcmc(new MCMC(module_name), true);

        // options passed to the progress setup
        const unsigned width = Rcpp::as<unsigned>(Rcpp::Function("getOption")("width"));
        const bool is_interactive = Rcpp::as<bool>(Rcpp::Function("interactive")());
        mcmc->setup.init_R_context(
            std::min(width, 78u) / 2,
            is_interactive,
            []
            { Rcpp::checkUserInterrupt(); });

        return mcmc;
    }

    /// Share the imported tables of `src` with `dest`, copy-on-write.
    void share_data_store(Rcpp::XPtr<MCMC> dest, Rcpp::XPtr<MCMC> src)
    {
        dest->data.share(src->data);
    }

    bool has_contact_risk(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return mcmc_ptr->model->has_contact_risk();
    }

    bool has_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        return mcmc_ptr->model->has_infec_profile();
    }

    double evaluate_log_lik(Rcpp::XPtr<MCMC> mcmc_ptr, std::string hhid)
    {
        try
        {
            return mcmc_ptr->evaluate_log_lik(hhid);
        }
        catch (const std::exception &e)
        {
            mcmc_ptr->setup.crash(e);
        }
    }
    double run_inference(Rcpp::XPtr<MCMC> mcmc_ptr)
    {
        try
        {
            return mcmc_ptr->run_inference();
        }
        catch (std::exception &e)
        {
            mcmc_ptr->setup.crash(e);
        }
    }

    /// Aggregate individuals by hhid, given as 1-based codes of the distinct ids.
    /// - missing status: NA (NA_INTEGER is INT_MIN).
    measure::cases::Table summarise_houses_from_df(Rcpp::DataFrame data_wide)
    {
        Rcpp::CharacterVector col_names = data_wide.names();

        std::vector<int> hhid = {};
        if (vect::contains(col_names, "hhid"))
            hhid = Rcpp::as<std::vector<int>>(data_wide["hhid"]);

        std::vector<int> status = {};
        if (vect::contains(col_names, "status"))
            status = Rcpp::as<std::vector<int>>(data_wide["status"]);

        return measure::summarise_houses(hhid, status);
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
        measure::cases::AttackRates res = measure::attack_rate_by_hhsize(summarise_houses_from_df(data_wide));
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
        measure::cases::SecAttackRates res = measure::sar_by_hhsize(summarise_houses_from_df(data_wide));
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
        measure::cases::Table res = summarise_houses_from_df(data_wide);
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

    double measure_geom_mean_hhsize(Rcpp::DataFrame data_wide)
    {
        measure::cases::Table table = summarise_houses_from_df(data_wide);
        return table.max_cases > 1 ? measure::geom_mean_hhsize(table) : NA_REAL;
    }
}

RCPP_MODULE(housepi)
{
    Rcpp::class_<ExposeBase>("ExposeBase")
        .method("args", &ExposeBase::args)
        .method("domain", &ExposeBase::domain)
        .method("bound_aliases", &ExposeBase::bound_aliases)
        .method("is_discrete", &ExposeBase::is_discrete)
        .method("is_positive", &ExposeBase::is_positive)
        .method("codomain", &ExposeBase::codomain)
        .method("loc", &ExposeBase::loc);

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
    Rcpp::function("key_Latent", &key_Latent);
    Rcpp::function("key_Param_discrete", &key_Param_discrete);
    Rcpp::function("key_Latent_discrete", &key_Latent_discrete);
    Rcpp::function("key_Data_discrete", &key_Data_discrete);
    Rcpp::function("risk_anchors", &risk_anchors);
    Rcpp::function("risk_status_stratified", &risk_status_stratified);

    Rcpp::function("measure_sar_by_hhsize", &measure_sar_by_hhsize);
    Rcpp::function("measure_ar_by_hhsize", &measure_ar_by_hhsize);
    Rcpp::function("measure_cases_by_hhsize", &measure_cases_by_hhsize);
    Rcpp::function("measure_geom_mean_hhsize", &measure_geom_mean_hhsize);
    Rcpp::function("count_houses_by_hhsize", &count_houses_by_hhsize);

    Rcpp::function("aux_chain_binomial_debugger", &housepi::chain_binom::aux_debugger);

    Rcpp::function("create_mcmc_runner", &create_mcmc_runner);
    Rcpp::function("run_inference", &run_inference);
    Rcpp::function("run_simulation", &run_simulation);
    Rcpp::function("evaluate_log_lik", &evaluate_log_lik);
    Rcpp::function("import_data_wide", &import_data_wide);
    Rcpp::function("import_data_long", &import_data_long);
    Rcpp::function("share_data_store", &share_data_store);
    Rcpp::function("import_group_risk", &import_group_risk);
    Rcpp::function("import_contact_risk", &import_contact_risk);
    Rcpp::function("import_setup", &import_setup);
    Rcpp::function("import_param", &import_param);
    Rcpp::function("import_infec_profile", &import_infec_profile);
    Rcpp::function("import_latent", &import_latent);
    Rcpp::function("import_param_table", &import_param_table);
    Rcpp::function("has_contact_risk", &has_contact_risk);
    Rcpp::function("has_infec_profile", &has_infec_profile);
    Rcpp::function("export_data_wide", &export_data_wide);
    Rcpp::function("export_data_long", &export_data_long);
    Rcpp::function("summarise_house_cases", &summarise_house_cases);
    Rcpp::function("export_setup", &export_setup);
    Rcpp::function("export_augmented_data_wide", &export_augmented_data_wide);
    Rcpp::function("flush_param_draws", &flush_param_draws);
    Rcpp::function("flush_param_log_prior_draws", &flush_param_log_prior_draws);
    Rcpp::function("flush_latent_draws", &flush_latent_draws);
    Rcpp::function("flush_house_log_lik_draws", &flush_house_log_lik_draws);
    Rcpp::function("export_house_snapshot", &export_house_snapshot);
}
#endif
