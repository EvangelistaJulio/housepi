/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef _STANDALONE
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include <Rcpp.h>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/progressbar.h"
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
    /// Exception raised by housepi code, its message written as cli markup.
    ///
    /// - any other exception: foreign, its message rendered literally in R.
    /// @note Keep its qualified name in sync with the R side.
    struct TrustedException : std::exception
    {
        std::string m;
        explicit TrustedException(std::string s) : m(std::move(s)) {}
        const char *what() const noexcept override { return m.c_str(); }
    };
}

namespace housepi::msg
{
    /// Throw `txt` as a TrustedException.
    [[noreturn]] void stop(const std::string &txt)
    {
        throw housepi::TrustedException(txt);
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
        if (WarningBuffer::hold(txt))
            return;
        Rcpp::List cond = Rcpp::List::create(Rcpp::_["message"] = txt,
                                             Rcpp::_["call"] = R_NilValue);
        cond.attr("class") = Rcpp::CharacterVector::create(
            "housepi_cpp_warning", "warning", "condition");
        Rcpp::Function("warning")(cond);
    }
}

//////////////////////////////  Rcpp helpers  //////////////////////////////

namespace housepi
{
    namespace
    {
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

        /// Export the named arguments of `func`.
        Rcpp::NumericVector rcpp_export_args(const Function *func)
        {
            Rcpp::NumericVector vals(func->n_args());
            Rcpp::CharacterVector keys(func->n_args());
            for (size_t k = 0; k < func->n_args(); ++k)
            {
                vals[k] = hsp_to_rcpp(func->arg(k));
                keys[k] = func->arg_name(k);
            }
            vals.names() = std::move(keys);
            return vals;
        }

        /// Convert `input` to an IntegerVector when `discr` is requested.
        /// Values are assumed whole.
        /// - value beyond the int range: `input` returned unchanged.
        SEXP rcpp_typed_column(const Rcpp::NumericVector &input, bool discr)
        {
            constexpr int machine_max_int = std::numeric_limits<int>::max();

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
        /// Hand-built so that column names pass verbatim.
        /// @note Rcpp's DataFrame would sanitise them.
        Rcpp::DataFrame rcpp_data_frame(Rcpp::List output,
                                        Rcpp::CharacterVector names)
        {
            output.names() = names;
            R_xlen_t n_rows = output.size() ? Rf_xlength(output[0]) : 0;
            // compact row names, the form R stores for the 1..n sequence
            output.attr("row.names") = Rcpp::IntegerVector::create(
                NA_INTEGER, -static_cast<int>(n_rows));
            output.attr("class") = "data.frame";
            return output;
        }

        /// Allocate `n_cols` numeric columns of `n_rows` NA values.
        std::vector<Rcpp::NumericVector> rcpp_columns(size_t n_cols, size_t n_rows)
        {
            std::vector<Rcpp::NumericVector> out;
            out.reserve(n_cols);
            for (size_t c = 0; c < n_cols; ++c)
                out.emplace_back(n_rows, NA_REAL);
            return out;
        }

        /// Fill a fixed number of named column slots, then stamp them as an R data.frame.
        /// - column: protected by its slot from the moment it is added;
        /// - repeated name: a second column, never a replacement.
        struct RcppFrame
        {
            explicit RcppFrame(size_t n_cols)
                : cols_(n_cols), names_(n_cols) {}

            /// Place column `col` named `name` in the next slot.
            /// - slot past the declared width: Rcpp index error.
            void add(const std::string &name, SEXP col)
            {
                names_(next_) = name;
                cols_(next_++) = col;
            }
            /// Place a numeric column of `n_rows` NA values in the next slot, returned for filling.
            Rcpp::NumericVector add_na(const std::string &name, size_t n_rows)
            {
                Rcpp::NumericVector col(n_rows, NA_REAL);
                add(name, col);
                return col;
            }
            /// Build the data.frame from the filled slots.
            Rcpp::DataFrame frame() { return rcpp_data_frame(cols_, names_); }

        private:
            Rcpp::List cols_;             /// Column slots.
            Rcpp::CharacterVector names_; /// Column names, one per slot.
            size_t next_ = 0;             /// Next slot to fill.
        };

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

            /// Arguments, named.
            Rcpp::NumericVector args() const { return rcpp_export_args(ptr_.get()); }

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
            double loc() const { return hsp_to_rcpp(ptr_->loc().value()); }

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
            static std::vector<registry::Kind> reg_kinds() { return registry::func_kinds(); }

            ExposeFunction(std::string name, bool /*discr*/)
                : ExposeBase(std::move(name), family(), reg_kinds()) {}

            MAKE_EXPOSE_METHOD_CALL(func, ptr_->func(v))
            MAKE_EXPOSE_METHOD_CALL(log_func, ptr_->log_func(v))
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

            ExposeDistribution(std::string name, bool discr)
                : ExposeBase(std::move(name), family(), reg_kinds(discr)) {}

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

            ExposeInfecProfile(std::string name, bool discr)
                : ExposeBase(std::move(name), family(), reg_kinds(discr))
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

        /// Convert an R data frame to a table of columns, by name.
        /// - NULL: empty table.
        template <class Table>
        Table rdataframe_to(SEXP input)
        {
            using Column = typename Table::mapped_type;
            if (Rf_isNull(input))
                return {};

            Rcpp::DataFrame df(input);
            const Rcpp::CharacterVector names = df.names();
            Table output;
            for (R_xlen_t i = 0; i < df.size(); ++i)
            {
                Rcpp::Vector<Rcpp::traits::r_sexptype_traits<typename Column::value_type>::rtype> col = df[i];
                output[Rcpp::as<std::string>(names[i])] = Column(col.begin(), col.end());
            }
            return output;
        }

        // import helpers run serially at the Rcpp boundary
        // their build() paths raise housepi::TrustedException through msg::stop
        // Rcpp forwards it to R with its class
        // no import needs its own try/catch conversion
        void import_setup(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->setup.build(rdataframe_to<str_table>(input));
        }

        void import_param(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->param.build(rdataframe_to<str_table>(input), mcmc_ptr->translator);
        }

        void import_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->detach_data_store();
            mcmc_ptr->data.build_wide(rdataframe_to<str_table>(input));
            mcmc_ptr->init_equivalence_class_compression();
        }

        void import_data_long(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            // @todo Decide whether long data should be copy-on-write too.
            mcmc_ptr->data.build_long(rdataframe_to<str_table>(input));
        }

        void import_group_risk(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->risks.build_group(mcmc_ptr->data.wide(), mcmc_ptr->param, rdataframe_to<str_table>(input));
            mcmc_ptr->init_equivalence_class_compression();
        }

        void import_contact_risk(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->risks.build_contact(mcmc_ptr->data.wide(), mcmc_ptr->param, rdataframe_to<str_table>(input));
            mcmc_ptr->init_equivalence_class_compression();
        }

        void import_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->infec_profile.build(mcmc_ptr->translator, rdataframe_to<str_table>(input));
        }

        void import_latent(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->augm.build(mcmc_ptr->infec_profile, mcmc_ptr->status_process, rdataframe_to<str_table>(input));
            mcmc_ptr->init_equivalence_class_compression();
        }

        void import_status_process(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input, size_t n_levels)
        {
            mcmc_ptr->status_process.build(mcmc_ptr->data.wide(), rdataframe_to<str_table>(input), n_levels);
        }

        void import_param_table(Rcpp::XPtr<MCMC> mcmc_ptr, SEXP input)
        {
            mcmc_ptr->param.build_for_predic(rdataframe_to<dbl_table>(input));
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
                Rcpp::_["id_draws"] = setup.id_draws,
                Rcpp::_["seed"] = setup.seed,
                Rcpp::_["time_step"] = setup.time_step);
        }

        //////  MCMC : member ranges  //////

        /// Count number of individuals in the data.
        size_t n_members(const data::Wide &data_wide) { return data_wide.report_indiv.total; }
        /// Count number of members of `house`.
        size_t n_members(const Household &house) { return house.size(); }

        /// Visit every member of `house`, in row order.
        template <class Fn>
        void for_each_member(const Household &house, Fn &&fn)
        {
            for (const Individual &indiv : house)
                fn(house, indiv);
        }
        /// Visit every member of every household, in row order.
        template <class Fn>
        void for_each_member(const data::Wide &data_wide, Fn &&fn)
        {
            for (const Household &house : data_wide)
                for_each_member(house, fn);
        }

        /// Export the member IDs of `scope`.
        /// @tparam Scope A Household or data::Wide.
        template <class Scope>
        Rcpp::CharacterVector rcpp_ids(const data::Wide &data_wide, const Scope &scope)
        {
            Rcpp::CharacterVector out(n_members(scope));
            size_t row = 0;
            for_each_member(scope, [&](const Household &, const Individual &indiv)
                            { out[row++] = data_wide.id(indiv); });
            return out;
        }

        //////  MCMC : data::Wide : exports  //////

        /// Number of columns `add_wide_columns` adds for `n_data` data columns.
        size_t n_wide_cols(size_t n_data, bool is_status) { return 2 + is_status + n_data; }

        /// Add id, hh_id, status when `is_status`, then data `columns`, on the members of `scope`.
        template <class Scope>
        void add_wide_columns(RcppFrame &frame,
                              const Data &data,
                              const Scope &scope,
                              const std::vector<std::string> &columns,
                              bool is_status)
        {
            const data::Wide &data_wide = data.wide();
            const size_t n_rows = n_members(scope);
            const std::vector<size_t> keys = data_wide.consult(columns);
            Rcpp::CharacterVector hh_id_col(n_rows);
            Rcpp::IntegerVector status_col(n_rows);

            frame.add("id", rcpp_ids(data_wide, scope));
            frame.add("hh_id", hh_id_col);
            if (is_status)
                frame.add("status", status_col);
            std::vector<Rcpp::NumericVector> cols;
            cols.reserve(columns.size());
            for (const std::string &column : columns)
                cols.push_back(frame.add_na(column, n_rows));

            size_t row = 0;
            for_each_member(scope, [&](const Household &house, const Individual &indiv)
                            {
                            hh_id_col[row] = data_wide.hh_id(house);
                            status_col[row] = hsp_to_rcpp(indiv);
                            for (size_t d = 0; d < cols.size(); ++d)
                                cols[d][row] = hsp_to_rcpp(data.at(keys[d], indiv).value());
                            ++row; });
        }

        /// Add exported latent columns on the members of `scope`, integer-typed when discrete.
        template <class Scope>
        void add_latent_columns(RcppFrame &frame, const Augmented &augm, const Scope &scope)
        {
            const std::vector<size_t> &list = augm.scheme().list_latent_export;
            std::vector<Rcpp::NumericVector> cols = rcpp_columns(list.size(), n_members(scope));

            size_t row = 0;
            for_each_member(scope, [&](const Household &, const Individual &indiv)
                            {
                            for (size_t c = 0; c < list.size(); ++c)
                                cols[c][row] = hsp_to_rcpp(augm.at(list[c], indiv));
                            ++row; });

            // discrete_export_cols starts with the list_latent_export columns
            const std::vector<bool> discr_cols = augm.discrete_export_cols();
            for (size_t c = 0; c < list.size(); ++c)
                frame.add(augm.consult(list[c]), rcpp_typed_column(cols[c], discr_cols[c]));
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

            const std::vector<std::string> columns = engine.data.wide().imported_cols();
            RcppFrame frame(n_wide_cols(columns.size(), true));
            add_wide_columns(frame, engine.data, engine.data.wide(), columns, true);
            return frame.frame();
        }

        /// Export data and latents at the current state to an R data frame.
        Rcpp::DataFrame export_augmented_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            const data::Wide &data_wide = mcmc_ptr->data.wide();
            const Augmented &augm = mcmc_ptr->augm;

            const std::vector<size_t> &list = augm.scheme().list_latent_export;
            const std::vector<std::string> columns = vect::difference(data_wide.all(), augm.consult(list));
            RcppFrame frame(n_wide_cols(columns.size(), true) + list.size());
            add_wide_columns(frame, mcmc_ptr->data, data_wide, columns, true);
            add_latent_columns(frame, augm, data_wide);
            return frame.frame();
        }

        //////  MCMC : data::Long : exports  //////

        /// Export the long data of the members of `scope`, one row per entry.
        /// - no long data: empty data frame.
        template <class Scope>
        Rcpp::DataFrame rcpp_long_frame(const Data &data, const Scope &scope)
        {
            const data::Wide &data_wide = data.wide();
            const data::Long &data_long = data.long_format();

            if (data_long.indiv_datalist.empty())
                return Rcpp::DataFrame::create();

            // entries of a member: length of its first column
            auto n_entries = [&](const Individual &indiv)
            {
                const data::Long::IndivLong *ptr = data_long.get_data_long_ptr(indiv);
                return ptr->empty() ? size_t(0) : ptr->begin()->size();
            };
            size_t n_rows = 0;
            for_each_member(scope, [&](const Household &, const Individual &indiv)
                            { n_rows += n_entries(indiv); });

            RcppFrame frame(1 + data_long.size());
            Rcpp::CharacterVector id_col(n_rows);
            frame.add("id", id_col);
            std::vector<Rcpp::NumericVector> cols;
            cols.reserve(data_long.size());
            for (size_t d = 0; d < data_long.size(); ++d)
                cols.push_back(frame.add_na(data_long.consult(d), n_rows));

            size_t row = 0;
            for_each_member(scope, [&](const Household &, const Individual &indiv)
                            {
                            const data::Long::IndivLong *ptr = data_long.get_data_long_ptr(indiv);
                            for (size_t k = 0; k < n_entries(indiv); ++k, ++row)
                            {
                                id_col[row] = data_wide.id(indiv);
                                for (size_t d = 0; d < cols.size(); ++d)
                                    cols[d][row] = hsp_to_rcpp(ptr->at(d)[k]);
                            } });
            return frame.frame();
        }

        /// Export the long data to an R data frame.
        Rcpp::DataFrame export_data_long(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            return rcpp_long_frame(mcmc_ptr->data, mcmc_ptr->data.wide());
        }

        //////  MCMCBase : report exports  //////

        /// Format `value` as a percentage of `total`.
        double format_percentage(size_t value, size_t total)
        {
            return std::round(10000. * static_cast<double>(value) / static_cast<double>(total)) / 100.;
        }

        /// Summarise household cases to an R data frame.
        Rcpp::DataFrame summarise_house_cases(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            const data::ReportIndividual &indiv = mcmc_ptr->data.wide().report_indiv;
            const data::ReportHousehold &house = mcmc_ptr->data.wide().report_house;

            auto pct_indiv = [&](size_t n)
            { return format_percentage(n, indiv.total); };
            auto pct_house = [&](size_t n)
            { return format_percentage(n, house.total); };

            struct Row
            {
                const char *level, *metric; /// Row labels.
                size_t count;               /// Individuals or households counted.
                double percent;             /// Share of the level total.
            };
            const std::vector<Row> rows = {
                {"individual", "total", indiv.total, 100.},
                {"individual", "missing", indiv.missing, pct_indiv(indiv.missing)},
                {"individual", "not infected", indiv.escaped, pct_indiv(indiv.escaped)},
                {"individual", "infected", indiv.infectious, pct_indiv(indiv.infectious)},
                {"household", "total", house.total, 100.},
                {"household", "no cases", house.no_cases, pct_house(house.no_cases)},
                {"household", "at least 1 case", house.at_least_1_case, pct_house(house.at_least_1_case)},
                {"household", "at least 2 cases", house.at_least_2_cases, pct_house(house.at_least_2_cases)},
                {"household", "classes", house.distinct, NA_REAL}};

            Rcpp::CharacterVector level(rows.size()), metric(rows.size());
            Rcpp::IntegerVector count(rows.size());
            Rcpp::NumericVector percent(rows.size());
            for (size_t r = 0; r < rows.size(); ++r)
            {
                level[r] = rows[r].level;
                metric[r] = rows[r].metric;
                count[r] = static_cast<int>(rows[r].count);
                percent[r] = rows[r].percent;
            }

            RcppFrame frame(4);
            frame.add("level", level);
            frame.add("metric", metric);
            frame.add("count", count);
            frame.add("percent", percent);
            return frame.frame();
        }

        //////  Posterior : Rcpp export methods  //////

        /// Flush `draws` as a data.frame, releasing its C++ storage column by column.
        /// - missing draw: NA;
        /// - table columns reading one stored column share its SEXP, until R modifies one;
        /// - flushed again: empty data.frame.
        Rcpp::DataFrame flush_draws(Posterior::Draws &draws)
        {
            Rcpp::List stored(draws.cols.size());
            for (size_t c = 0; c < draws.cols.size(); ++c)
            {
                Rcpp::NumericVector col(draws.cols[c].size());
                std::transform(draws.cols[c].begin(), draws.cols[c].end(), col.begin(),
                               [](double x)
                               { return hsp_to_rcpp(x); });
                stored[c] = col;
                std::vector<double>().swap(draws.cols[c]);
            }

            RcppFrame frame(draws.names.size());
            for (size_t k = 0; k < draws.names.size(); ++k)
                frame.add(draws.names[k], VECTOR_ELT(stored, draws.source.empty() ? k : draws.source[k]));
            draws = {};
            return frame.frame();
        }

        /// Flush posterior draws `Member` as a data.frame.
        template <Posterior::Draws Posterior::*Member>
        Rcpp::DataFrame flush_posterior_(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            return flush_draws(mcmc_ptr->posterior.*Member);
        }

        //////  MCMC : house state snapshot  //////

        /// Snapshot the current parameter values, shared by every household.
        Rcpp::NumericVector rcpp_snapshot_param(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            const Parameter &param = mcmc_ptr->param;

            if (!param.is_allocated())
                return Rcpp::NumericVector::create();

            Rcpp::NumericVector output(param.size());
            for (R_xlen_t p = 0; p < output.size(); ++p)
                output[p] = hsp_to_rcpp(param.at(p).value());
            output.names() = Rcpp::wrap(param.all());

            return output;
        }

        /// Snapshot the data of `house`.
        Rcpp::DataFrame rcpp_snapshot_data_wide(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                const Household &house)
        {
            const std::vector<std::string> &columns = mcmc_ptr->data.wide().all();
            RcppFrame frame(n_wide_cols(columns.size(), true));
            add_wide_columns(frame, mcmc_ptr->data, house, columns, true);
            return frame.frame();
        }

        /// Snapshot the group risks of the members of `house`.
        Rcpp::DataFrame rcpp_snapshot_group_risks(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                  const Household &house)
        {
            const risks::Group &group_risk = mcmc_ptr->risks.group();
            const std::vector<std::string> &headers = group_risk.all();

            if (headers.empty() || group_risk.is_empty())
                return Rcpp::DataFrame::create();

            RcppFrame frame(1 + headers.size());
            frame.add("id", rcpp_ids(mcmc_ptr->data.wide(), house));
            for (size_t c = 0; c < headers.size(); ++c)
            {
                Rcpp::NumericVector col = frame.add_na(headers[c], house.size());
                for (size_t i = 0; i < house.size(); ++i)
                    col[i] = hsp_to_rcpp(group_risk.at(c, house.at(i)));
            }
            return frame.frame();
        }

        /// Snapshot the contact risks of `house` as a from-by-to matrix, with member ids as dimnames.
        /// - no contact risk defined: empty 0x0 matrix.
        Rcpp::NumericMatrix rcpp_snapshot_contact_risks(Rcpp::XPtr<MCMC> mcmc_ptr,
                                                        const Household &house)
        {
            const risks::Contact &contact_risk = mcmc_ptr->risks.contact();

            if (contact_risk.is_empty())
                return Rcpp::NumericMatrix(0, 0);

            const size_t n = house.size();
            Rcpp::NumericMatrix mtx(n, n);
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j < n; ++j)
                    mtx(i, j) = hsp_to_rcpp(contact_risk.value(house, i, j));

            const Rcpp::CharacterVector ids = rcpp_ids(mcmc_ptr->data.wide(), house);
            mtx.attr("dimnames") = Rcpp::List::create(ids, ids);
            return mtx;
        }

        /// Snapshot the infectivity profile of `house`, as resolved for each infected member.
        Rcpp::List rcpp_snapshot_infec_profile(Rcpp::XPtr<MCMC> mcmc_ptr,
                                               const Household &house)
        {
            const InfecProfile &infec_profile = mcmc_ptr->infec_profile;
            const bool is_latent = infec_profile.is_latent();

            if (!infec_profile.generic && !is_latent)
                return Rcpp::List::create();

            size_t n_rows = house.size();
            Rcpp::CharacterVector infprof_name(n_rows, NA_STRING);
            Rcpp::LogicalVector infprof_discr(n_rows, NA_LOGICAL);
            Rcpp::List infprof_args(n_rows, Rcpp::NumericVector(0));
            Rcpp::NumericVector infprof_min(n_rows, NA_REAL);
            Rcpp::NumericVector infprof_max(n_rows, NA_REAL);

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
                Rcpp::_["id"] = rcpp_ids(mcmc_ptr->data.wide(), house),
                Rcpp::_["discr"] = infprof_discr,
                Rcpp::_["func"] = infprof_name,
                Rcpp::_["func_args"] = infprof_args,
                Rcpp::_["min"] = infprof_min,
                Rcpp::_["max"] = infprof_max);
        }

        /// Snapshot the log-density of each key of `list`, per member of `house`.
        /// - key with no distribution for a member: NA.
        Rcpp::DataFrame rcpp_snapshot_log_lik(Rcpp::XPtr<MCMC> mcmc_ptr,
                                              const Household &house,
                                              const std::vector<size_t> &list)
        {
            const Augmented &augm = mcmc_ptr->augm;

            if (list.empty() || !augm.is_allocated())
                return Rcpp::DataFrame::create();

            RcppFrame frame(1 + list.size());
            frame.add("id", rcpp_ids(mcmc_ptr->data.wide(), house));
            for (size_t c = 0; c < list.size(); ++c)
            {
                Rcpp::NumericVector col = frame.add_na(augm.consult(list[c]), house.size());
                for (size_t i = 0; i < house.size(); ++i)
                {
                    const Variable &avar = augm.at(list[c], house.at(i));
                    if (avar.has_distr())
                        col[i] = hsp_to_rcpp(avar.guarded_log_pdf());
                }
            }
            return frame.frame();
        }

        /// Snapshot the latents of `house`.
        Rcpp::DataFrame rcpp_snapshot_latent(Rcpp::XPtr<MCMC> mcmc_ptr,
                                             const Household &house)
        {
            const Augmented &augm = mcmc_ptr->augm;

            if (augm.scheme().list_latent_export.empty() || !augm.is_allocated())
                return Rcpp::DataFrame::create();

            RcppFrame frame(1 + augm.scheme().list_latent_export.size());
            frame.add("id", rcpp_ids(mcmc_ptr->data.wide(), house));
            add_latent_columns(frame, augm, house);
            return frame.frame();
        }

        /// Export every snapshot of household `hh_id`, for debugging.
        Rcpp::List export_house_snapshot(Rcpp::XPtr<MCMC> mcmc_ptr, const std::string &hh_id)
        {
            try
            {
                if (mcmc_ptr->data.wide().is_empty())
                    return Rcpp::List::create();

                const Household &house = mcmc_ptr->data.wide().at(hh_id);
                const Augmented &augm = mcmc_ptr->augm;

                return Rcpp::List::create(
                    Rcpp::_["data_wide"] = rcpp_snapshot_data_wide(mcmc_ptr, house),
                    Rcpp::_["data_long"] = rcpp_long_frame(mcmc_ptr->data, house),
                    Rcpp::_["parameters"] = rcpp_snapshot_param(mcmc_ptr),
                    Rcpp::_["group_risks"] = rcpp_snapshot_group_risks(mcmc_ptr, house),
                    Rcpp::_["contact_risks"] = rcpp_snapshot_contact_risks(mcmc_ptr, house),
                    Rcpp::_["latent"] = rcpp_snapshot_latent(mcmc_ptr, house),
                    Rcpp::_["infectivity_profiles"] = rcpp_snapshot_infec_profile(mcmc_ptr, house),
                    Rcpp::_["latent_log_lik"] = rcpp_snapshot_log_lik(mcmc_ptr, house, augm.scheme().list_latent_export),
                    Rcpp::_["observation_log_lik"] = rcpp_snapshot_log_lik(mcmc_ptr, house, augm.scheme().list_obs_member),
                    Rcpp::_["transmission_log_lik"] = hsp_to_rcpp(mcmc_ptr->data.state.transm_log_lik(house).value()));
            }
            catch (const housepi::TrustedException &e)
            {
                msg::stop(msg::join(e.what(), "for household " + msg::val(hh_id)));
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

                parallel::for_each(
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
                                    latent_cols[c++][row] = hsp_to_rcpp(augm.at(a, indiv));
                                // observation nodes: generated for infected members, NA otherwise
                                for (size_t a : augm.scheme().list_obs_member)
                                    latent_cols[c++][row] = indiv.is_infected_status()
                                                                ? hsp_to_rcpp(augm.at(a, indiv))
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
                                    const size_t a = augm.consult(augm.def()[surr.row].name);
                                    size_t r = row_offsets[H];
                                    for (const Individual &indiv : house)
                                        latent_cols[c][r++] = hsp_to_rcpp(augm.at(a, indiv));
                                }
                                ++c;
                            }
                        }
                        catch (const housepi::TrustedException &e)
                        {
                            msg::stop(msg::join(e.what(), "for " + data_wide.str(data_wide.at(H))));
                        }
                    },
                    engine.setup.budget.inner);
            }

            /// Record parameter values per repetition, one column per parameter.
            struct ParamTrace
            {
                ParamTrace(const Parameter &param, size_t n_rep)
                    : param_(param)
                {
                    for (const std::string &name : param.all())
                        draws_.add(name, n_rep);
                }

                /// Record the current parameter values as repetition `r`.
                void record(size_t r)
                {
                    for (size_t p = 0; p < draws_.cols.size(); ++p)
                        draws_.cols[p][r] = param_.at(p).value();
                }
                /// Flush the recorded values as a data.frame, releasing C++ storage.
                Rcpp::DataFrame flush() { return flush_draws(draws_); }

            private:
                const Parameter &param_; /// Parameters traced.
                Posterior::Draws draws_; /// Values per parameter, one row per repetition.
            };

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
                    sim_list[n_shared + 1 + c] = rcpp_typed_column(latent_cols[c], discr_cols[c]);
                return rcpp_data_frame(sim_list, sim_names);
            }

            /// Restore observed status and point every read-cell back at its imported value.
            ///
            /// - post-simulation reads of data_wide return the input;
            /// - generated values leave through the simulated frames;
            /// - wired cells refer to Variables the next run clears.
            void restore(MCMC &engine)
            {
                engine.augm.release_simulation();
                engine.data.wide().recover_all_status();
                engine.data.reset_readcells();
            }
        }

        //////  Progress bar  //////

        /// Report the progress of a run to the R session.
        ///
        /// - each report checks for user interrupts;
        /// - R `NULL` handle: no bar to fill.
        struct Progress
        {
            ProgressBar *bar; /// Bar filled, null for none.

            explicit Progress(SEXP handle)
                : bar(Rf_isNull(handle) ? nullptr : Rcpp::XPtr<ProgressBar>(handle).get()) {}

            /// Report `done` units of the current run.
            void operator()(double done) const
            {
                Rcpp::checkUserInterrupt();
                if (bar)
                    bar->update(done);
            }
            /// Add `units` filled by a completed run.
            void advance(double units) const
            {
                if (bar)
                    bar->advance(units);
            }
        };

        /// Open progress bar `name` over `total` units, printed every `pace` seconds.
        /// @param width Symbols in a full bar.
        /// @param symbol Character filling the bar.
        Rcpp::XPtr<ProgressBar> progress_open(std::string name, double total, double pace, unsigned width,
                                              bool is_interactive, std::string symbol)
        {
            return Rcpp::XPtr<ProgressBar>(
                new ProgressBar(name, total, pace, width, is_interactive, symbol), true);
        }

        /// Close progress bar `ptr` with its closing line.
        void progress_close(Rcpp::XPtr<ProgressBar> ptr) { ptr->close(); }
        /// End the line of progress bar `ptr` left open by an error.
        void progress_crash(Rcpp::XPtr<ProgressBar> ptr) { ptr->crash(); }

        /// Run the simulation, returning wide and long results after the reporting process.
        Rcpp::List run_simulation(Rcpp::XPtr<MCMC> mcmc_ptr,
                                  Rcpp::Nullable<Rcpp::Function> reporting_wide,
                                  Rcpp::Nullable<Rcpp::Function> apply_wide,
                                  Rcpp::Nullable<Rcpp::Function> reporting_and_apply_long,
                                  SEXP progress)
        {
            const Progress report(progress);
            MCMC &engine = *mcmc_ptr;
            Setup &setup = engine.setup;
            const Augmented &augm = engine.augm;
            const data::Wide &data_wide = engine.data.wide();

            // every exit restores the engine
            try
            {
                engine.augm.init_house_for_simulation(); // build household-level observations
                data_wide.assert_required_everyone(engine.data.state, augm.get_live_cols());
                engine.model->init_model_for_simulation(); // model-specific initialization
                engine.risks.assert_identifiability(engine.param);
                engine.augm.init_members_for_simulation(engine.status_process); // bind member definitions once per run

                // @todo Assert the infectivity profile here.
                // latents are not built yet, so a latent profile cannot be checked
                // engine.assert_valid_infectivity_profile(); // assert profile is valid

                //////  Reporting Process & Measures  //////

                // @todo Hold nullable R functions without a pointer work-around.
                std::unique_ptr<Rcpp::Function> f_reporting = nullptr, f_wide = nullptr, f_long = nullptr;

                if (!reporting_wide.isNull())
                    f_reporting = std::make_unique<Rcpp::Function>(reporting_wide);
                if (!apply_wide.isNull())
                    f_wide = std::make_unique<Rcpp::Function>(apply_wide);
                if (!reporting_and_apply_long.isNull())
                    f_long = std::make_unique<Rcpp::Function>(reporting_and_apply_long);

                //////  Initialization  //////

                Rcpp::List out_wide(setup.n_steps); // allocate memory for wide results
                Rcpp::List out_long;                // long results are optional
                if (f_long)
                    out_long = Rcpp::List(setup.n_steps);

                const std::vector<std::string> static_cols = engine.get_static_data_names();
                RcppFrame shared_frame(n_wide_cols(static_cols.size(), false));
                add_wide_columns(shared_frame, engine.data, data_wide, static_cols, false);
                const Rcpp::List shared_df = shared_frame.frame();

                // simulated data column names
                Rcpp::CharacterVector sim_names = Rcpp::wrap(engine.get_col_names());
                // row offsets for house parallelization
                const std::vector<size_t> row_offsets = engine.get_row_offsets();

                // flags for discrete columns, set once variables are bound after the first repetition
                std::vector<bool> discr_cols(augm.n_export_cols(), false);

                simulate::ParamTrace param_trace(engine.param, setup.n_steps);

                //////  Simulations  //////

                for (size_t r = 0; r < setup.n_steps; ++r)
                {
                    engine.augm.clear_all();
                    engine.data.wide().recover_all_status();
                    engine.param.load_sample_values(r);
                    engine.param.draw_sample(setup.prng()); // draw parameters with no pinned value
                    param_trace.record(r);

                    // this repetition's own columns: status, exported latents, observations and surrogates
                    Rcpp::IntegerVector status_col(data_wide.report_indiv.total);
                    std::vector<Rcpp::NumericVector> latent_cols = rcpp_columns(augm.n_export_cols(), data_wide.report_indiv.total);

                    // run one repetition
                    simulate::single_repetition(
                        engine, status_col, latent_cols, row_offsets);

                    if (r == 0)
                        discr_cols = augm.discrete_export_cols();

                    // build simulated list
                    Rcpp::DataFrame simulated(simulate::build_simulated(
                        shared_df, status_col, latent_cols, discr_cols, sim_names));

                    // wide reporting process
                    if (f_reporting)
                        simulated = (*f_reporting)(simulated);

                    // apply the measure to data-wide, or store it
                    if (f_wide)
                        out_wide[r] = (*f_wide)(simulated);
                    else
                        out_wide[r] = simulated;
                    // long reporting process and apply measure
                    if (f_long)
                        out_long[r] = (*f_long)(simulated);

                    report((r + 1.) / setup.n_steps);
                }
                report.advance(1.);
                simulate::restore(engine);

                return Rcpp::List::create(Rcpp::_["wide"] = out_wide,
                                          Rcpp::_["long"] = out_long,
                                          Rcpp::_["parameters"] = param_trace.flush());
            }
            catch (...)
            {
                simulate::restore(engine);
                throw;
            }
        }

        //////////////////////////////  Rcpp exports  //////////////////////////////

        /// Export the keys `expr` reads on the engine.
#define MAKE_EXPORT_KEYS(name, expr)                      \
    Rcpp::CharacterVector name(Rcpp::XPtr<MCMC> mcmc_ptr) \
    {                                                     \
        return Rcpp::wrap(mcmc_ptr->expr);                \
    }

        // keys the model declares through HOUSEPI_DECLARE_REGISTER_KEYS
        MAKE_EXPORT_KEYS(key_Data, model->key_Data())
        MAKE_EXPORT_KEYS(key_Long, model->key_Long())
        MAKE_EXPORT_KEYS(key_Param, model->key_Param())
        MAKE_EXPORT_KEYS(key_GroupRisk, model->key_GroupRisk())
        MAKE_EXPORT_KEYS(key_Latent, model->key_Latent())
        // keys the model requests as discrete
        MAKE_EXPORT_KEYS(key_Param_discrete, param.discrete_keys())
        MAKE_EXPORT_KEYS(key_Latent_discrete, augm.discrete_keys())
        MAKE_EXPORT_KEYS(key_Data_discrete, data.wide().discrete_keys())

        /// Baselines the model anchors to risks, keyed by parameter name.
        Rcpp::List risk_anchors(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            const std::vector<risks::Anchor> &anchors = mcmc_ptr->risks.anchors();

            Rcpp::List out(anchors.size());
            Rcpp::CharacterVector names(anchors.size());
            for (size_t a = 0; a < anchors.size(); ++a)
            {
                names[a] = mcmc_ptr->param.consult(anchors[a].baseline);
                out[a] = Rcpp::List::create(
                    Rcpp::_["risks"] = Rcpp::wrap(mcmc_ptr->risks.group().consult(anchors[a].risks)),
                    Rcpp::_["contacts"] = anchors[a].contacts);
            }
            out.names() = names;
            return out;
        }

        /// Group risks the model reads on members with infected status only.
        Rcpp::CharacterVector risk_status_stratified(Rcpp::XPtr<MCMC> mcmc_ptr)
        {
            const risks::Group &group = mcmc_ptr->risks.group();
            return Rcpp::wrap(group.consult(group.status_stratified()));
        }

        Rcpp::XPtr<MCMC> create_mcmc_runner(std::string module_name)
        {
            return Rcpp::XPtr<MCMC>(new MCMC(module_name), true);
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

        double evaluate_log_lik(Rcpp::XPtr<MCMC> mcmc_ptr, std::string hh_id)
        {
            return mcmc_ptr->evaluate_log_lik(hh_id);
        }
        /// Run inference on every engine of `engine_list`, `workers` of them at once.
        ///
        /// - main thread: checks user interrupts and refreshes `progress`;
        /// - first error or user interrupt: every engine stops at its next iteration.
        /// @param progress Progress handle filled one unit per engine, or R `NULL`.
        /// @return List of `seconds` per engine, the `error` that stopped the run and the engine that `failed`.
        Rcpp::List run_inference_list(Rcpp::List engine_list, int workers, SEXP progress)
        {
            const Progress report(progress);
            const size_t n = engine_list.size();
            std::vector<MCMC *> mcmc_vec(n);

            for (size_t i = 0; i < n; ++i)
                mcmc_vec[i] = Rcpp::XPtr<MCMC>(engine_list[i]).get();
            const InferenceRun run = run_inference(mcmc_vec, std::max(1, workers), report);
            report.advance(n);

            Rcpp::NumericVector seconds(n);
            for (size_t i = 0; i < n; ++i)
                seconds[i] = hsp_to_rcpp(run.seconds[i]);

            // the error as the R condition Rcpp builds, classed by exception type
            // an exception outside std::exception (user interrupt) escapes to R
            Rcpp::RObject error;
            int failed = NA_INTEGER;
            if (run.failure.error)
                try
                {
                    std::rethrow_exception(run.failure.error);
                }
                catch (const std::exception &e)
                {
                    error = ::exception_to_r_condition(e);
                    failed = static_cast<int>(run.failure.index) + 1;
                }

            return Rcpp::List::create(Rcpp::_["seconds"] = seconds, Rcpp::_["error"] = error,
                                      Rcpp::_["failed"] = failed);
        }

        /// Aggregate individuals by hh_id, given as 1-based codes of the distinct ids.
        /// - missing status: NA (NA_INTEGER is INT_MIN).
        measure::cases::Table summarise_houses_from_df(Rcpp::DataFrame data_wide)
        {
            Rcpp::CharacterVector col_names = data_wide.names();

            std::vector<int> hh_id = {};
            if (vect::contains(col_names, "hh_id"))
                hh_id = Rcpp::as<std::vector<int>>(data_wide["hh_id"]);

            std::vector<int> status = {};
            if (vect::contains(col_names, "status"))
                status = Rcpp::as<std::vector<int>>(data_wide["status"]);

            return measure::summarise_houses(hh_id, status);
        }

        /// Export `n` consecutive household sizes from `grid_min`.
        Rcpp::IntegerVector rcpp_hh_size_grid(unsigned grid_min, size_t n)
        {
            Rcpp::IntegerVector out(n);
            std::iota(out.begin(), out.end(), static_cast<int>(grid_min));
            return out;
        }

        /// Export a mean rate per household size, in column `name`.
        /// - size with no contributing household: NA.
        Rcpp::DataFrame rcpp_rate_by_hh_size(unsigned grid_min,
                                             const std::vector<double> &mean,
                                             const std::vector<unsigned> &count,
                                             const char *name)
        {
            Rcpp::NumericVector rate(mean.size());
            for (size_t i = 0; i < mean.size(); ++i)
                rate[i] = count[i] > 0 ? mean[i] : NA_REAL;

            RcppFrame frame(2);
            frame.add("hh_size", rcpp_hh_size_grid(grid_min, mean.size()));
            frame.add(name, rate);
            return frame.frame();
        }

        Rcpp::DataFrame count_houses_by_hh_size(Rcpp::DataFrame data_wide)
        {
            const measure::cases::HouseCounts res =
                measure::house_counts(summarise_houses_from_df(data_wide));

            RcppFrame frame(3);
            frame.add("hh_size", rcpp_hh_size_grid(res.grid_min, res.total.size()));
            frame.add("total", Rcpp::IntegerVector(res.total.begin(), res.total.end()));
            frame.add("at_least_1_case", Rcpp::IntegerVector(res.at_least_1_case.begin(), res.at_least_1_case.end()));
            return frame.frame();
        }

        Rcpp::DataFrame measure_ar_by_hh_size(Rcpp::DataFrame data_wide)
        {
            const measure::cases::AttackRates res =
                measure::attack_rate_by_hh_size(summarise_houses_from_df(data_wide));
            return rcpp_rate_by_hh_size(res.grid_min, res.ar_mean, res.count, "ar");
        }

        Rcpp::DataFrame measure_sar_by_hh_size(Rcpp::DataFrame data_wide)
        {
            const measure::cases::SecAttackRates res =
                measure::sar_by_hh_size(summarise_houses_from_df(data_wide));
            return rcpp_rate_by_hh_size(res.grid_min, res.sar_mean, res.count, "sar");
        }

        Rcpp::DataFrame measure_cases_by_hh_size(Rcpp::DataFrame data_wide)
        {
            const measure::cases::Table res = summarise_houses_from_df(data_wide);
            const size_t n_rows = res.counts.size();

            RcppFrame frame(2 + res.max_size);
            frame.add("hh_size", rcpp_hh_size_grid(res.min_size, n_rows));
            for (unsigned c = 0; c <= res.max_size; ++c)
            {
                // more cases than members: NA
                Rcpp::IntegerVector col(n_rows);
                for (size_t r = 0; r < n_rows; ++r)
                    col[r] = (c > res.min_size + r) ? NA_INTEGER : static_cast<int>(res.counts[r][c]);
                frame.add(std::to_string(c), col);
            }
            return frame.frame();
        }

        double measure_geom_mean_hh_size(Rcpp::DataFrame data_wide)
        {
            measure::cases::Table table = summarise_houses_from_df(data_wide);
            return table.max_cases > 1 ? measure::geom_mean_hh_size(table) : NA_REAL;
        }
    }
}

RCPP_MODULE(housepi)
{
    using namespace housepi;

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

    Rcpp::function("measure_sar_by_hh_size", &measure_sar_by_hh_size);
    Rcpp::function("measure_ar_by_hh_size", &measure_ar_by_hh_size);
    Rcpp::function("measure_cases_by_hh_size", &measure_cases_by_hh_size);
    Rcpp::function("measure_geom_mean_hh_size", &measure_geom_mean_hh_size);
    Rcpp::function("count_houses_by_hh_size", &count_houses_by_hh_size);

    Rcpp::function("aux_chain_binomial_debugger", &chain_binom::aux_debugger);

    Rcpp::function("create_mcmc_runner", &create_mcmc_runner);
    Rcpp::function("run_inference_list", &run_inference_list);
    Rcpp::function("run_simulation", &run_simulation);
    Rcpp::function("progress_open", &progress_open);
    Rcpp::function("progress_close", &progress_close);
    Rcpp::function("progress_crash", &progress_crash);
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
    Rcpp::function("import_status_process", &import_status_process);
    Rcpp::function("import_param_table", &import_param_table);
    Rcpp::function("has_contact_risk", &has_contact_risk);
    Rcpp::function("has_infec_profile", &has_infec_profile);
    Rcpp::function("export_data_wide", &export_data_wide);
    Rcpp::function("export_data_long", &export_data_long);
    Rcpp::function("summarise_house_cases", &summarise_house_cases);
    Rcpp::function("export_setup", &export_setup);
    Rcpp::function("export_augmented_data_wide", &export_augmented_data_wide);
    Rcpp::function("flush_param_draws", &flush_posterior_<&Posterior::param_values>);
    Rcpp::function("flush_param_log_prior_draws", &flush_posterior_<&Posterior::param_priors>);
    Rcpp::function("flush_latent_draws", &flush_posterior_<&Posterior::latent_values>);
    Rcpp::function("flush_house_log_lik_draws", &flush_posterior_<&Posterior::house_log_lik>);
    Rcpp::function("export_house_snapshot", &export_house_snapshot);
}
#endif
