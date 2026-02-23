/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifdef _HOUSEPI_RCPP

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"

#include "distribution.h"
#include "function.h"

#include <Rcpp.h>

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
        if (is_new_housepi_warning(txt))
        {
            Rf_warningcall(R_NilValue, "%s", txt.c_str());
            add_new_housepi_warning(txt);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  *** Housepi Rcpp Base : class ***  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//' @name ExposeFunction
//' @title Rcpp Function from Housepi
//' @description Call housepi registered function.
//' @field initialize \itemize{
//'   \item Parameter: name `character(1)` with function name.
//'   \item Parameter: args A named `numeric()` with arguments.
//'   \item Parameter: min `numeric(1)` Minimum value of the function.
//'   \item Parameter: max `numeric(1)` Maximum value of the function.
//'   \item Parameter: is_discr `logical(1)` Wheter discretized.
//' }
//' @field func \itemize{
//'   \item Parameter: x `numeric()` Input vector.
//' }
//' @return `numeric()`
//' @field prim \itemize{
//'   \item Parameter: x `numeric()` Input vector.
//' }
//' @return `numeric()`
//' @keywords internal

//' @name ExposeDistribution
//' @title Rcpp Distribution from Housepi
//' @description Call housepi registered distribution.
//' @field initialize \itemize{
//'   \item Parameter: name `character(1)` Name of the distribution.
//'   \item Parameter: args `numeric()` Custom arguments.
//'   \item Parameter: min `numeric(1)` Minimum value of the distribution.
//'   \item Parameter: max `numeric(1)` Maximum value of the distribution.
//'   \item Parameter: is_discr `logical(1)` Wheter discretized.
//' }
//' @field pdf \itemize{
//'   \item Parameter: x `numeric()` Input vector.
//' }
//' @return `numeric()`
//' @field cdf \itemize{
//'   \item Parameter: x `numeric()` Input vector.
//' }
//' @return `numeric()`
//' @keywords internal

//' @name ExposeInfecProfile
//' @title Rcpp Infectivity Profile from Housepi
//' @description Call housepi registered infectivity profile.
//' @field initialize \itemize{
//'   \item Parameter: name `character(1)` Name of the profile.
//'   \item Parameter: args `numeric()` Custom arguments.
//'   \item Parameter: min `numeric(1)` Minimum value of the profile.
//'   \item Parameter: max `numeric(1)` Maximum value of the profile.
//'   \item Parameter: is_discr `logical(1)` Wheter discretized.
//' }
//' @field inst \itemize{
//'   \item Parameter: x `numeric()` Input vector.
//' }
//' @return `numeric()`
//' @field cumul \itemize{
//'   \item Parameter: x `numeric()` Input vector.
//' }
//' @return `numeric()`
//' @keywords internal

namespace
{
    // Base class to wrap internal hosuepi functions
    class ExposeBase
    {
    public:
        ExposeBase(std::string name, const std::string &family,
                   SEXP args_, double min, double max, bool is_discr,
                   std::vector<std::string> registry_types)
            : family_(family)
        {
            try
            {
                ptr_ = compose_func_registries(name, family_, registry_types, false);

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

                if (!hsp::na::is(min))
                    ptr_->init_attr("min", Attribute(min));
                if (!hsp::na::is(max))
                    ptr_->init_attr("max", Attribute(max));

                if (is_discr)
                    ptr_->init_discretization();
                ptr_->init_truncation();
                ptr_->assert_is_callable_explicitly();
            }
            catch (const std::exception &e)
            {
                hsp::msg::error(
                    "While initializing " + family_ + " '" + name + "':\n" + e.what());
            }
        }
        virtual ~ExposeBase() = default;

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
                hsp::msg::error(
                    "While calling " + family_ + " '" + std::string(ptr_->name()) + "':\n" + e.what());
            }
        }
    };

    class ExposeFunction : public ExposeBase
    {
    public:
        ExposeFunction(std::string name_, SEXP args_,
                       double min, double max, bool is_discr)
            : ExposeBase(name_, "function", args_, min, max, is_discr,
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
        ExposeDistribution(std::string name_, SEXP args_,
                           double min, double max, bool is_discr)
            : ExposeBase(name_, "distribution", args_, min, max, is_discr,
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
        ExposeInfecProfile(std::string name_, SEXP args_,
                           double min, double max, bool is_discr)
            : ExposeBase(name_, "infectivity profile", args_, min, max, is_discr,
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
}

RCPP_MODULE(housepi)
{
    Rcpp::class_<ExposeBase>("ExposeBase")
        .method("args", &ExposeBase::args);

    Rcpp::class_<ExposeFunction>("ExposeFunction")
        .constructor<std::string, SEXP, double, double, bool>()
        .derives<ExposeBase>("ExposeBase")
        .method("func", &ExposeFunction::func)
        .method("prim", &ExposeFunction::prim);

    Rcpp::class_<ExposeDistribution>("ExposeDistribution")
        .constructor<std::string, SEXP, double, double, bool>()
        .derives<ExposeBase>("ExposeBase")
        .method("pdf", &ExposeDistribution::pdf)
        .method("cdf", &ExposeDistribution::cdf);

    Rcpp::class_<ExposeInfecProfile>("ExposeInfecProfile")
        .constructor<std::string, SEXP, double, double, bool>()
        .derives<ExposeBase>("ExposeBase")
        .method("inst", &ExposeInfecProfile::inst)
        .method("cumul", &ExposeInfecProfile::cumul);
}
#endif
