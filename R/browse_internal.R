## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

rd_test_name_ <- paste("`character(1)` Name of the %1$s.",
                       "If `NA` print the list of available %1$ss.")
rd_test_args_ <- "named arguments. User argument values for the %s"
rd_test_discr_ <- "`logical(1)` Treat definition as discrete."
rd_test_return_ <- paste("`NULL` if the available names were printed.",
                         "Otherwise [%1$s][housepi::%1$s] object.")
rd_test_x_ <- "`numeric()` Points of evaluation."

# Show registered names of family, as `show_model_names()` does for templates.
# No `{.run}` links: there is no per-name description to open.
show_registered_ <- function(names, family) {
  fmt_names <- vapply(sort(names),
                      \(nm) inline_cli_("{.val {nm}}"),
                      character(1), USE.NAMES = FALSE)
  cli::cli_inform(c(i = "Available {family}s:",
                    paste(cli::symbol$bullet, fmt_names)))
  invisible()
}

# Interval notation of a `min`, `max`, `closed` list shifted by `loc`,
# e.g. "[0, Inf)".
format_interval_ <- function(x, loc = 0) {
  sprintf("%s%s, %s%s",
          if (x$closed[["min"]]) "[" else "(",
          format(x$min + loc),
          format(x$max + loc),
          if (x$closed[["max"]]) "]" else ")")
}

############################################################################
############################################################################
############################# Housepi Exposed ##############################
############################################################################
############################################################################

# dictionary of exposed object
.family <- list(func = "function",
                distr = "distribution",
                infec = "infectivity profile")

#' @title Exposed housepi definition
#' @description
#' Base of the objects returned by [browse_functions()],
#' [browse_distributions()] and [browse_infectivity_profiles()].
#' It carries one definition built by the C++ kernel, reports the user
#' argument values it was built with, and the domain it is defined on.
#' @examples
#' dst <- browse_distributions("gamma", shape = 2, scale = 1)
#' dst$args
#' dst$domain
#' @docType class
HousepiExposed <- R6::R6Class(
  classname = "HousepiExposed",
  cloneable = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @field
    obj_ = NULL,
    # family word, set by each derived class
    family_ = NA_character_,
    name_ = NA_character_,
    discr_ = FALSE,
    loc_ = NA_real_,

    # Call a C++ member through the shared interceptor.
    call_ = function(method, ...) {
      with_cpp_conditions(cpp_member_(private$obj_, method)(...))
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field args `numeric()`
    #' Named argument values the definition was built with.
    args = function(value) {
      if (!missing(value)) stop_read_only("args")
      private$call_("args")
    },

    #' @field domain `list()`
    #' Values the definition is evaluated on: bounds `min` and `max`, and
    #' `closed`, whether each bound belongs to the interval. The bounds are
    #' read before the location shift `loc`.
    domain = function(value) {
      if (!missing(value)) stop_read_only("domain")
      private$call_("domain")
    },

    #' @field is_discrete `logical(1)`
    #' Whether the values the definition takes lie on a lattice of
    #' integers. A distribution answers for the values it draws.
    is_discrete = function(value) {
      if (!missing(value)) stop_read_only("is_discrete")
      private$discr_
    },

    #' @field is_positive `logical(1)`
    #' Whether the values the definition takes are known to be strictly
    #' positive. A distribution answers for the values it draws.
    #' Note that `FALSE` reads as "not known to be positive".
    is_positive = function(value) {
      if (!missing(value)) stop_read_only("is_positive")
      private$call_("is_positive")
    },

    #' @field codomain `list()`
    #' Values evaluating the definition returns, in the same form as `domain`.
    #' A distribution returns a density, a mass or a cumulative probability,
    #' never a draw: see `domain` for the values it is defined on.
    codomain = function(value) {
      if (!missing(value)) stop_read_only("codomain")
      private$call_("codomain")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Internal: the `browse_*()` functions build this object,
    #' don't call `$new()` directly.
    #' @param obj Internal engine handle.
    #' @param name Registered name of the definition.
    initialize = function(obj, name) {
      private$obj_ <- obj
      private$name_ <- name
      private$discr_ <- private$call_("is_discrete")
      private$loc_ <- private$call_("loc")
      invisible(self)
    },

    #' @description Override the default printing behavior.
    print = function() {
      args <- self$args
      items <- vapply(names(args),
                      \(arg) inline_cli_("{.arg {arg}} = {args[[arg]]}"),
                      character(1), USE.NAMES = FALSE)
      if (!is.na(private$loc_)) {
        items <- c(items, "{.arg loc} = {format(private$loc_)}")
      }
      loc <- if (is.na(private$loc_)) 0 else private$loc_
      print_bullets(c(
        i = "Exposed {private$family_} {.val {private$name_}}.",
        items,
        paste(if (private$discr_) "Discrete" else "Continuous",
              "on {format_interval_(self$domain, loc)}.")
      ))
      invisible(self)
    }
  )
)

############################################################################
############################################################################
######################### Housepi Exposed Function #########################
############################################################################
############################################################################

#' @title Exposed housepi function
#' @description
#' A housepi function and its primitive on user-supplied points.
#' @examples
#' fnc <- browse_functions("pow", x2 = 3)
#' fnc$func(seq(0, 5))
#' fnc$prim(seq(0, 5))
#' @docType class
HousepiExposedFunction <- R6::R6Class(
  classname = "HousepiExposedFunction",
  inherit = HousepiExposed,
  cloneable = FALSE,

  private = list(family_ = .family$func),

  public = list(
    #' @description Evaluate the function.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    func = function(x) private$call_("func", x),

    #' @description Evaluate the natural logarithm of the function.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    log_func = function(x) private$call_("log_func", x),

    #' @description Evaluate the primitive of the function.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    prim = function(x) private$call_("prim", x)
  )
)

############################################################################
############################################################################
####################### Housepi Exposed Distribution #######################
############################################################################
############################################################################

#' @title Exposed housepi distribution
#' @description
#' A housepi distribution, its density/mass, cumulative, quantiles,
#' and random generation.
#' @examples
#' dst <- browse_distributions("gamma", shape = 2, scale = 1, max = 10)
#' dst$pdf(seq(0, 5))
#' dst$quantile(c(0.05, 0.5, 0.95))
#' dst$rand(5)
#' @docType class
HousepiExposedDistribution <- R6::R6Class(
  classname = "HousepiExposedDistribution",
  inherit = HousepiExposed,
  cloneable = FALSE,

  private = list(family_ = .family$distr),

  public = list(
    #' @description Evaluate the probability density/mass.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    pdf = function(x) private$call_("pdf", x),

    #' @description
    #' Evaluate the natural logarithm of the probability density/mass.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    log_pdf = function(x) private$call_("log_pdf", x),

    #' @description
    #' Evaluate the natural logarithm of the ratio between the probability
    #' density at two points.
    #' @param xo `numeric()` Points the ratio is taken from.
    #' @param xn `numeric()` Points the ratio is taken to.
    #' @return `numeric()`
    log_ratio = function(xo, xn) private$call_("log_ratio", xo, xn),

    #' @description Evaluate the cumulative distribution.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    cdf = function(x) private$call_("cdf", x),

    #' @description Evaluate the quantiles.
    #' @param p `numeric()` Probabilities, `NA` outside \eqn{[0, 1]}.
    #' @return `numeric()`
    quantile = function(p) private$call_("quantile", p),

    #' @description Draw from the distribution.
    #' @param n `integer(1)` Number of draws.
    #' @return `integer()` for a discrete distribution, `numeric()` otherwise.
    rand = function(n) private$call_("rand", n)
  )
)

############################################################################
############################################################################
######################### Housepi Exposed Profile ##########################
############################################################################
############################################################################

#' @title Exposed housepi infectivity profile
#' @description
#' A housepi infectivity profile, its instantaneous and cumulative profile.
#' @examples
#' ipf <- browse_infectivity_profiles("gamma", mean = 4, sd = 1.6)
#' ipf$inst(seq(0, 10))
#' ipf$cumul(seq(0, 10))
#' @docType class
HousepiExposedProfile <- R6::R6Class(
  classname = "HousepiExposedProfile",
  inherit = HousepiExposed,
  cloneable = FALSE,

  private = list(family_ = .family$infec),

  public = list(
    #' @description Evaluate the instantaneous infectivity profile.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    inst = function(x) private$call_("inst", x),

    #' @description Evaluate the cumulative infectivity profile.
    #' @param x `r rd_test_x_`
    #' @return `numeric()`
    cumul = function(x) private$call_("cumul", x)
  )
)

############################################################################
############################################################################
############################################################################

#' @title Navigate through housepi functions
#' @description
#' Print the list of available functions,
#' or expose one to access user defined functions.
#' @param name `r sprintf(rd_test_name_, .family$func)`
#' @param ... `r sprintf(rd_test_args_, .family$func)`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `r sprintf(rd_test_return_, "HousepiExposedFunction")`
#' @examples
#' # list all available functions
#' browse_functions()
#'
#' fnc <- browse_functions("pow", x2 = 3)
#' fnc$func(seq(0, 5))
#' fnc$prim(seq(0, 5))
#' @export
browse_functions <- function(name = NA_character_, ...,
                             min = NA_real_,
                             max = NA_real_,
                             discr = FALSE) {
  if (is.na(name)) {
    return(show_registered_(.engine$list_functions(), .family$func))
  }
  assert_unique_arg_names(names(list(...)), .family$func, name)
  obj <- .engine$expose_function(name, unlist(list(...)), min, max, discr)
  HousepiExposedFunction$new(obj, name)
}

#' @title Navigate through housepi distributions
#' @description
#' Print the list of available distributions,
#' or expose one to access user defined distributions.
#' @param name `r sprintf(rd_test_name_, .family$distr)`
#' @param ... `r sprintf(rd_test_args_, .family$distr)`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `r sprintf(rd_test_return_, "HousepiExposedDistribution")`
#' @examples
#' # list all available distributions
#' browse_distributions()
#'
#' dst <- browse_distributions("weibull", shape = 2, scale = 1)
#' dst$pdf(seq(0, 5))
#' dst$cdf(seq(0, 5))
#' dst$domain
#' @export
browse_distributions <- function(name = NA_character_, ...,
                                 min = NA_real_,
                                 max = NA_real_,
                                 discr = FALSE) {
  if (is.na(name)) {
    return(show_registered_(.engine$list_distributions(), .family$distr))
  }
  assert_unique_arg_names(names(list(...)), .family$distr, name)
  obj <- .engine$expose_distribution(name, unlist(list(...)), min, max, discr)
  HousepiExposedDistribution$new(obj, name)
}

#' @title Navigate through housepi infectivity profiles
#' @description
#' Print the list of available infectivity profiles,
#' or expose one to access user defined infectivity profiles.
#' @param name `r sprintf(rd_test_name_, .family$infec)`
#' @param ... `r sprintf(rd_test_args_, .family$infec)`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `r sprintf(rd_test_return_, "HousepiExposedProfile")`
#' @examples
#' # list all available infectivity profiles
#' browse_infectivity_profiles()
#'
#' ipf <- browse_infectivity_profiles("gamma", shape = 2, scale = 1)
#' ipf$inst(seq(0, 5))
#' ipf$cumul(seq(0, 5))
#' @export
browse_infectivity_profiles <- function(name = NA_character_, ...,
                                        min = NA_real_,
                                        max = NA_real_,
                                        discr = FALSE) {
  if (is.na(name)) {
    return(show_registered_(.engine$list_infec_profiles(), .family$infec))
  }
  # FIXME: this is unaware of loc !
  if (is.na(min)) min <- if (discr) 1 else 0
  assert_unique_arg_names(names(list(...)), .family$infec, name)
  obj <- .engine$expose_infec_profile(name, unlist(list(...)), min, max, discr)
  HousepiExposedProfile$new(obj, name)
}
