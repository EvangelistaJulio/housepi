## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

rd_test_name_ <- paste("`character(1)` Name of the %1$s.",
                       "If `NA` print the list of available %1$ss.")
rd_test_args_ <- "named arguments. User argument values for the %s"
rd_test_discr_ <- "`logical(1)` Treat definition as discrete."
rd_test_return_ <- paste("`NULL` if the available names were printed.",
                         "Otherwise [%1$s][housepi::%1$s] object.")
rd_test_x_ <- "`numeric()` Points of evaluation."

# Show the registered names of a family, as `show_model_names()` does for
# templates. No `{.run}` links: there is no per-name description to open.
show_registered_ <- function(names, family) {
  cli::cli_inform(c(
    stats::setNames(sprintf("Available %ss:", family), "i"),
    paste0(cli::symbol$bullet, " ", sprintf("{.val %s}", sort(names)))
  ))
  invisible()
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

    #' @field range `character(1)`
    #' Domain of the definition, in interval notation.
    range = function(value) {
      if (!missing(value)) stop_read_only("range")
      private$call_("range")
    },

    #' @field is_discrete `logical(1)`
    #' Whether the definition is evaluated on a discrete domain.
    is_discrete = function(value) {
      if (!missing(value)) stop_read_only("is_discrete")
      private$discr_
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' New object of class [HousepiExposed][housepi::HousepiExposed].
    #' @param obj Definition built by the C++ kernel.
    #' @param name `character(1)` Registered name of the definition.
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
      items <- if (length(args) == 0) character() else
        sprintf("{.arg %s} = %s", names(args), format(unname(args)))
      if (!is.na(private$loc_))
        items <- c(items, sprintf("{.arg loc} = %s", format(private$loc_)))
      print_bullets(c(
        i = sprintf("Exposed %s {.val %s}.", private$family_, private$name_),
        items,
        sprintf("%s on %s.", ifelse(private$discr_, "Discrete", "Continuous"),
                self$range)
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

    #' @description
    #' Evaluate the natural logarithm of the ratio between the function at two
    #' points.
    #' @param xo `numeric()` Points the ratio is taken from.
    #' @param xn `numeric()` Points the ratio is taken to.
    #' @return `numeric()`
    log_ratio = function(xo, xn) private$call_("log_ratio", xo, xn),

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
#' dst$range
#' @export
browse_distributions <- function(name = NA_character_, ...,
                                 min = NA_real_,
                                 max = NA_real_,
                                 discr = FALSE) {
  if (is.na(name)) {
    return(show_registered_(.engine$list_distributions(), .family$distr))
  }
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
  if (is.na(min)) min <- ifelse(discr, 1, 0)
  obj <- .engine$expose_infec_profile(name, unlist(list(...)), min, max, discr)
  HousepiExposedProfile$new(obj, name)
}
