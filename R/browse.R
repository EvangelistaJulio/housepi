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
# No `{.run}` links, since no name has a description to open.
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

########  HousepiExposed  ########

# family words of exposed objects
.family <- list(func = "function",
                distr = "distribution",
                infec = "infectivity profile")

#' @title Exposed housepi definition
#' @description
#' Base of the objects returned by [browse_functions()],
#' [browse_distributions()] and [browse_infectivity_profiles()].
#' It holds one definition, the argument values it was built with,
#' and the domain it is defined on.
#' @examples
#' dst <- browse_distributions("gamma", shape = 2, scale = 1)
#' dst$args
#' dst$domain
#' @docType class
HousepiExposed <- R6::R6Class(
  classname = "HousepiExposed",
  cloneable = FALSE,

  ####  Private  ####

  private = list(
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

  ####  Active  ####

  active = list(
    #' @field args `numeric()`
    #' Named argument values the definition was built with.
    args = function(value) {
      assert_read_only_(value, "args")
      private$call_("args")
    },

    #' @field domain `list()`
    #' Values the definition is evaluated on: bounds `min` and `max`,
    #' and `closed`, whether each bound belongs to the interval.
    #' The bounds are read before the location shift `loc`.
    domain = function(value) {
      assert_read_only_(value, "domain")
      private$call_("domain")
    },

    #' @field is_discrete `logical(1)`
    #' Whether the values the definition takes lie on a lattice of integers.
    #' A distribution answers for the values it draws.
    is_discrete = function(value) {
      assert_read_only_(value, "is_discrete")
      private$discr_
    },

    #' @field is_positive `logical(1)`
    #' Whether the values the definition takes are known to be above zero.
    #' A distribution answers for the values it draws.
    #' `FALSE` means not known to be above zero.
    is_positive = function(value) {
      assert_read_only_(value, "is_positive")
      private$call_("is_positive")
    },

    #' @field codomain `list()`
    #' Values evaluating the definition returns, in the same form as `domain`.
    #' A distribution returns a density, a mass or a cumulative probability,
    #' never a draw: see `domain` for the values it is defined on.
    codomain = function(value) {
      assert_read_only_(value, "codomain")
      private$call_("codomain")
    }
  ),

  ####  Public  ####

  public = list(
    #' @description
    #' Create an exposed definition.
    #' Use the `browse_*()` functions instead of calling `$new()`.
    #' @param obj Internal engine handle.
    #' @param name Registered name of the definition.
    initialize = function(obj, name) {
      private$obj_ <- obj
      private$name_ <- name
      private$discr_ <- private$call_("is_discrete")
      private$loc_ <- private$call_("loc")
      invisible(self)
    },

    #' @description Print the family, name, arguments and domain.
    print = function() {
      args <- vapply(self$args, \(x) escape_cli_(format(x)), character(1))
      if (!is.na(private$loc_)) args <- c(args, loc = format(private$loc_))
      loc <- if (is.na(private$loc_)) 0 else private$loc_
      print_rows_("Exposed {private$family_} {.val {private$name_}}.", c(
        args,
        Domain = paste(if (private$discr_) "discrete" else "continuous",
                       "on {format_interval_(self$domain, loc)}")
      ))
      invisible(self)
    }
  )
)

########  HousepiExposedFunction  ########

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

########  HousepiExposedDistribution  ########

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
    #' Evaluate the log-ratio of the density/mass between two points.
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
    #' @return `numeric()` Draws, `integer()` for a discrete distribution.
    rand = function(n) private$call_("rand", n)
  )
)

########  HousepiExposedProfile  ########

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

########  Browse  ########

#' @title Browse housepi functions
#' @description
#' Print the list of available functions,
#' or expose one to evaluate it.
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
  if (is.na(name))
    return(show_registered_(.engine$list_functions(), .family$func))

  assert_unique_arg_names(names(list(...)), .family$func, name)
  obj <- .engine$expose_function(name, unlist(list(...)), min, max, discr)
  HousepiExposedFunction$new(obj, name)
}

#' @title Browse housepi distributions
#' @description
#' Print the list of available distributions,
#' or expose one to evaluate it.
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
  if (is.na(name))
    return(show_registered_(.engine$list_distributions(), .family$distr))

  assert_unique_arg_names(names(list(...)), .family$distr, name)
  obj <- .engine$expose_distribution(name, unlist(list(...)), min, max, discr)
  HousepiExposedDistribution$new(obj, name)
}

#' @title Browse housepi infectivity profiles
#' @description
#' Print the list of available infectivity profiles,
#' or expose one to evaluate it.
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
  if (is.na(name))
    return(show_registered_(.engine$list_infec_profiles(), .family$infec))

  # @todo the default `min` ignores `loc`
  if (is.na(min)) min <- if (discr) 1 else 0
  assert_unique_arg_names(names(list(...)), .family$infec, name)
  obj <- .engine$expose_infec_profile(name, unlist(list(...)), min, max, discr)
  HousepiExposedProfile$new(obj, name)
}

# Directory of all housepi template models.
template_models_dir <- function() {
  system.file("models", package = "housepi")
}

#' Directory of template model `x`.
#' @param x `character(1)` Template name.
#' @noRd
get_template_dir <- function(x) {
  flat_choice <- spaced_choice(x)
  path_to_options <- list.files(template_models_dir(), full.names = TRUE)
  flat_options <- spaced_choice(basename(path_to_options))
  # extract path if present (unique since folder name)
  mask <- flat_options == flat_choice
  if (any(mask)) return(path_to_options[mask])
  NA_character_
}

# Hint pointing to the template listing.
browse_templates_hint_ <- paste(
  "Show the available templates with",
  "{.run [browse_templates()](housepi::browse_templates())}."
)

#' Suggest the closest template name to a mistyped `x`, as "did you mean".
#' @noRd
find_closest_template_name <- function(x) {
  flat_choice <- spaced_choice(x)
  options <- list.files(template_models_dir(), full.names = FALSE)
  flat_options <- spaced_choice(options)
  # no suggestion for trivial mistakes
  sms <- find_closest_match(flat_choice, flat_options)
  if (length(sms) && nchar(flat_choice) > 3) return(sms)
  browse_templates_hint_
}

#' Show the names of the available template models.
#' @noRd
show_model_names <- function() {
  model_names <- list.files(template_models_dir())
  links <- vapply(
    tolower(gsub("-| |_", "_", model_names)),
    \(m) inline_cli_("{.run [{m}](housepi::browse_templates(\"{m}\"))}"),
    character(1), USE.NAMES = FALSE
  )
  cli::cli_inform(c(i = "Available template models:",
                    paste(cli::symbol$bullet, links)))
}

#' Copy template model `name` to directory `path`.
#' @noRd
copy_template_to_path <- function(name, path) {
  # avoid current directory
  if (path == ".") {
    path <- name
    warn_h("Using current directory is not allowed for copying a template.",
           `!` = "Creating folder {.file {path}} instead.")
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) raise_permissions_error()

  # source cpp is always present, definition templates are optional
  for (f in c("model.cpp", unname(template_def_files))) {
    from_file <- file.path(get_template_dir(name), f)
    if (file.exists(from_file)) overwrite(from_file, file.path(path, f))
  }
}

#' Open the help page of template model `name`.
#' @noRd
open_model_helper_ <- function(name) {
  print(utils::help(paste0("template_", name), package = "housepi"))
  invisible()
}

#' @title Browse household model templates
#' @description Print the list of available templates,
#' get their descriptions or copy their source code to a custom path.
#' @param name `character(1)` Template name.
#' If `NA`, print the list of available templates.
#' @param path `character(1)` Path to copy the template scripts to.
#' If `NA` and `name` is not, open the template help page.
#' @param replace `logical(1)` Whether to replace the content of `path`.
#' @return `NULL`, invisibly.
#' The function prints the templates, opens a description or copies the files.
#' @examples
#' # list all available templates
#' browse_templates()
#' @export
browse_templates <- function(name = NA_character_,
                             path = NA_character_,
                             replace = FALSE) {
  if (is.na(name)) show_model_names()
  else if (!is.character(name)) stop_h("Model name is not a character.")
  else if (is.na(get_template_dir(name))) stop_h(
    "Model {.val {name}} was not a known template.",
    find_closest_template_name(name)
  )
  # asking for template model description
  else if (is.na(path)) open_model_helper_(basename(get_template_dir(name)))
  # asking for copying template to folder with model
  else if (file.exists(file.path(path, "model.cpp")) && !replace) warn_h(
    "A model already exists in the specified folder.",
    "To replace it, run again with option {.code replace = TRUE}."
  )
  # asking for copying template to empty folder
  else copy_template_to_path(name, path)
  invisible()
}
