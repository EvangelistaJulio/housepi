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
rd_test_return_ <- paste(
  "An exposed %s object with methods `$%s()` and `$%s()`,",
  "or `NULL` if the available names were printed."
)

# Show the registered names of a family, as `show_model_names()` does for
# templates. No `{.run}` links: there is no per-name description to open.
show_registered_names_ <- function(names, family) {
  cli::cli_inform(c(
    stats::setNames(sprintf("Available %s:", family), "i"),
    paste0(cli::symbol$bullet, " ", sprintf("{.val %s}", sort(names)))
  ))
  invisible()
}

# Proxy over an exposed C++ object (function / distribution / infectivity
# profile). Its `$` routes every method call through the shared condition
# interceptor (reusing `wrap_cpp_member_`), and it prints legibly. Construction
# and `init()` stay wrapped explicitly in the `browse_*()` helpers, since the
# object does not exist yet at that point.
housepi_exposed_ <- function(obj, family, name) {
  structure(list(obj = obj, family = family, name = name),
            class = "housepi_exposed")
}

#' @noRd
`$.housepi_exposed` <- function(x, name) {
  wrap_cpp_member_(.subset2(x, "obj"), name)
}

#' @noRd
print.housepi_exposed <- function(x, ...) {
  args <- x$args()
  header <- c(i = sprintf("Exposed %s {.val %s}.",
                          .subset2(x, "family"), .subset2(x, "name")))
  items <- if (length(args) == 0) character() else
    sprintf("{.arg %s} = %s", names(args), format(unname(args)))
  print_bullets(c(header, items))
  invisible(x)
}

#' @title Navigate through housepi functions
#' @description
#' Print the list of available functions,
#' or expose one to access user defined functions.
#' @param name `r sprintf(rd_test_name_, "function")`
#' @param ... `r sprintf(rd_test_args_, "function")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `r sprintf(rd_test_return_, "function", "func", "prim")`
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
  if (is.na(name)) return(show_registered_names_(.engine$list_functions(),
                                                 "functions"))
  # build routes through the `.engine` proxy, so construction/init errors are
  # already intercepted; the proxy below only adds method interception + print.
  obj <- .engine$expose_function(name, unlist(list(...)), min, max, discr)
  housepi_exposed_(obj, "function", name)
}

#' @title Navigate through housepi distributions
#' @description
#' Print the list of available distributions,
#' or expose one to access user defined distributions.
#' @param name `r sprintf(rd_test_name_, "distribution")`
#' @param ... `r sprintf(rd_test_args_, "distribution")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `r sprintf(rd_test_return_, "distribution", "pdf", "cdf")`
#' @examples
#' # list all available distributions
#' browse_distributions()
#'
#' dst <- browse_distributions("weibull", shape = 2, scale = 1)
#' dst$pdf(seq(0, 5))
#' dst$cdf(seq(0, 5))
#' @export
browse_distributions <- function(name = NA_character_, ...,
                                 min = NA_real_,
                                 max = NA_real_,
                                 discr = FALSE) {
  if (is.na(name)) return(show_registered_names_(.engine$list_distributions(),
                                                 "distributions"))
  obj <- .engine$expose_distribution(name, unlist(list(...)), min, max, discr)
  housepi_exposed_(obj, "distribution", name)
}

#' @title Navigate through housepi infectivity profiles
#' @description
#' Print the list of available infectivity profiles,
#' or expose one to access user defined infectivity profiles.
#' @param name `r sprintf(rd_test_name_, "infectivity profile")`
#' @param ... `r sprintf(rd_test_args_, "infectivity profile")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `r sprintf(rd_test_return_, "infectivity profile", "inst", "cumul")`
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
  if (is.na(name)) return(show_registered_names_(.engine$list_infec_profiles(),
                                                 "infectivity profiles"))
  if (is.na(min)) min <- ifelse(discr, 1, 0)
  obj <- .engine$expose_infec_profile(name, unlist(list(...)), min, max, discr)
  housepi_exposed_(obj, "infectivity profile", name)
}
