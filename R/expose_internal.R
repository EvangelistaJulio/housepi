## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

rd_test_name_ <- "`character(1)` Name of the %s."
rd_test_x_ <- "`numeric()` Points where to evaluate."
rd_test_args_ <- "named arguments. User argument values for the %s"
rd_test_discr_ <- "`logical(1)` Treat definition as discrete."

# Proxy over an exposed C++ object (function / distribution / infectivity
# profile). Its `$` routes every method call through the shared condition
# interceptor (reusing `wrap_cpp_member_`), and it prints legibly. Construction
# and `init()` stay wrapped explicitly in the `expose_*()` helpers, since the
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

#' @title Expose housepi function
#' @description
#' Expose an object from housepi registered function.
#' Used to access user defined functions.
#' @param name `r sprintf(rd_test_name_, "function")`
#' @param ... `r sprintf(rd_test_args_, "function")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return An exposed function object with methods `$func()` and `$prim()`.
#' @examples
#' fnc <- expose_function("pow", x2 = 3)
#' fnc$func(seq(0, 5))
#' fnc$prim(seq(0, 5))
#' @export
expose_function <- function(name, ...,
                            min = NA_real_,
                            max = NA_real_,
                            discr = FALSE) {
  # build routes through the `.engine` proxy, so construction/init errors are
  # already intercepted; the proxy below only adds method interception + print.
  obj <- .engine$expose_function(name, unlist(list(...)), min, max, discr)
  housepi_exposed_(obj, "function", name)
}

#' @title Expose housepi distribution
#' @description
#' Expose an object from housepi registered distribution.
#' Used to access user defined distributions.
#' @param name `r sprintf(rd_test_name_, "distribution")`
#' @param ... `r sprintf(rd_test_args_, "distribution")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return An exposed distribution object with methods `$pdf()` and `$cdf()`.
#' @examples
#' dst <- expose_distribution("weibull", shape = 2, scale = 1)
#' dst$pdf(seq(0, 5))
#' dst$cdf(seq(0, 5))
#' @export
expose_distribution <- function(name, ...,
                                min = NA_real_,
                                max = NA_real_,
                                discr = FALSE) {
  obj <- .engine$expose_distribution(name, unlist(list(...)), min, max, discr)
  housepi_exposed_(obj, "distribution", name)
}

#' @title Expose housepi infectivity profile
#' @description
#' Expose an object from housepi registered infectivity profile.
#' Used to access user defined infectivity profiles.
#' @param name `r sprintf(rd_test_name_, "infectivity profile")`
#' @param ... `r sprintf(rd_test_args_, "infectivity profile")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return An exposed infectivity profile with methods `$inst()` and `$cumul()`.
#' @examples
#' ipf <- expose_infectivity_profile("gamma", shape = 2, scale = 1)
#' ipf$inst(seq(0, 5))
#' ipf$cumul(seq(0, 5))
#' @export
expose_infectivity_profile <- function(name, ...,
                                       min = NA_real_,
                                       max = NA_real_,
                                       discr = FALSE) {
  if (is.na(min)) min <- ifelse(discr, 1, 0)
  obj <- .engine$expose_infec_profile(name, unlist(list(...)), min, max, discr)
  housepi_exposed_(obj, "infectivity profile", name)
}
