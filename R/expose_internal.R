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

#' @title Expose housepi function
#' @description
#' Expose an object from housepi registered function.
#' Used to access user defined functions.
#' @param name `r sprintf(rd_test_name_, "function")`
#' @param ... `r sprintf(rd_test_args_, "function")`
#' @param min `r rd_var_min_`
#' @param max `r rd_var_max_`
#' @param discr `r rd_test_discr_`
#' @return `ExposeFunction` object.
#' @examples
#' fnc <- expose_function("pow", x2 = 3)
#' fnc$func(seq(0, 5))
#' fnc$prim(seq(0, 5))
#' @export
expose_function <- function(name, ...,
                            min = NA_real_,
                            max = NA_real_,
                            discr = FALSE) {
  new_obj <- .engine$ExposeFunction$new(name)
  new_obj$init(args = unlist(list(...)), min = min, max = max, discr = discr)
  new_obj
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
#' @return `ExposeDistribution` object.
#' @examples
#' dst <- expose_distribution("weibull", shape = 2, scale = 1)
#' dst$pdf(seq(0, 5))
#' dst$cdf(seq(0, 5))
#' @export
expose_distribution <- function(name, ...,
                                min = NA_real_,
                                max = NA_real_,
                                discr = FALSE) {
  new_obj <- .engine$ExposeDistribution$new(name)
  new_obj$init(args = unlist(list(...)), min = min, max = max, discr = discr)
  new_obj
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
#' @return `ExposeInfecProfile` object.
#' @examples
#' ipf <- expose_infec_profile("gamma", shape = 2, scale = 1)
#' ipf$inst(seq(0, 5))
#' ipf$cumul(seq(0, 5))
#' @export
expose_infec_profile <- function(name, ...,
                                 min = NA_real_,
                                 max = NA_real_,
                                 discr = FALSE) {
  if (is.na(min)) min <- ifelse(discr, 1, 0)
  new_obj <- .engine$ExposeInfecProfile$new(name)
  new_obj$init(args = unlist(list(...)), min = min, max = max, discr = discr)
  new_obj
}
