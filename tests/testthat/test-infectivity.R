library(testthat)

aux_add <- function(self, func, status = NA_integer_,
                    discr = FALSE, ..., max = NULL) {
  self$add_quo(func = func, status = status, discr = discr,
               rlang::enquos(...), max = rlang::enquo(max))
  invisible(self)
}

test_that("Default infectivity profile definitions are correct", {
  infec <- HousepiInfecProfile$new()

  # Add a profile with default arguments
  aux_add(infec, func = "dummy", discr = FALSE, a = 1, b = 2, max = 10)

  expect_equal(infec$defs$func, "dummy|a:b|{1|N}:{2|N}")
  expect_equal(infec$defs$status, NA_integer_)
  expect_equal(infec$defs$discr, FALSE)
  expect_equal(infec$defs$max, "10|N")
})

test_that("Infectivity profile with discrete values is correct", {
  infec <- HousepiInfecProfile$new()

  # Add a discrete profile
  aux_add(infec, func = "discrete_profile", discr = TRUE,
          x = 5, max = 20)

  expect_equal(infec$defs$func, "discrete_profile|x|{5|N}")
  expect_equal(infec$defs$discr, TRUE)
  expect_equal(infec$defs$max, "20|N")
})

test_that("Error handling in infectivity profile definitions", {
  infec <- HousepiInfecProfile$new()

  # Missing function name
  expect_error(aux_add(infec, func = NA, discr = FALSE, a = 1),
               regexp = "Incomplete definition")
})

test_that("Status specific definition work fine", {
  infec <- HousepiInfecProfile$new()

  aux_add(infec, func = "discrete_profile", discr = TRUE,
          x = 5, max = 20)

  expect_warning(aux_add(infec, func = "another_generic"),
                 regexp = "Overwriting")

  helper_test_aux_load(infec,
                       list(func = "lognormal", status = 1),
                       list(func = "ferretti20", status = 2),
                       v1 = "uniform")
  expect_equal(infec$defs$func, c("lognormal", "ferretti20", "uniform"))
  expect_equal(infec$defs$status, c(1, 2, NA_integer_))
})

test_that("Spurious names are dropped", {
  external_list <- list(
    ip1 = list(func = "gamma", status = NA,
               discr = TRUE, shape = ~ shape, scale = ~ scale),
    ip2 = list(func = "gamma", status = 3,
               discr = TRUE, shape = ~ shape, scale = ~ scale)
  )

  infec <- HousepiInfecProfile$new()
  expect_no_error(helper_test_aux_load(infec, !!!external_list))
  defs <- infec$defs
  defs$status <- NULL
  expect_equal(nrow(unique(defs)), 1)
})
