library(testthat)

# Model carrying no infectivity profile of its own.
mk_infec_model <- function() new_housepi_model(model = "test")

## INFECTIVITY PROFILE : DEFINITIONS ##

test_that("Default infectivity profile definitions are correct", {
  # add a profile with default arguments
  infec <- mk_infec_model() |>
    add_infectivity_profile(func = "gamma", discr = FALSE,
                            shape = 1, scale = 2, max = 10) |>
    get_infectivity_profile()

  expect_equal(infec$func, "gamma")
  expect_equal(infec$func_args[[1]], c(shape = 1, scale = 2))
  expect_equal(infec$status, NA_real_)
  expect_equal(infec$discr, FALSE)
  expect_equal(infec$max, 10)
})

test_that("Infectivity profile with discrete values is correct", {
  # add a discrete profile
  infec <- mk_infec_model() |>
    add_infectivity_profile("ferretti20", discr = TRUE,
                            incub_period = 5, max = 20) |>
    get_infectivity_profile()

  expect_equal(infec$func, "ferretti20")
  expect_equal(infec$func_args[[1]], c(incub_period = 5))
  expect_equal(infec$discr, TRUE)
  expect_equal(infec$max, 20)
})

test_that("Error handling in infectivity profile definitions", {
  # missing function name
  mk_infec_model() |>
    add_infectivity_profile(func = NA, discr = FALSE, a = 1) |>
    expect_error("[Ii]ncomplete definition")
  # unknown family
  mk_infec_model() |>
    add_infectivity_profile(func = "fulton", x = 5) |>
    expect_error("[Uu]nknown infectivity profile")
})

test_that("Status specific definition work fine", {
  infec <- mk_infec_model() |>
    add_infectivity_profile("ferretti20", discr = TRUE,
                            incub_period = 5, max = 20) |>
    add_infectivity_profile(
      import = list(
        list(func = "lognormal", status = 1),
        list(func = "ferretti20", incub_period = 5, status = 2),
        v1 = "uniform"
      ),
      append = FALSE
    ) |>
    get_infectivity_profile()

  expect_equal(infec$func, c("lognormal", "ferretti20", "uniform"))
  expect_equal(infec$func_args[[2]], c(incub_period = 5))
  expect_equal(infec$status, c(1, 2, NA_real_))
})

test_that("Spurious names are dropped", {
  external_list <- list(
    ip1 = list(func = "gamma", status = NA,
               discr = TRUE, shape = ~ shape, scale = ~ scale),
    ip2 = list(func = "gamma", status = 2,
               discr = TRUE, shape = ~ shape, scale = ~ scale)
  )

  infec <- mk_infec_model() |>
    add_augmented("shape", model = "uniform", min = 1, max = 5) |>
    add_augmented("scale", model = "uniform", min = 1, max = 5) |>
    add_infectivity_profile(import = external_list) |>
    get_infectivity_profile()

  infec$status <- NULL
  expect_equal(nrow(unique(infec)), 1)
})
