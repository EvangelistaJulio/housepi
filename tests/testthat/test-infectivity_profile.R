library(testthat)

# model carrying no infectivity profile of its own
mk_infec <- local({
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template")
  list(model = hhmodel)
})

# same model, but in discrete-time
mk_infec_discr <- local({
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", discr = TRUE, quiet = TRUE) |>
    add_parameter(import = "template") |>
    # UPDATE: 'infec_time' is declared discrete here because
    # template leaves it continuous, which a discrete-time model refuses
    add_augmented(
      infec_time = list(model = "uniform", discr = TRUE,
                        min = start_date@D - 10,  max = end_date@D)
    )
  list(model = hhmodel)
})

############################################################################
################ INFECTIVITY PROFILE : DEFINITIONS #########################
############################################################################

test_that("error handling in infectivity profile definitions", {
  # missing function name
  mk_infec$model |>
    add_infectivity_profile(func = NA, discr = FALSE, a = 1) |>
    expect_error("[Ii]ncomplete definition")
  # unknown family
  mk_infec$model |>
    add_infectivity_profile(func = "fulton", x = 5) |>
    expect_error("[Uu]nknown infectivity profile")
  # mixing definition modes
  mk_infec$model |>
    add_infectivity_profile("sokolov", import = "template") |>
    expect_error("[Mm]ixing")
  # importing from unknown source
  mk_infec$model |>
    add_infectivity_profile(import = "templte") |>
    expect_error("[Uu]nknown import source")
})

test_that("infectivity profile definitions are correct", {
  defs <- mk_infec$model |>
    add_infectivity_profile(func = "gamma", shape = 1, scale = 2,
                            min = 2, max = 10) |>
    get_infectivity_profile()

  expect_equal(defs$func, "gamma")
  expect_equal(defs$func_args[[1]], c(shape = 1, scale = 2))
  expect_equal(defs$status, NA_real_)
  expect_equal(defs$discr, FALSE)
  expect_equal(defs$min, 2)
  expect_equal(defs$max, 10)
})

test_that("infectivity profile with discrete passes on continuos time", {
  # add a discrete profile
  defs <- mk_infec$model |>
    add_infectivity_profile("ferretti20", discr = TRUE,
                            incub_period = 5, max = 20) |>
    get_infectivity_profile()

  expect_equal(defs$func, "ferretti20")
  expect_equal(defs$func_args[[1]], c(incub_period = 5))
  expect_equal(defs$discr, TRUE)
  expect_equal(defs$max, 20)
})

test_that("status specific definition work fine", {
  defs <- mk_infec$model |>
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

  expect_equal(defs$func, c("lognormal", "ferretti20", "uniform"))
  expect_equal(defs$func_args[[2]], c(incub_period = 5))
  expect_equal(defs$status, c(1, 2, NA_real_))
})

test_that("spurious names are dropped", {
  external_list <- list(
    ip1 = list(func = "gamma", status = NA,
               discr = TRUE, shape = ~ shape@P, scale = ~ scale@P),
    ip2 = list(func = "gamma", status = 2,
               discr = TRUE, shape = ~ shape@P, scale = ~ scale@P)
  )

  infec <- mk_infec$model |>
    add_parameter(shape = 4, scale = 2) |>
    add_infectivity_profile(import = external_list) |>
    get_infectivity_profile()

  infec$status <- NULL
  expect_equal(nrow(unique(infec)), 1)
})

test_that("range minimum defaults to the first transmission day", {
  profile_min <- function(hhmodel, ...) {
    hhmodel |>
      add_infectivity_profile(...) |>
      get_infectivity_profile() |>
      getElement("min")
  }

  # discrete profile change it even in continuous time
  profile_min(mk_infec$model, "gamma", mean = 4, var = 3, discr = TRUE) |>
    expect_equal(1)
  # a user-supplied minimum wins in either build
  profile_min(mk_infec_discr$model, "gamma", mean = 4, var = 3, min = 2) |>
    expect_equal(2)
  # template definitions return minimum according to model time choice
  profile_min(mk_infec$model, import = "template") |>
    expect_equal(0)
  profile_min(mk_infec_discr$model, import = "template") |>
    expect_equal(1)
})

############################################################################
################ INFECTIVITY PROFILE : END-TO-END ##########################
############################################################################

test_that("a model needing a profile throws", {
  mk_infec$model |>
    get_infectivity_profile() |>
    expect_null()
  mk_infec$model |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing infectivity profile")
  mk_infec$model |>
    run_inference() |>
    expect_error("[Mm]issing infectivity profile")
  mk_infec$model |>
    run_simulation() |>
    expect_error("[Mm]issing infectivity profile")
})

test_that("a profile transmitting at generation time 0 is refused", {
  # distribution with non-zero at 0
  mk_infec$model |>
    add_infectivity_profile("gamma", mean = 4, var = 3, loc = -1) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation")

  # FIXME: disabled assert in simulation because of clash with status-specific
  # same behavior must occur in simulation
  # mk_infec$model |>
  #   add_infectivity_profile("exponential") |>
  #   run_simulation() |>
  #   expect_error("non-zero at generation")

  # discretization resolves it automatically
  mk_infec$model |>
    add_infectivity_profile("exponential", discr = TRUE) |>
    evaluate_log_lik() |>
    expect_silent()

  # a registered profile considered as authored
  # NOTE: 'ferretti20' is a custom with non-zero pdf at origin
  mk_infec$model |>
    add_infectivity_profile("ferretti20", incub_period = 5) |>
    evaluate_log_lik() |>
    expect_silent()
})

test_that("non-zero profile for day-0 is always refused in discrete time", {
  mk_infec_discr$model |>
    add_infectivity_profile("gamma", mean = 4, var = 3, min = 0) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation day 0")
  mk_infec_discr$model |>
    add_infectivity_profile("ferretti20", incub_period = 5, min = 0) |>
    evaluate_log_lik() |>
    expect_error("non-zero at generation day 0")
})
