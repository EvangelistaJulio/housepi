library(testthat)

# model carrying no infectivity profile of its own
mk_infec <- local({
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template")
  list(model = hhmodel)
})

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
