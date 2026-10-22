library(testthat)

####  Helpers  ####

# one short inference run, quiet, holding the warnings it raised
mk_draws <- local({
  results <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L,
                  warn = FALSE)
  list(results = results)
})

# two parameters and a latent traced on one member
mk_latent_draws <- local({
  results <- data.frame(hh_id = 1, id = c("m1", "m2"),
                        status = c(1, 0), obs = 0) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform"),
                  rho_2 = list(prior = "uniform")) |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  id_draws = "m1") |>
    suppressWarnings()
  list(results = results)
})

####  Definitions  ####

test_that("Invalid posterior summary settings are rejected", {
  set_housepi_summary(cri_level = 0.5) |>
    expect_error("between 0.5 and 1")
  set_housepi_summary(cri_level = 1) |>
    expect_error("between 0.5 and 1")
  set_housepi_summary(cri_level = NA_real_) |>
    expect_error("non-NA")
  set_housepi_summary(cri_level = "0.9") |>
    expect_error("numeric")
  set_housepi_summary(cri_level = c(0.9, 0.95)) |>
    expect_error("single")

  set_housepi_summary(point_estimate = "mode") |>
    expect_error("one of")
  # an abbreviation is not a choice
  set_housepi_summary(point_estimate = "med") |>
    expect_error("one of")
  set_housepi_summary(point_estimate = NA_character_) |>
    expect_error("non-NA")
  set_housepi_summary(point_estimate = 0.5) |>
    expect_error("character")
  set_housepi_summary(point_estimate = c("mean", "median")) |>
    expect_error("single")
})

test_that("A rejected setting leaves both values in place", {
  withr::with_options(
    list(housepi.point_estimate = "median", housepi.cri_level = 0.95), {
      # the point estimate is valid, the credible level is not
      set_housepi_summary(point_estimate = "mean", cri_level = 2) |>
        expect_error("between 0.5 and 1")
      get_global_point_estimate() |> expect_identical("median")
      get_global_cri_level() |> expect_equal(0.95)
    }
  )
})

test_that("The global posterior summary is set and read back", {
  withr::with_options(
    list(housepi.point_estimate = NULL, housepi.cri_level = NULL), {
      get_global_point_estimate() |> expect_identical("median")
      get_global_cri_level() |> expect_equal(0.9)
    }
  )
  withr::with_options(
    list(housepi.point_estimate = "median", housepi.cri_level = 0.95), {
      # an argument left out keeps its value
      old <- set_housepi_summary(cri_level = 0.90)
      old |> expect_equal(list(point_estimate = "median", cri_level = 0.95))
      get_global_cri_level() |> expect_equal(0.90)
      get_global_point_estimate() |> expect_identical("median")

      # the returned values restore what was there
      do.call(set_housepi_summary, old)
      get_global_cri_level() |> expect_equal(0.95)
    }
  )
})

test_that("Summarise-draws default columns track cri_level", {
  plain <- posterior::example_draws()
  draws <- plain
  class(draws) <- unique(c("housepi_draws", class(draws)))

  # spelling alias dispatches to the same handler
  posterior::summarize_draws(draws) |>
    expect_identical(posterior::summarise_draws(draws))

  # by default the method reproduces posterior::summarise_draws
  withr::with_options(list(housepi.cri_level = NULL), {
    posterior::summarise_draws(draws) |>
      colnames() |>
      expect_identical(colnames(posterior::summarise_draws(plain)))
  })

  # housepi method matches posterior::summarise_draws on every column
  # except the credible-interval quantile pair, which tracks cri_level.
  ref_cols <- setdiff(colnames(posterior::summarise_draws(plain)),
                      c("q5", "q95"))
  withr::with_options(list(housepi.cri_level = 0.95), {
    summary <- posterior::summarise_draws(draws)
    colnames(summary) |> expect_contains(c(ref_cols, "q2.5", "q97.5"))
  })
  withr::with_options(list(housepi.cri_level = 0.80), {
    summary <- posterior::summarise_draws(draws)
    colnames(summary) |> expect_contains(c(ref_cols, "q10", "q90"))
  })

  # custom summary functions fall through to vanilla posterior behaviour
  posterior::summarise_draws(draws, mean, sd) |>
    colnames() |>
    expect_setequal(c("variable", "mean", "sd"))

  # plain draws remain untouched by housepi (vanilla q5 / q95 quantiles)
  plain_cols <- colnames(posterior::summarise_draws(plain))
  plain_cols |> expect_contains(c("q5", "q95"))
  (c("q2.5", "q97.5", "q10", "q90") %in% plain_cols) |> expect_all_false()
})

####  End-to-end  ####

test_that("Active fields are read-only", {
  mk_draws$results |>
    expect_read_only_field(
      "parameter_draws", "latent_draws",
      "log_lik_draws", "lprior_draws", "complete_draws",
      "warnings", "inference_secs"
    )
})

test_that("Draws active fields are correct", {
  results <- mk_draws$results
  posterior::is_draws(results$parameter_draws) |> expect_true()
  posterior::niterations(results$parameter_draws) |> expect_identical(100L)

  # no latent is traced
  posterior::is_draws(results$latent_draws) |> expect_true()
  posterior::nvariables(results$latent_draws) |> expect_identical(0L)

  posterior::is_draws(results$log_lik_draws) |> expect_true()
  posterior::nvariables(results$log_lik_draws) |> expect_gt(0)

  posterior::is_draws(results$lprior_draws) |> expect_true()
  ("lprior" %in% posterior::variables(results$lprior_draws)) |> expect_true()

  # parameter names, log_lik[...] and lprior (loo / priorsense names)
  posterior::is_draws(results$complete_draws) |> expect_true()
  vars <- posterior::variables(results$complete_draws)
  vars |> expect_match("log_lik\\[", all = FALSE)
  ("lprior" %in% vars) |> expect_true()
})

test_that("Exported getter functions work correctly", {
  err <- "HousepiPosteriorDraws"
  get_parameter_draws("not draws") |> expect_error(err)
  get_latent_draws(42) |> expect_error(err)
  get_log_lik_draws(list()) |> expect_error(err)
  get_complete_draws(new_housepi_model(model = "test")) |> expect_error(err)

  # the run has parameters but no latent draws
  get_latent_draws(mk_draws$results) |>
    expect_warning("id_draws")
})

test_that("Print behavior is stable", {
  local_results <- mk_draws$results
  print(local_results) |> expect_snapshot(transform = helper_scrub_secs)
})

test_that("A short run records one warning, under its sample size", {
  flagged <- mk_draws$results$warnings
  flagged |> expect_nrow(1L)
  flagged$message |> expect_match("[Ss]ample size was below")
  flagged$count |> expect_identical(1L)
})

test_that("A quiet run holds back the warnings it records", {
  run <- function(warn) {
    flu_1918_frederick |>
      new_housepi_model("final-size") |>
      add_parameter(import = "template") |>
      run_inference(steps = 20L, chains = 1L, adapt = 0L, thin = 1L,
                    warn = warn)
  }
  # the same run is silent or not, and records the same warnings either way
  expect_no_warning(quiet <- run(FALSE))
  expect_warning(loud <- run(TRUE), "[Ss]ample size was below")
  quiet$warnings |> expect_identical(loud$warnings)
})

test_that("A quiet run holds back its build warnings too", {
  # more threads than the machine has, which the engine reduces and reports
  threads <- RcppParallel::defaultNumThreads() + 4L
  run <- function(warn) {
    flu_1918_frederick |>
      new_housepi_model("final-size") |>
      add_parameter(import = "template") |>
      run_inference(steps = 20L, chains = 1L, adapt = 0L, thin = 1L,
                    threads = threads, warn = warn)
  }
  # the run is too short to diagnose, so it warns twice
  run(TRUE) |>
    expect_warning("[Rr]educing .threads.") |>
    expect_warning("[Ss]ample size was below")
  expect_no_warning(results <- run(FALSE))
  # a build warning is held back but never recorded: only diagnostics are
  results$warnings$message |> expect_no_match("threads")
})

test_that("The sampling loop reports its own cost", {
  secs <- mk_draws$results$inference_secs
  secs |> expect_length(1L)
  secs |> expect_all_finite()
  secs |> expect_gt(0)
})

test_that("Diagnostic warning fires on poor mixing", {
  expect_warning(
    results <- flu_1918_frederick |>
      new_housepi_model("final-size") |>
      add_parameter(
        # explicitly define micro-steps
        alpha = list(prior = "uniform", proposal = "rw_log",
                     proposal_scale = 1e-4),
        beta  = list(prior = "uniform", proposal = "rw_log",
                     proposal_scale = 1e-4),
        gamma = 1, hh_size_ref = 5
      ) |>
      run_inference(steps = 1000L, thin = 1L, adapt = 0L, chains = 2L),
    "alpha.*beta.*not mixed"
  )
  # the recorded line states what was flagged without the findings under it
  results$warnings$message |> expect_match("[Ff]lagged parameter draws")
})

test_that("Pointwise log-lik is carried into the complete-draws bundle", {
  # complete_draws renames the pointwise columns to the loo / priorsense
  # convention; their values must be the same pointwise log-likelihood
  expect_same_pointwise <- function(results) {
    pointwise <- draws_to_frame(get_log_lik_draws(results))
    bundle <- draws_to_frame(get_complete_draws(results))
    ll_cols <- grep("^log_lik\\[", colnames(bundle))
    length(ll_cols) |> expect_identical(ncol(pointwise))
    unname(bundle[, ll_cols]) |> expect_equal(unname(pointwise))
  }
  expect_same_pointwise(mk_draws$results)

  # the pointwise store is filled over two threads alike
  flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference(steps = 10L, thin = 2L, adapt = 2L, burn = 2L,
                  chains = 1L, threads = 2L) |>
    suppressWarnings() |>
    expect_same_pointwise()
})

test_that("A model without parameters returns empty parameter draws", {
  results <- data.frame(hh_id = 1, id = "m1", status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings()

  print(results) |> expect_message("Latent: 1")
  # nothing carries a prior: no lprior column, in its field or the bundle
  posterior::nvariables(results$lprior_draws) |> expect_identical(0L)
  posterior::variables(get_complete_draws(results)) |>
    expect_identical("log_lik[1]")
  get_parameter_draws(results) |>
    expect_warning("[Ee]mpty draws: no parameter")
})

test_that("Draws accessors subset the variables they return", {
  results <- mk_latent_draws$results

  # asking for nothing returns everything
  posterior::variables(get_parameter_draws(results)) |>
    expect_setequal(c("rho", "rho_2"))
  posterior::variables(get_parameter_draws(results, "rho")) |>
    expect_identical("rho")
  posterior::variables(get_parameter_draws(results, "rho_2")) |>
    expect_identical("rho_2")

  # a latent trace is named `<id>_<latent>`: the latent name selects it across
  # the stored members, the full name selects one member
  posterior::variables(get_latent_draws(results)) |>
    expect_identical("m1_hidden")
  posterior::variables(get_latent_draws(results, "hidden")) |>
    expect_identical("m1_hidden")
  posterior::variables(get_latent_draws(results, "m1_hidden")) |>
    expect_identical("m1_hidden")
  get_latent_draws(results, "idden") |>
    expect_error("[Nn]o draws of variable.*idden")

  # a base name selects its indexed elements
  posterior::variables(get_complete_draws(results, "log_lik")) |>
    expect_identical("log_lik[1]")
  get_parameter_draws(results, c("rho", "hidden")) |>
    expect_error("[Nn]o draws of variable.*hidden")

  # the subset stays a housepi draws object, so summarise_draws still dispatches
  get_parameter_draws(results, "rho") |> expect_s3_class("housepi_draws")
  get_complete_draws(results, "lprior") |> expect_s3_class("housepi_draws")
})
