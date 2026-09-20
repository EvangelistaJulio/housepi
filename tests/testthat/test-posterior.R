library(testthat)

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
  results <- data.frame(hhid = 1, id = c("m1", "m2"),
                        status = c(1, 0), obs = 0) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform"),
                  rho_2 = list(prior = "uniform")) |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden@A) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  id_draws = "m1") |>
    suppressWarnings()
  list(results = results)
})

##############################
## POSTERIOR : DEFINITIONS  ##
##############################

test_that("Invalid posterior summary settings are rejected", {
  set_global_posterior_summary(cri_level = 0.5) |>
    expect_error("between 0.5 and 1")
  set_global_posterior_summary(cri_level = 1) |>
    expect_error("between 0.5 and 1")
  set_global_posterior_summary(cri_level = NA_real_) |>
    expect_error("non-NA")
  set_global_posterior_summary(cri_level = "0.9") |>
    expect_error("numeric")
  set_global_posterior_summary(cri_level = c(0.9, 0.95)) |>
    expect_error("single")

  set_global_posterior_summary(point_estimate = "mode") |>
    expect_error("one between")
  # an abbreviation is not a choice
  set_global_posterior_summary(point_estimate = "med") |>
    expect_error("one between")
  set_global_posterior_summary(point_estimate = NA_character_) |>
    expect_error("non-NA")
  set_global_posterior_summary(point_estimate = 0.5) |>
    expect_error("character")
  set_global_posterior_summary(point_estimate = c("mean", "median")) |>
    expect_error("single")
})

test_that("A rejected setting leaves both values in place", {
  withr::with_options(list(housepi.point_estimate = "median",
                           housepi.cri_level = 0.95), {
    # the point estimate is valid, the credible level is not
    set_global_posterior_summary(point_estimate = "mean", cri_level = 2) |>
      expect_error("between 0.5 and 1")
    expect_equal(get_global_point_estimate(), "median")
    expect_equal(get_global_cri_level(), 0.95)
  })
})

test_that("The global posterior summary is set and read back", {
  withr::with_options(list(housepi.point_estimate = NULL,
                           housepi.cri_level = NULL), {
    expect_equal(get_global_point_estimate(), "median")
    expect_equal(get_global_cri_level(), 0.95)
  })
  withr::with_options(list(housepi.point_estimate = "median",
                           housepi.cri_level = 0.95), {
    # an argument left out keeps its value
    old <- set_global_posterior_summary(cri_level = 0.90)
    expect_equal(old, list(point_estimate = "median", cri_level = 0.95))
    expect_equal(get_global_cri_level(), 0.90)
    expect_equal(get_global_point_estimate(), "median")

    # the returned values restore what was there
    do.call(set_global_posterior_summary, old)
    expect_equal(get_global_cri_level(), 0.95)
  })
})

test_that("Summarise-draws default columns track cri_level", {
  plain <- posterior::example_draws()
  draws <- plain
  class(draws) <- unique(c("housepi_draws", class(draws)))

  # spelling alias dispatches to the same handler
  expect_identical(posterior::summarize_draws(draws),
                   posterior::summarise_draws(draws))

  # housepi method matches posterior::summarise_draws on every column
  # except the credible-interval quantile pair, which tracks cri_level.
  ref_cols <- setdiff(colnames(posterior::summarise_draws(plain)),
                      c("q5", "q95"))
  withr::with_options(list(housepi.cri_level = 0.95), {
    summary <- posterior::summarise_draws(draws)
    expect_true(all(c(ref_cols, "q2.5", "q97.5") %in% colnames(summary)))
  })
  withr::with_options(list(housepi.cri_level = 0.80), {
    summary <- posterior::summarise_draws(draws)
    expect_true(all(c(ref_cols, "q10", "q90") %in% colnames(summary)))
  })

  # custom summary functions fall through to vanilla posterior behaviour
  posterior::summarise_draws(draws, mean, sd) |>
    colnames() |>
    expect_setequal(c("variable", "mean", "sd"))

  # plain draws remain untouched by housepi (vanilla q5 / q95 quantiles)
  plain_cols <- colnames(posterior::summarise_draws(plain))
  expect_true(all(c("q5", "q95") %in% plain_cols))
  expect_false(any(c("q2.5", "q97.5", "q10", "q90") %in% plain_cols))
})

#############################
## POSTERIOR : END-TO-END  ##
#############################

test_that("Active fields are read-only", {
  mk_draws$results |>
    helper_test_read_only(
      "parameter_draws", "latent_draws",
      "log_lik_draws", "lprior_draws", "complete_draws",
      "warnings", "inference_seconds"
    )
})

test_that("Draws active fields are correct", {
  results <- mk_draws$results
  expect_true(posterior::is_draws(results$parameter_draws))
  expect_equal(posterior::niterations(results$parameter_draws), 100L)

  # no latent is traced
  expect_true(posterior::is_draws(results$latent_draws))
  expect_equal(posterior::nvariables(results$latent_draws), 0)

  expect_true(posterior::is_draws(results$log_lik_draws))
  expect_gt(posterior::nvariables(results$log_lik_draws), 0)

  expect_true(posterior::is_draws(results$lprior_draws))
  expect_true("lprior" %in% posterior::variables(results$lprior_draws))

  # parameter names, log_lik[...] and lprior (loo / priorsense names)
  expect_true(posterior::is_draws(results$complete_draws))
  vars <- posterior::variables(results$complete_draws)
  expect_true(any(grepl("log_lik\\[", vars)))
  expect_true("lprior" %in% vars)
})

test_that("Exported getter functions work correctly", {
  err <- "HousepiPosteriorDraws"
  expect_error(get_parameter_draws("not draws"), err)
  expect_error(get_latent_draws(42), err)
  expect_error(get_log_lik_draws(list()), err)
  expect_error(get_complete_draws(new_housepi_model(model = "test")), err)

  # the run has parameters but no latent draws
  get_latent_draws(mk_draws$results) |>
    expect_warning("id_draws")
})

test_that("Print behavior is stable", {
  local_results <- mk_draws$results
  expect_snapshot(print(local_results))
})

test_that("A short run records the diagnostic it skipped", {
  flagged <- mk_draws$results$warnings
  expect_equal(nrow(flagged), 1L)
  expect_match(flagged$message, "[Ss]ample size was below")
  expect_equal(flagged$count, 1L)
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
  expect_equal(quiet$warnings, loud$warnings)
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
  expect_warning(
    expect_warning(run(TRUE), "[Rr]educing .threads."),
    "[Ss]ample size was below"
  )
  expect_no_warning(results <- run(FALSE))
  # a build warning is held back but never recorded: only diagnostics are
  expect_false(any(grepl("threads", results$warnings$message)))
})

test_that("The sampling loop reports its own cost", {
  seconds <- mk_draws$results$inference_seconds
  expect_length(seconds, 1L)
  expect_true(is.finite(seconds) && seconds > 0)
})

test_that("Diagnostic warning fires on poor mixing", {
  flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(
      # explicitly define micro-steps
      alpha = list(prior = "uniform", proposal = "rw_log",
                   proposal_scale = 1e-4),
      beta  = list(prior = "uniform", proposal = "rw_log",
                   proposal_scale = 1e-4),
      gamma = 1, hhsize_ref = 5
    ) |>
    run_inference(steps = 1000L, thin = 1L, adapt = 0L, chains = 2L) |>
    expect_warning("alpha.*beta.*not mixed")
})

test_that("Pointwise log-lik is carried into the complete-draws bundle", {
  # complete_draws renames the pointwise columns to the loo / priorsense
  # convention; their values must be the same pointwise log-likelihood
  expect_same_pointwise <- function(results) {
    pointwise <- posterior::as_draws_matrix(get_log_lik_draws(results))
    bundle <- posterior::as_draws_matrix(get_complete_draws(results))
    ll_cols <- grep("^log_lik\\[", colnames(bundle))
    expect_equal(length(ll_cols), ncol(pointwise))
    expect_equal(unname(bundle[, ll_cols]), unname(pointwise))
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
  results <- data.frame(hhid = 1, id = "m1", status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden@A) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings()

  expect_output(print(results), "1 latent")
  # nothing carries a prior: the log-prior is zero at every draw
  posterior::extract_variable(results$lprior_draws, "lprior") |>
    expect_equal(rep(0, 5))
  get_parameter_draws(results) |>
    expect_warning("[Ee]mpty draws: no parameter")
})

test_that("Draws accessors subset the variables they return", {
  results <- mk_latent_draws$results

  # asking for nothing returns everything
  expect_setequal(posterior::variables(get_parameter_draws(results)),
                  c("rho", "rho_2"))
  expect_equal(posterior::variables(get_parameter_draws(results, "rho")),
               "rho")
  expect_equal(posterior::variables(get_parameter_draws(results, "rho_2")),
               "rho_2")

  # a latent trace is named `<id>_<latent>`: the latent name selects it across
  # the stored members, the full name selects one member
  expect_equal(posterior::variables(get_latent_draws(results)), "m1_hidden")
  expect_equal(posterior::variables(get_latent_draws(results, "hidden")),
               "m1_hidden")
  expect_equal(posterior::variables(get_latent_draws(results, "m1_hidden")),
               "m1_hidden")
  get_latent_draws(results, "idden") |>
    expect_error("[Nn]o draws of variable.*idden")

  # a base name selects its indexed elements
  expect_equal(posterior::variables(get_complete_draws(results, "log_lik")),
               "log_lik[1]")
  get_parameter_draws(results, c("rho", "hidden")) |>
    expect_error("[Nn]o draws of variable.*hidden")

  # the subset stays a housepi draws object, so summarise_draws still dispatches
  expect_s3_class(get_parameter_draws(results, "rho"), "housepi_draws")
  expect_s3_class(get_complete_draws(results, "lprior"), "housepi_draws")
})
