library(testthat)

####  Helpers  ####

# one short inference run, quiet, holding the warnings it raised
mk_draws <- local({
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  hh_fit <- hh_model |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L,
                  warn = FALSE)
  list(model = hh_model, fit = hh_fit)
})

# two parameters and a latent traced on one member
mk_latent_draws <- local({
  hh_model <- data.frame(hh_id = 1, id = c("m1", "m2"),
                         status = c(1, 0), obs = 0) |>
    new_housepi_model("test") |>
    # declared out of order: draws follow get_parameters(), not declaration
    add_parameter(rho_2 = list(prior = "uniform"),
                  rho = list(prior = "uniform")) |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden)
  hh_fit <- hh_model |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  id_draws = "m1") |>
    suppressWarnings()
  list(model = hh_model, fit = hh_fit)
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
  mk_draws$fit |>
    expect_read_only_field("warnings", "inference_secs")
})

test_that("Draws getters return each block", {
  hh_fit <- mk_draws$fit
  param <- get_parameter_draws(hh_fit)
  param |> expect_s3_class("draws")
  posterior::niterations(param) |> expect_identical(100L)

  # no latent is traced
  get_latent_draws(hh_fit) |> expect_warning("[Ee]mpty draws")

  log_lik <- get_log_lik_draws(hh_fit)
  log_lik |> expect_s3_class("draws")
  posterior::nvariables(log_lik) |> expect_gt(0)

  # parameter names, log_lik[...] and lprior (loo / priorsense names)
  complete <- get_complete_draws(hh_fit)
  complete |> expect_s3_class("draws")
  n_house <- posterior::nvariables(log_lik)
  posterior::variables(complete) |>
    expect_identical(c(posterior::variables(param),
                       paste0("log_lik[", seq_len(n_house), "]"), "lprior"))
  posterior::variables(get_complete_draws(hh_fit, "lprior")) |>
    expect_identical("lprior")
})

test_that("Exported getter functions work correctly", {
  err <- "HousepiPosteriorDraws"
  get_parameter_draws("not draws") |> expect_error(err)
  get_latent_draws(42) |> expect_error(err)
  get_log_lik_draws(list()) |> expect_error(err)
  get_complete_draws(new_housepi_model(model = "test")) |> expect_error(err)

  # the run has parameters but no latent draws
  get_latent_draws(mk_draws$fit) |>
    expect_warning("id_draws")
})

test_that("Print behavior is stable", {
  local_results <- mk_draws$fit
  print(local_results) |> expect_snapshot(transform = helper_scrub_secs)
})

test_that("A short run records one warning per flagged draws block", {
  # 100 draws on one chain: the ESS stays below 100
  flagged <- mk_draws$fit$warnings
  flagged |> expect_nrow(1L)
  flagged$message |> expect_match("[Ff]lagged parameter draws")
  flagged$count |> expect_identical(1L)
  # the print closes on the warning count
  helper_print_lines(print(mk_draws$fit)) |>
    expect_match("^! Warning: 1 issue flagged, see .\\$warnings", all = FALSE)
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
  expect_warning(loud <- run(TRUE), "ESS.*below 100 per chain")
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
  # 20 draws cap the ESS below 100, so the run warns twice
  run(TRUE) |>
    expect_warning("[Rr]educing .threads.") |>
    expect_warning("ESS.*below 100 per chain")
  expect_no_warning(hh_fit <- run(FALSE))
  # a build warning is held back but never recorded: only diagnostics are
  hh_fit$warnings$message |> expect_no_match("threads")
})

test_that("The sampling loop reports its own cost", {
  secs <- mk_draws$fit$inference_secs
  secs |> expect_length(1L)
  secs |> expect_all_finite()
  secs |> expect_gt(0)
})

test_that("Diagnostic warning fires on poor mixing", {
  expect_warning(
    hh_fit <- flu_1918_frederick |>
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
  hh_fit$warnings$message |> expect_match("[Ff]lagged parameter draws")
})

test_that("Pointwise log-lik is carried into the complete-draws bundle", {
  # get_complete_draws() renames the pointwise columns to the loo / priorsense
  # convention; their values must be the same pointwise log-likelihood
  expect_same_pointwise <- function(hh_fit) {
    pointwise <- helper_draws_to_frame(get_log_lik_draws(hh_fit))
    bundle <- helper_draws_to_frame(get_complete_draws(hh_fit))
    ll_cols <- grep("^log_lik\\[", colnames(bundle))
    ll_cols |> expect_length(ncol(pointwise))
    unname(bundle[, ll_cols]) |> expect_equal(unname(pointwise))
  }
  expect_same_pointwise(mk_draws$fit)

  # the pointwise store is filled over two threads alike
  flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference(steps = 10L, thin = 2L, adapt = 2L,
                  chains = 1L, threads = 2L) |>
    suppressWarnings() |>
    expect_same_pointwise()
})

test_that("A model without parameters returns empty parameter draws", {
  hh_fit <- data.frame(hh_id = 1, id = "m1", status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings()

  print(hh_fit) |> expect_message("Latent: 1")
  # nothing carries a prior: no lprior column in the bundle
  posterior::variables(get_complete_draws(hh_fit)) |>
    expect_identical("log_lik[1]")
  get_parameter_draws(hh_fit) |>
    expect_warning("[Ee]mpty draws: no parameter")
})

test_that("Draws accessors subset the variables they return", {
  hh_fit <- mk_latent_draws$fit

  # asking for nothing returns everything, in the order of get_parameters()
  posterior::variables(get_parameter_draws(hh_fit)) |>
    expect_identical(get_parameters(mk_latent_draws$model)$name)
  posterior::variables(get_parameter_draws(hh_fit, "rho")) |>
    expect_identical("rho")
  posterior::variables(get_parameter_draws(hh_fit, "rho_2")) |>
    expect_identical("rho_2")

  # a latent trace is named `<id>_<latent>`: the latent name selects it across
  # the stored members, the full name selects one member
  posterior::variables(get_latent_draws(hh_fit)) |>
    expect_identical("m1_hidden")
  posterior::variables(get_latent_draws(hh_fit, "hidden")) |>
    expect_identical("m1_hidden")
  posterior::variables(get_latent_draws(hh_fit, "m1_hidden")) |>
    expect_identical("m1_hidden")
  get_latent_draws(hh_fit, "idden") |>
    expect_error("[Nn]o draws of variable.*idden")

  # a base name selects its indexed elements
  posterior::variables(get_complete_draws(hh_fit, "log_lik")) |>
    expect_identical("log_lik[1]")
  get_parameter_draws(hh_fit, c("rho", "hidden")) |>
    expect_error("[Nn]o draws of variable.*hidden")

  # the subset stays a housepi draws object, so summarise_draws still dispatches
  get_parameter_draws(hh_fit, "rho") |> expect_s3_class("housepi_draws")
  get_complete_draws(hh_fit, "lprior") |> expect_s3_class("housepi_draws")
})

test_that("Runs of one model bind their chains", {
  hh_fit <- mk_draws$fit
  # the model built anew, as in another session
  hh_fit_2 <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference(steps = 100L, chains = 2L, adapt = 0L, thin = 1L,
                  warn = FALSE)
  # diagnostics run again over the chains bound, and warn as a run does
  combined <- c(hh_fit, hh_fit_2) |> suppressWarnings()

  combined |> expect_s3_class("HousepiPosteriorDraws")
  param <- get_parameter_draws(combined)
  posterior::nchains(param) |> expect_identical(3L)
  posterior::niterations(param) |> expect_identical(100L)
  # the chains of the second run follow those of the first
  posterior::subset_draws(param, chain = 2:3)$alpha |>
    expect_identical(get_parameter_draws(hh_fit_2)$alpha)
  posterior::nchains(get_log_lik_draws(combined)) |> expect_identical(3L)
  combined$inference_secs |>
    expect_equal(hh_fit$inference_secs + hh_fit_2$inference_secs)
  combined$warnings$message |> expect_match("[Ff]lagged parameter draws")
  # the runs combined are left as they were
  posterior::nchains(get_parameter_draws(hh_fit)) |> expect_identical(1L)
})

test_that("Runs of different settings, data or members do not combine", {
  hh_fit <- mk_draws$fit
  run <- function(hh_model, steps = 100L, ...) {
    run_inference(hh_model, steps = steps, chains = 1L, adapt = 0L,
                  thin = 1L, warn = FALSE, ...)
  }
  c(hh_fit, hh_fit) |> expect_error("[Rr]epeat a seed")
  c(hh_fit, run(mk_draws$model, steps = 50L)) |>
    expect_error("[Dd]iffer in.*steps")
  hh_model <- flu_1918_frederick[-1L, ] |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  c(hh_fit, run(hh_model)) |> expect_error("[Dd]iffer in.*dataset")
  hh_model <- mk_draws$model |>
    add_parameter("beta", prior = "uniform", min = 0, max = 2, append = TRUE) |>
    suppressWarnings()
  c(hh_fit, run(hh_model)) |> expect_error("[Dd]iffer in.*parameters")
  # `[^.]` keeps to the fields, the hint names the model too
  c(hh_fit, mk_latent_draws$fit) |>
    expect_error("[Dd]iffer in[^.]*model")
  hh_fit_2 <- run(mk_latent_draws$model, steps = 5L, id_draws = "m2")
  c(mk_latent_draws$fit, hh_fit_2) |>
    expect_error("[Dd]iffer in.*id_draws")
  # a default run traces one household at random, among two here:
  # one of twenty runs traces the other one but once in 500,000 suites
  hh_model <- data.frame(hh_id = c("h1", "h1", "h2", "h2"),
                         id = c("m1", "m2", "m3", "m4"),
                         status = c(1, 0, 1, 0), obs = 0) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform")) |>
    add_latent(hidden = list(model = "uniform")) |>
    helper_flat_sink(loc = hidden)
  traced <- \(x) posterior::variables(get_latent_draws(x))
  runs <- replicate(20L, run(hh_model, steps = 5L), simplify = FALSE)
  other <- Find(\(x) !identical(traced(x), traced(runs[[1L]])), runs)
  c(runs[[1L]], other) |> expect_error("[Dd]ifferent variables")
})
