library(testthat)

######################################
###  PARAM RECOVERY : DEFINITIONS  ###
######################################

# households of two sizes, the covariate constant within each one
mk_recov_data <- local({
  hhdata <- data.frame(hh_id = c("h1", rep("h2", 3)),
                       id = paste0("m", seq(4)),
                       datum = c("a", rep("b", 3)))
  list(hhdata = hhdata)
})

# one replicate of two parameters at two household counts
# NOTE: values chosen so every statistic is readable by hand
mk_recov_frame <- local({
  estimates <- data.frame(
    houses = rep(c(10L, 20L), each = 2L), rep = 1L,
    variable = rep(c("rho", "rho_2"), 2L),
    truth = c(1, 0, 1, 0),
    estimate = c(2, 1, 0, 3),
    cri_lower = c(0, 2, -1, 1), cri_upper = c(4, 3, 1, 5),
    covers = c(TRUE, FALSE, TRUE, FALSE),
    rhat = 1, ess_bulk = c(100, 200, 300, 400),
    flagged = c(FALSE, TRUE, FALSE, FALSE)
  )
  list(estimates = estimates)
})

test_that("The household count defaults to the observed one", {
  sample_house_structure(data.frame(datum = 1)) |>
    expect_error("[Cc]olumn.*hh_id.*is missing")
  sample_house_structure(mk_recov_data$hhdata, houses = 0) |>
    expect_error("[Nn]umber of households")
  # the sampler draws one structure, so a sweep of counts is refused here
  sample_house_structure(mk_recov_data$hhdata, houses = c(2, 3)) |>
    expect_error("[Nn]umber of households")

  # default sampling targets the observed household count
  struct <- sample_house_structure(mk_recov_data$hhdata)
  length(unique(struct$hh_id)) |> expect_identical(2L)
  # without a sampler the observed count stands, whatever count is asked
  expect_warning(
    houses <- sweep_houses_(mk_recov_data$hhdata, houses = 6, verbatim = TRUE),
    "[Aa]rgument.*houses.*is ignored"
  )
  houses |> expect_identical(2L)
})

test_that("Households are resampled whole", {
  struct <- sample_house_structure(mk_recov_data$hhdata, houses = 6)
  # every household is renumbered
  sort(unique(struct$hh_id)) |> expect_identical(1:6)
  # an identifier carried by the input is regenerated
  struct$id |> expect_identical(seq_len(nrow(struct)))
  # household sizes are only the observed ones
  sizes <- as.integer(table(struct$hh_id))
  sizes |> expect_in(c(1L, 3L))
  # covariate datum is still unique per household and not mixed
  levels <- tapply(struct$datum, struct$hh_id, \(x) length(unique(x)))
  levels |> as.vector() |> expect_all_equal(1L)
})

test_that("Recovery statistics summarise the replicates", {
  summary <- recovery_summary_(mk_recov_frame$estimates)
  summary$houses |> expect_identical(c(10L, 10L, 20L, 20L))
  summary$variable |> expect_identical(c("rho", "rho_2", "rho", "rho_2"))
  summary$bias |> expect_equal(c(1, 1, -1, 3))
  summary$rmse |> expect_equal(c(1, 1, 1, 3))
  # relative bias is undefined against a zero truth
  summary$rel_bias |> expect_equal(c(1, NA, -1, NA))
  # the interval misses the truth of `rho_2` at either count
  summary$coverage |> expect_equal(c(1, 0, 1, 0))
  summary$cri_width |> expect_equal(c(4, 1, 2, 4))
  summary$n_flagged |> expect_identical(c(0L, 1L, 0L, 0L))
})

test_that("Sweep cost is averaged over the replicates", {
  replicates <- data.frame(
    houses = c(10L, 10L, 20L), rep = c(1L, 2L, 1L),
    inference_secs = c(1, 3, 5), min_ess_bulk = c(100, 600, NA),
    flagged = c(FALSE, TRUE, FALSE)
  )
  overview <- sweep_overview_(replicates)
  overview$nrep |> expect_identical(c(2L, 1L))
  overview$n_flagged |> expect_identical(c(1L, 0L))
  overview$mean_inference_secs |> expect_equal(c(2, 5))
  # pooled over the replicates: (1 + 3) / (100 + 600)
  # a fit without effective sample size carries no cost per draw
  overview$secs_per_ess |> expect_equal(c(4 / 700, NaN))
})

test_that("A table of values collapses to its point estimate", {
  truth <- data.frame(rho = c(1, 2, 6), rho_2 = c(0, 0, 3))
  recovery_collapse_truth_(truth, "median") |>
    expect_equal(c(rho = 2, rho_2 = 0))
  recovery_collapse_truth_(truth, "mean") |>
    expect_equal(c(rho = 3, rho_2 = 1))
})

#####################################
###  PARAM RECOVERY : END-TO-END  ###
#####################################

# a sweep short enough for the suite: every fit is too short to diagnose,
# so the whole run is flagged and its one aggregated warning is asserted below
mk_recov <- local({
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  results <- hh_model |>
    check_recovery(nrep = 2, houses = c(20, 40),
                   steps = 50, chains = 1, progress = FALSE) |>
    suppressWarnings()
  list(model = hh_model, results = results)
})

# no fit raised a warning, no fit is flagged and every point estimate lies
# within four posterior sd of its simulating value, which a calibrated
# posterior strays past about once in 16,000 fits
# the check is run under `suppressWarnings()`: its warnings are read back here
expect_recovered <- function(results) {
  results$warnings$message |> expect_identical(character())
  estimates <- results$estimates
  estimates$flagged |> expect_all_false()
  (abs(estimates$estimate - estimates$truth) < 4 * estimates$sd) |>
    expect_all_true()
}

test_that("Recovery results are read only", {
  mk_recov$results |>
    expect_read_only_field("estimates", "replicates", "summary",
                           "overview", "truth", "houses", "nrep",
                           "point_estimate", "cri_level", "warnings")
})

test_that("Recovery results print their settings", {
  results <- mk_recov$results
  lines <- helper_print_lines(printed <- print(results))
  # printing returns the object itself
  printed |> expect_identical(results)
  lines[1] |> expect_match("[Pp]arameter recovery check.*final size")
  lines |> expect_contains(c("Households: 20 and 40", "Replicates each: 2",
                             "Point estimate: median",
                             "Credible interval: 90%"))
  # the count is that of the flagged replicates
  lines |>
    expect_contains(paste("Flagged by diagnostics:",
                          sum(results$replicates$flagged)))
  # durations read as the progress bar prints them
  lines |> expect_match("^Run time: [0-9]+s$", all = FALSE)
  lines |> expect_match("^Inference time per fit: [0-9]+s$", all = FALSE)
})

test_that("Recovery sweeps every household count", {
  results <- mk_recov$results
  results$houses |> expect_identical(c(20L, 40L))
  results$truth |> expect_equal(c(alpha = 0.15, beta = 0.20, gamma = 1))
  # two counts, two replicates, three estimated parameters
  results$estimates |> expect_nrow(12L)
  sort(unique(results$estimates$variable)) |>
    expect_identical(c("alpha", "beta", "gamma"))
  results$replicates |> expect_nrow(4L)
  results$summary |> expect_nrow(6L)
  results$overview$houses |> expect_identical(c(20L, 40L))
  results$nrep |> expect_identical(2L)
  results$point_estimate |> expect_identical("median")
  results$cri_level |> expect_equal(0.9)
})

test_that("Each replicate reports its cost and its survey", {
  replicates <- mk_recov$results$replicates
  (replicates$inference_secs > 0) |> expect_all_true()
  estimates <- mk_recov$results$estimates
  replicates$min_ess_bulk |>
    expect_equal(as.vector(tapply(estimates$ess_bulk,
                                  estimates[c("rep", "houses")], min)))
  (mk_recov$results$overview$secs_per_ess > 0) |> expect_all_true()
  mk_recov$results$summary$coverage |> expect_in(c(0, 0.5, 1))
  replicates$seed |> expect_type("integer")
  replicates$seed |> expect_no_na()
  estimates$covers |>
    expect_identical(estimates$cri_lower <= estimates$truth &
                       estimates$truth <= estimates$cri_upper)
})

test_that("A failed replicate is named in the error", {
  hhdata <- sample_house_structure(mk_recov$model$data_wide, houses = 20)
  # every status missing passes the model build and stops the engine before
  # any fit samples
  hhdata_2 <- hhdata
  hhdata_2$status <- NA_integer_
  fit_args <- list(steps = 50, burn = 0, thin = 1, adapt = 0, chains = 1,
                   threads = 1L, id_draws = ".none")
  sweep_fit_wave_(mk_recov$model, list(hhdata, hhdata_2), reps = 1:2,
                  houses = 20L, seeds = 1:2, fit_args = fit_args,
                  outer = 2L, bar = NULL) |>
    expect_error("[Rr]eplicate 2 at.*20.*households.*missing status")
})

test_that("Recovery fits the survey the wide reporting process returns", {
  # the wide process empties every status and fails the fit, the long one
  # would fail the simulation if it ran
  mk_recov$model |>
    add_reporting_process(wide = \(x) transform(x, status = NA_integer_),
                          long = \(x) stop("long reporting process ran")) |>
    check_recovery(nrep = 1, houses = 20, steps = 50, chains = 1,
                   progress = FALSE) |>
    expect_error("[Rr]eplicate 1 at.*20.*households.*missing status")
})

test_that("A failed replicate stops the fits sampled with it", {
  # `record = 5` fails at chain initialization, `record = 0.5` samples
  hh_model <- data.frame(hh_id = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    add_observation("record", model = "uniform", loc = hidden,
                    min = -0.1, max = 0.1)
  hhdata <- hh_model$data_wide
  hhdata_2 <- transform(hhdata, record = 5)
  # the fit of `hhdata` alone samples for minutes
  fit_args <- list(steps = 400000000L, burn = 0, thin = 100000L, adapt = 0,
                   chains = 1, threads = 1L, id_draws = ".none")
  fit_wave <- function(wides, outer, failed) {
    system.time(
      sweep_fit_wave_(hh_model, wides, reps = 1:2, houses = 1L, seeds = 1:2,
                      fit_args = fit_args, outer = outer, bar = NULL) |>
        expect_error(paste0("[Rr]eplicate ", failed,
                            " at.*[Ee]xceeded attempts"))
    )[["elapsed"]]
  }
  # sampled together: the running fit stops at its next iteration
  fit_wave(list(hhdata, hhdata_2), outer = 2L, failed = 2L) |> expect_lt(10)
  # sampled in turn: the next fit never starts
  fit_wave(list(hhdata_2, hhdata), outer = 1L, failed = 1L) |> expect_lt(10)
})

test_that("A run too short to diagnose is flagged whole", {
  expect_warning(
    results <- check_recovery(
      mk_recov$model, nrep = 2, houses = 20, steps = 50, chains = 1,
      progress = FALSE
    ),
    "raised 1 distinct warning"
  )
  results$warnings$message |> expect_match("[Tt]otal sample size")
  results$warnings$count |> expect_identical(2L)
  # diagnosed all the same, and always flagged
  results$estimates$rhat |> expect_no_na()
  results$estimates$ess_bulk |> expect_no_na()
  results$estimates$flagged |> expect_all_true()
  results$overview$n_flagged |> expect_identical(2L)
})

test_that("Recovery leaves the model it was given untouched", {
  hh_model <- mk_recov$model
  params <- get_parameters(hh_model)
  stats::setNames(params$value, params$name)[c("alpha", "beta")] |>
    expect_equal(c(alpha = 0.15, beta = 0.20))
  hh_model$data_wide |> expect_nrow(nrow(flu_1918_frederick))
})

test_that("Every estimated parameter needs a simulation value", {
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = list(value = 0.15, prior = "uniform"),
                  beta = list(prior = "uniform"),
                  gamma = list(value = 1, prior = "uniform", max = 3),
                  hh_size_ref = 4)

  check_recovery(hh_model, nrep = 1, progress = FALSE) |>
    expect_error("[Nn]o simulation value.*beta")
  # a value given for a parameter the model does not carry is reported; the
  # fit is too short to pass the diagnostics
  check_recovery(hh_model, nrep = 1, steps = 50, thin = 1, chains = 1,
                 progress = FALSE, draws = data.frame(beta = 0.2, rho = 1)) |>
    expect_warning("[Nn]o parameter named") |>
    suppressWarnings()
})

test_that("A collapsed baseline takes no simulation value", {
  # `community` names no reference, so `alpha` is the scale of its own group
  expect_warning(
    hh_model <- data.frame(
      hh_id = rep(1:4, each = 3),
      status = c(1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1),
      member = rep(c("child", "adult", "adult"), 4)
    ) |>
      new_housepi_model("final-size") |>
      add_parameter(import = "template") |>
      add_group_risk("community", by = "member"),
    "alpha.*constant reference"
  )

  check_recovery(hh_model, nrep = 1, progress = FALSE,
                 draws = data.frame(alpha = 0.3)) |>
    expect_error("draws.*alpha.*other than 1")
})

test_that("Recovery refuses a sweep it cannot run", {
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = list(value = 0.15, const = TRUE),
                  beta = list(value = 0.20, const = TRUE),
                  gamma = list(value = 1, const = TRUE), hh_size_ref = 4)

  check_recovery(hh_model, nrep = 1, progress = FALSE) |>
    expect_error("[Nn]o parameter to estimate")
  new_housepi_model(model = "final-size") |>
    add_parameter(import = "template") |>
    check_recovery(nrep = 1, progress = FALSE) |>
    expect_error("[Nn]o data is attached")
  # a setting no fit can run is refused up front, naming no replicate
  failed <- check_recovery(mk_recov$model, nrep = 1, houses = 20, steps = 50,
                           thin = 100, chains = 1, progress = FALSE) |>
    expect_error("[Tt]hin.*cannot be greater than.*steps")
  conditionMessage(failed) |> expect_no_match("[Rr]eplicate")
})

test_that("Latent inference recovers the simulating parameter", {
  skip_on_cran()
  # @note test model simulates only on the statuses it is given
  # `rho` mixes as slowly whatever the household count: few households and a
  # long thinned chain keep R-hat clear of its limit
  data.frame(hh_id = seq(50), status = 1) |>
    new_housepi_model("test") |>
    add_parameter("rho", value = 2, prior = "uniform", min = -5, max = 5) |>
    add_latent("hidden", model = "normal", mean = rho, sd = 1) |>
    add_observation("record", model = "normal", mean = hidden, sd = 1) |>
    check_recovery(nrep = 1, steps = 4000, thin = 20, progress = FALSE) |>
    suppressWarnings() |>
    expect_recovered()
})

test_that("Final-size inference recovers the simulating parameters", {
  skip_on_cran()
  flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    check_recovery(nrep = 1, progress = FALSE) |>
    suppressWarnings() |>
    expect_recovered()
})

test_that("Case-ascertained inference recovers the simulating parameters", {
  skip_on_cran()
  # observed onsets date the infection times the transmission reads
  sarscov2_2021_icofs[c("id", "hh_id")] |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_infectivity_profile(import = "template") |>
    add_latent(
      infec_time = list(value = 0, min = start_date - 14, max = end_date,
                        proposal = "rw"),
      incub_period = list(model = "gamma", mean = 2, sd = 1, max = 14),
      onset_time = infec_time + incub_period
    ) |>
    add_observation(
      onset_date = list(discr = TRUE, model = "normal",
                        mean = onset_time, sd = 1, min = floor(onset_time)),
      start_date = list(value = min(onset_date), level = "house"),
      end_date = start_date + 21
    ) |>
    check_recovery(nrep = 1, houses = 200, progress = FALSE) |>
    suppressWarnings() |>
    expect_recovered()
})
