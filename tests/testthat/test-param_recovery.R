library(testthat)

######################################
###  PARAM RECOVERY : DEFINITIONS  ###
######################################

# households of two sizes, the covariate constant within each one
mk_recov_data <- local({
  hhdata <- data.frame(hhid = c("h1", rep("h2", 3)),
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

test_that("Default household resampling behaves as expected", {
  sample_house_structure(data.frame(datum = 1)) |>
    expect_error("[Cc]olumn.*hhid.*is missing")
  sample_house_structure(mk_recov_data$hhdata, houses = 0) |>
    expect_error("[Nn]umber of households")
  # the sampler draws one structure, so a sweep of counts is refused here
  sample_house_structure(mk_recov_data$hhdata, houses = c(2, 3)) |>
    expect_error("[Nn]umber of households")

  # default sampling targets the observed household count
  struct <- sample_house_structure(mk_recov_data$hhdata)
  length(unique(struct$hhid)) |> expect_identical(2L)
})

test_that("Households are resampled whole", {
  struct <- sample_house_structure(mk_recov_data$hhdata, houses = 6)
  # every household is renumbered
  sort(unique(struct$hhid)) |> expect_identical(1:6)
  # an identifier carried by the input is regenerated
  struct$id |> expect_identical(seq_len(nrow(struct)))
  # household sizes are only the observed ones
  sizes <- as.integer(table(struct$hhid))
  sizes |> expect_in(c(1L, 3L))
  # covariate datum is still unique per household and not mixed
  levels <- tapply(struct$datum, struct$hhid, \(x) length(unique(x)))
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
    inference_secs = c(1, 3, NA), min_ess_bulk = c(100, 600, NA),
    flagged = c(FALSE, TRUE, NA), error = c(NA, NA, "crashed")
  )
  overview <- recovery_overview_(replicates)
  overview$nrep |> expect_identical(c(2L, 1L))
  overview$n_crashed |> expect_identical(c(0L, 1L))
  overview$n_flagged |> expect_identical(c(1L, 0L))
  # a crashed replicate carries no cost to average
  overview$mean_inference_secs |> expect_equal(c(2, NaN))
  # pooled over the replicates: (1 + 3) / (100 + 600)
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
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  results <- hhmodel |>
    run_parameter_recovery_check(nrep = 2, houses = c(20, 40),
                                 steps = 50, chains = 1, progress = FALSE) |>
    suppressWarnings()
  list(model = hhmodel, results = results)
})

# no fit is flagged and every point estimate lies within four posterior sd
# of its simulating value, which a calibrated posterior strays past about once
# in 16,000 fits
expect_recovered <- function(results) {
  estimates <- results$estimates
  estimates$flagged |> expect_all_false()
  (abs(estimates$estimate - estimates$truth) < 4 * estimates$sd) |>
    expect_all_true()
}

test_that("Recovery results are read only", {
  mk_recov$results |>
    helper_test_read_only("estimates", "replicates", "summary",
                          "overview", "truth", "houses", "nrep",
                          "point_estimate", "cri_level", "warnings")
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
  results$overview$n_crashed |> expect_identical(c(0L, 0L))
  results$nrep |> expect_identical(2L)
  results$point_estimate |> expect_identical("median")
  results$cri_level |> expect_equal(0.95)
})

test_that("Each replicate reports its cost and its survey", {
  replicates <- mk_recov$results$replicates
  (replicates$inference_secs > 0) |> expect_all_true()
  estimates <- mk_recov$results$estimates
  replicates$min_ess_bulk |>
    expect_equal(as.vector(tapply(estimates$ess_bulk,
                                  estimates[c("rep", "houses")], min)))
  (mk_recov$results$overview$secs_per_ess > 0) |> expect_all_true()
  replicates$error |> expect_all_na()
  estimates$covers |>
    expect_identical(estimates$cri_lower <= estimates$truth &
                       estimates$truth <= estimates$cri_upper)
  mk_recov$results$summary$coverage |> expect_in(c(0, 0.5, 1))
})

test_that("A run too short to diagnose is flagged whole", {
  expect_warning(
    results <- run_parameter_recovery_check(
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
  hhmodel <- mk_recov$model
  params <- get_parameters(hhmodel)
  stats::setNames(params$value, params$name)[c("alpha", "beta")] |>
    expect_equal(c(alpha = 0.15, beta = 0.20))
  hhmodel$data_wide |> expect_nrow(nrow(flu_1918_frederick))
})

test_that("Every estimated parameter needs a simulation value", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = list(value = 0.15, prior = "uniform"),
                  beta = list(prior = "uniform"),
                  gamma = list(value = 1, prior = "uniform", max = 3),
                  hhsize_ref = 4)

  run_parameter_recovery_check(hhmodel, nrep = 1, progress = FALSE) |>
    expect_error("[Nn]o simulation value.*beta")
  # a value given for a parameter the model does not carry is reported
  run_parameter_recovery_check(hhmodel, nrep = 1, progress = FALSE,
                               truth = data.frame(beta = 0.2, rho = 1)) |>
    expect_warning("[Nn]o parameter named")
})

test_that("A collapsed baseline takes no simulation value", {
  # `community` names no reference, so `alpha` is the scale of its own group
  expect_warning(
    hhmodel <- data.frame(
      hhid = rep(1:4, each = 3),
      status = c(1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1),
      member = rep(c("child", "adult", "adult"), 4)
    ) |>
      new_housepi_model("final-size") |>
      add_parameter(import = "template") |>
      add_group_risk("community", by = "member"),
    "alpha.*constant reference"
  )

  run_parameter_recovery_check(hhmodel, nrep = 1, progress = FALSE,
                               truth = data.frame(alpha = 0.3)) |>
    expect_error("[Tt]ruth.*alpha.*other than 1")
})

test_that("Recovery needs a parameter to recover", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = list(value = 0.15, const = TRUE),
                  beta = list(value = 0.20, const = TRUE),
                  gamma = list(value = 1, const = TRUE), hhsize_ref = 4)

  run_parameter_recovery_check(hhmodel, nrep = 1, progress = FALSE) |>
    expect_error("[Nn]o parameter to recover")
})

test_that("Latent inference recovers the simulating parameter", {
  skip_on_cran()
  # @note test model simulates only on the statuses it is given
  data.frame(hhid = seq(200), status = 1) |>
    new_housepi_model("test") |>
    add_parameter("rho", value = 2, prior = "uniform", min = -5, max = 5) |>
    add_latent("hidden", model = "normal", mean = rho@P, sd = 1) |>
    add_observation("record", model = "normal", mean = hidden@A, sd = 1) |>
    run_parameter_recovery_check(nrep = 1, progress = FALSE) |>
    expect_recovered()
})

test_that("Final-size inference recovers the simulating parameters", {
  skip_on_cran()
  flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_parameter_recovery_check(nrep = 1, progress = FALSE) |>
    expect_recovered()
})

test_that("Case-ascertained inference recovers the simulating parameters", {
  skip_on_cran()
  # observed onsets date the infection times the transmission reads
  sarscov2_2021_icofs[c("id", "hhid")] |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_infectivity_profile(import = "template") |>
    add_latent(
      infec_time = list(value = 0, model = "uniform",
                        min = start_date - 14, max = end_date),
      incub_period = list(model = "gamma", mean = 2, sd = 1, max = 14),
      onset_time = infec_time + incub_period
    ) |>
    add_observation(
      onset_date = list(discr = TRUE, model = "normal",
                        mean = onset_time, sd = 1, min = floor(onset_time)),
      start_date = list(value = min(onset_date), level = "house"),
      end_date = start_date + 21
    ) |>
    run_parameter_recovery_check(nrep = 1, houses = 200, progress = FALSE) |>
    expect_recovered()
})
