library(testthat)

####  Definitions  ####

# ranks of four replicates of two parameters, readable by hand
mk_calib_frame <- local({
  estimates <- data.frame(
    houses = 10L, rep = rep(1:4, each = 2L),
    variable = rep(c("rho", "rho_2"), 4L),
    rank = c(0L, 3L, 1L, 3L, 1L, 3L, 3L, 3L),
    flagged = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  )
  list(estimates = estimates)
})

test_that("A rank counts the thinned draws below the truth", {
  calibration_rank_(seq(10), 5.5, 10) |> expect_identical(5L)
  # ten draws evenly thinned out of a hundred: the rank tops at ten
  calibration_rank_(seq(100), 101, 10) |> expect_identical(10L)
  calibration_rank_(seq(100), 0, 10) |> expect_identical(0L)
  # thinned to 1, 12, 23, 34, 45, 56, ...
  calibration_rank_(seq(100), 50.5, 10) |> expect_identical(5L)
})

test_that("Ties with the truth are split at random", {
  # one draw below and three tied: every rank in 1..4 is reachable
  ranks <- replicate(200L, calibration_rank_(c(1, 2, 2, 2, 3), 2, 5))
  ranks |> expect_in(1:4)
  range(ranks) |> expect_identical(c(1L, 4L))
})

test_that("The uniformity statistic separates spread ranks from piled ones", {
  threshold <- calibration_log_threshold_(100, 9, 0.05)
  # every rank of 0..9 ten times: the ECDF of the uniform exactly
  calibration_log_gamma_(rep(0:9, 10L), 9) |> expect_gte(threshold)
  # posterior too narrow: ranks pile at both ends
  calibration_log_gamma_(rep(c(0L, 9L), 50L), 9) |> expect_lt(threshold)
  # posterior biased: ranks pile on one side
  calibration_log_gamma_(rep(0:4, 20L), 9) |> expect_lt(threshold)
})

test_that("The threshold is the largest level whose band holds at 0.95", {
  for (n in c(10, 50, 200)) {
    threshold <- calibration_log_threshold_(n, 99, 0.05)
    # uniform ranks rejected with probability at most 0.05
    calibration_inside_(n, 99, exp(threshold)) |> expect_gte(0.95)
    # any higher level drops the band below 0.95
    calibration_inside_(n, 99, exp(threshold + 1e-8)) |> expect_lt(0.95)
  }
})

test_that("Uniform ranks are rejected at the level of the test", {
  skip_on_cran()
  threshold <- calibration_log_threshold_(50, 99, 0.05)
  log_gamma <- replicate(2000L, {
    calibration_log_gamma_(sample.int(100L, 50L, replace = TRUE) - 1L, 99)
  })
  # four Monte Carlo standard errors around 0.05
  mean(log_gamma < threshold) |> expect_equal(0.05, tolerance = 0.02 / 0.05)
})

test_that("The ECDF of the ranks reads with its band", {
  threshold <- \(n) log(0.5)
  ecdf <- calibration_ecdf_(mk_calib_frame$estimates, 3, threshold)
  ecdf$variable |> expect_identical(rep(c("rho", "rho_2"), each = 3L))
  ecdf$z |> expect_equal(rep(1:3 / 4, 2L))
  # rho ranks 0, 1, 1, 3; rho_2 ranks all 3
  ecdf$ecdf |> expect_equal(c(0.25, 0.75, 0.75, 0, 0, 0))
  # pointwise band at level 0.5 for four ranks
  ecdf$lower |> expect_equal(rep(stats::qbinom(0.25, 4, 1:3 / 4) / 4, 2L))
  ecdf$upper |> expect_equal(rep(stats::qbinom(0.75, 4, 1:3 / 4) / 4, 2L))
})

test_that("The calibration summary reads one verdict per parameter", {
  threshold <- \(n) calibration_log_threshold_(n, 3, 0.05)
  summary <- calibration_summary_(mk_calib_frame$estimates, 3, threshold)
  summary$variable |> expect_identical(c("rho", "rho_2"))
  summary$nrep |> expect_identical(c(4L, 4L))
  summary$n_flagged |> expect_identical(c(1L, 0L))
  summary$log_gamma_threshold |> expect_equal(rep(threshold(4), 2L))
  # every rank of rho_2 at the top
  summary$uniform |> expect_identical(c(TRUE, FALSE))
})

####  End-to-end  ####

# a sweep short enough for the suite, every fit too short to pass the
# diagnostics, and the attached households are surveyed as they are
mk_calib <- local({
  hh_model <- data.frame(hh_id = seq(30), status = 1) |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform", min = -5, max = 5) |>
    add_latent("hidden", model = "normal", mean = rho, sd = 1) |>
    add_observation("record", model = "normal", mean = hidden, sd = 1)
  calibration <- hh_model |>
    check_calibration(nrep = 3, structure_sampler = NULL, rank_draws = 9,
                      steps = 400, thin = 1, chains = 1, progress = FALSE) |>
    suppressWarnings()
  list(model = hh_model, calibration = calibration)
})

test_that("Calibration calibration are read only", {
  mk_calib$calibration |>
    expect_read_only_field("estimates", "replicates",
                           "overview", "houses", "nrep", "rank_draws",
                           "warnings")
})

test_that("Every replicate is ranked at a value drawn from the prior", {
  calibration <- mk_calib$calibration
  estimates <- calibration$estimates
  estimates |> expect_nrow(3L)
  estimates$rep |> expect_identical(1:3)
  estimates$rank |> expect_in(0:9)
  # a fresh prior draw per replicate, not the pin
  unique(estimates$truth) |> expect_length(3L)
  (estimates$truth > -5 & estimates$truth < 5 & estimates$truth != 2) |>
    expect_all_true()
  calibration$evaluate_ecdf() |> expect_nrow(9L)
  calibration$summarise()$nrep |> expect_identical(3L)
  calibration$summarize() |> expect_identical(calibration$summarise())
  calibration$overview$nrep |> expect_identical(3L)
  calibration$houses |> expect_identical(30L)
  calibration$rank_draws |> expect_identical(9L)
})

test_that("Calibration calibration print their settings", {
  calibration <- mk_calib$calibration
  lines <- helper_print_lines(printed <- print(calibration))
  # printing returns the object itself
  printed |> expect_identical(calibration)
  lines[1] |> expect_match("[Ss]imulation calibration check.*\"test\"")
  lines |> expect_contains(c("Households: 30", "Replicates: 3",
                             "Posterior draws ranked among: 9"))
  # the count is that of the flagged replicates
  lines |>
    expect_contains(paste("Flagged by diagnostics:",
                          sum(calibration$replicates$flagged)))
  # durations read as the progress bar prints them
  lines |> expect_match("^Run time: [0-9]+s$", all = FALSE)
  lines |> expect_match("^Inference time per fit: [0-9]+s$", all = FALSE)
})

test_that("A sweep reproduces from a seed at any thread count", {
  sweep_at <- function(threads) {
    withr::with_seed(1L, check_calibration(
      mk_calib$model, nrep = 3, rank_draws = 9, steps = 400, thin = 1,
      chains = 1, threads = threads, progress = FALSE
    )) |>
      suppressWarnings()
  }
  # one fit per wave, then a single wave of three: the ranks break ties
  # from the same stream as the engine seeds
  calibration <- sweep_at(1L)
  calibration_2 <- sweep_at(3L)
  untimed <- \(replicates) replicates[names(replicates) != "inference_secs"]
  calibration_2$estimates |> expect_identical(calibration$estimates)
  untimed(calibration_2$replicates) |>
    expect_identical(untimed(calibration$replicates))
  calibration_2$warnings |> expect_identical(calibration$warnings)
})

test_that("Calibration refuses a sweep it cannot rank", {
  hh_model <- mk_calib$model
  check_calibration(hh_model, steps = 100L, chains = 1L, progress = FALSE) |>
    expect_error("[Ss]tores 10 draws.*fewer than.*rank_draws")
  check_calibration(hh_model, rank_draws = 0, progress = FALSE) |>
    expect_error("[Nn]umber of ranked draws")
  # one household count per check
  check_calibration(hh_model, houses = c(10, 20), progress = FALSE) |>
    expect_error("[Nn]umber of households")
  data.frame(hh_id = seq(30), status = 1) |>
    new_housepi_model("test") |>
    add_parameter("rho", const = 2) |>
    check_calibration(progress = FALSE) |>
    expect_error("[Nn]o parameter to estimate")
})

test_that("Calibration checks of one model pool their replicates", {
  rerun <- function(nrep = 2, rank_draws = 9, structure_sampler = NULL, ...) {
    check_calibration(mk_calib$model, nrep = nrep,
                      structure_sampler = structure_sampler,
                      rank_draws = rank_draws, steps = 400, thin = 1,
                      chains = 1, progress = FALSE, ...) |>
      suppressWarnings()
  }
  calibration <- mk_calib$calibration
  combined <- c(calibration, rerun())

  combined$nrep |> expect_identical(5L)
  combined$estimates$rep |> expect_identical(1:5)
  # the uniformity test reads the pooled replicates
  combined$summarise()$nrep |> expect_identical(5L)
  combined$summarise()$log_gamma_threshold |>
    expect_equal(calibration_log_threshold_(5, 9, 0.05))
  combined$evaluate_ecdf() |> expect_nrow(9L)
  calibration$nrep |> expect_identical(3L)
  c(calibration, calibration) |> expect_error("[Rr]epeat a replicate.*30")
  c(calibration, rerun(nrep = 1, rank_draws = 8)) |>
    expect_error("[Dd]iffer in.*rank_draws")

  # one count per check, several once pooled
  resampled <- \(nrep, houses) {
    rerun(nrep, structure_sampler = sample_house_structure, houses = houses)
  }
  calibration_2 <- resampled(2, 30)
  combined <- c(calibration_2, resampled(1, 40), resampled(1, 30))
  combined$houses |> expect_identical(c(30L, 40L))
  combined$nrep |> expect_identical(c(3L, 1L))
  # the third check follows the first at their shared count
  with(combined$replicates, rep[houses == 30L]) |> expect_identical(1:3)
  # each count is tested at its own number of replicates
  combined$summarise()$log_gamma_threshold |>
    expect_equal(c(calibration_log_threshold_(3, 9, 0.05),
                   calibration_log_threshold_(1, 9, 0.05)))
  combined$evaluate_ecdf() |> expect_nrow(18L)
  c(calibration, calibration_2) |>
    expect_error("[Dd]iffer in.*structure_sampler")
})

test_that("Latent inference is calibrated against the prior", {
  skip_on_cran()
  # a replicate flagged by convergence diagnostics warns
  calibration <- mk_calib$model |>
    check_calibration(threads = 4, progress = FALSE) |>
    suppressWarnings()
  # a test at level 0.001, so a calibrated sampler fails one run in 1000
  nrep <- calibration$summarise()$nrep
  threshold <- calibration_log_threshold_(nrep, 99, 0.001)
  (calibration$summarise()$log_gamma >= threshold) |> expect_all_true()
})
