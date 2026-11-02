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
# records of their own feed the fit the posterior check is conditioned on
mk_calib <- local({
  hh_model <- data.frame(hh_id = seq(30), status = 1,
                         record = stats::rnorm(30, 1, sqrt(2))) |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform", min = -5, max = 5) |>
    add_latent("hidden", model = "normal", mean = rho, sd = 1) |>
    add_observation("record", model = "normal", mean = hidden, sd = 1)
  prior_calibration <- hh_model |>
    check_calibration(nrep = 3, structure_sampler = NULL, rank_levels = 10,
                      steps = 400, thin = 1, chains = 1, progress = FALSE) |>
    suppressWarnings()
  hh_fit <- hh_model |>
    run_inference(steps = 1000, thin = 1, chains = 1, id_draws = ".none",
                  progress = FALSE) |>
    suppressWarnings()
  post_calibration <- hh_model |>
    check_calibration(hh_fit, nrep = 3, structure_sampler = NULL,
                      rank_levels = 10, steps = 400, thin = 1, chains = 1,
                      progress = FALSE) |>
    suppressWarnings()
  list(model = hh_model, fit = hh_fit,
       prior_calib = prior_calibration,
       post_calib = post_calibration)
})

test_that("Calibration calibration are read only", {
  mk_calib$prior_calib |>
    expect_read_only_field("estimates", "replicates",
                           "overview", "houses", "nrep", "source",
                           "rank_levels", "warnings")
})

test_that("Every replicate is ranked at a value drawn from the prior", {
  calibration <- mk_calib$prior_calib
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
  calibration$rank_levels |> expect_identical(10L)
  calibration$source |> expect_identical("prior")
})

test_that("Every replicate is ranked at a distinct posterior draw of the run", {
  calibration <- mk_calib$post_calib
  estimates <- calibration$estimates
  calibration$source |> expect_identical("posterior")
  estimates |> expect_nrow(3L)
  estimates$rank |> expect_in(0:9)
  unique(estimates$truth) |> expect_length(3L)
  estimates$truth |>
    expect_in(as.numeric(get_parameter_draws(mk_calib$fit)$rho))
  calibration$summarise()$nrep |> expect_identical(3L)
})

test_that("Calibration calibration print their settings", {
  calibration <- mk_calib$prior_calib
  lines <- helper_print_lines(printed <- print(calibration))
  # printing returns the object itself
  printed |> expect_identical(calibration)
  lines[1] |> expect_match("[Pp]rior calibration check.*\"test\"")
  lines |> expect_contains(c("Households: 30", "Replicates: 3",
                             "Rank levels: 10"))
  lines |> expect_no_match("Conditioned on")
  # a posterior check names the households it is conditioned on
  lines <- helper_print_lines(print(mk_calib$post_calib))
  lines[1] |> expect_match("[Pp]osterior calibration check.*\"test\"")
  lines |> expect_contains("Conditioned on: 30 observed households")
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
      mk_calib$model, nrep = 3, rank_levels = 10, steps = 400, thin = 1,
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
  hh_model |>
    check_calibration(steps = 100L, chains = 1L, progress = FALSE) |>
    expect_error("[Ss]tores 10 draws.*too few to rank among 99")
  hh_model |>
    check_calibration(rank_levels = 0, progress = FALSE) |>
    expect_error("[Nn]umber of rank levels")
  # a rank needs one posterior draw at least
  hh_model |>
    check_calibration(rank_levels = 1, progress = FALSE) |>
    expect_error("[Nn]umber of rank levels.*at least 2")
  # one household count per check
  hh_model |>
    check_calibration(houses = c(10, 20), progress = FALSE) |>
    expect_error("[Nn]umber of households")
  data.frame(hh_id = seq(30), status = 1) |>
    new_housepi_model("test") |>
    add_parameter("rho", const = 2) |>
    check_calibration(progress = FALSE) |>
    expect_error("[Nn]o parameter to estimate")
})

test_that("Posterior calibration refuses draws it cannot condition on", {
  hh_model <- mk_calib$model
  hh_fit <- mk_calib$fit
  hh_model |>
    check_calibration(get_parameter_draws(hh_fit), progress = FALSE) |>
    expect_error("[Dd]raws.*should come from.*run_inference")
  hh_model |>
    check_calibration(data.frame(rho = 1:3), progress = FALSE) |>
    expect_error("[Dd]raws.*should come from.*run_inference")
  # same definitions, other records
  get_dataset(hh_model) |>
    transform(record = record + 1) |>
    add_dataset(hh_model, data = _) |>
    check_calibration(hh_fit, progress = FALSE) |>
    expect_error("[Nn]ot fitted to this model.*[Dd]iffer in.*dataset")
  # a run of 1000 draws holds at most 1000 independent ones
  check_calibration(hh_model, hh_fit, nrep = 1001, progress = FALSE) |>
    expect_error("[Ii]ndependent draws, fewer than.*nrep")
})

test_that("Posterior calibration checks pool across runs only", {
  hh_model <- mk_calib$model
  calibration <- mk_calib$post_calib
  rerun <- function(hh_fit) {
    hh_model |>
      check_calibration(hh_fit, nrep = 2, structure_sampler = NULL,
                        rank_levels = 10, steps = 400, thin = 1, chains = 1,
                        progress = FALSE) |>
      suppressWarnings()
  }
  hh_fit_2 <- hh_model |>
    run_inference(steps = 1000, thin = 1, chains = 1,
                  id_draws = ".none", progress = FALSE) |>
    suppressWarnings()
  combined <- c(calibration, rerun(hh_fit_2))
  combined$nrep |> expect_identical(5L)
  combined$source |> expect_identical("posterior")
  # a run pooled into another is still the same run
  c(calibration, rerun(c(mk_calib$fit, hh_fit_2))) |>
    expect_error("[Rr]ank draws of the same run")
  c(calibration, calibration) |>
    expect_error("[Rr]ank draws of the same run")
  c(calibration, mk_calib$prior_calib) |>
    expect_error("[Dd]iffer in.*source")
})

test_that("Calibration checks of one model pool their replicates", {
  rerun <- function(nrep = 2, rank_levels = 10, structure_sampler = NULL,
                    ...) {
    mk_calib$model |>
      check_calibration(nrep = nrep,
                        structure_sampler = structure_sampler,
                        rank_levels = rank_levels, steps = 400, thin = 1,
                        chains = 1, progress = FALSE, ...) |>
      suppressWarnings()
  }
  calibration <- mk_calib$prior_calib
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
  c(calibration, rerun(nrep = 1, rank_levels = 9)) |>
    expect_error("[Dd]iffer in.*rank_levels")

  # one count per check, several once pooled
  resampled <- function(nrep, houses) {
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

test_that("Latent inference is calibrated on the data", {
  skip_on_cran()
  # records far from a prior they outweigh: a fit of the few simulated
  # households alone would shrink to the prior, and the ranks pile high
  hh_model <- data.frame(hh_id = seq(30), status = 1,
                         record = stats::rnorm(30, 3, sqrt(2))) |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "normal", mean = 0, sd = 1) |>
    add_latent("hidden", model = "normal", mean = rho, sd = 1) |>
    add_observation("record", model = "normal", mean = hidden, sd = 1)
  hh_fit <- hh_model |>
    run_inference(threads = 4, progress = FALSE) |>
    suppressWarnings()
  calibration <- hh_model |>
    check_calibration(hh_fit, nrep = 100, houses = 5, threads = 4,
                      progress = FALSE) |>
    suppressWarnings()
  # a test at level 0.001, so a calibrated sampler fails one run in 1000
  threshold <- calibration_log_threshold_(100, 99, 0.001)
  calibration$summarise()$log_gamma |> expect_gte(threshold)
})
