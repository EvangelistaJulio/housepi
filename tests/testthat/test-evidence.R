library(testthat)

test_that("Power-posterior log-evidence matches a closed-form marginal", {
  x <- 2
  a <- 1
  b <- 3
  evidence <- log((x + b) / (x + a)) / (b - a) - x / (b + x) / (a + x)

  hhmodel <- data.frame(house_id = 1, id = 1, status = 1, x = x) |>
    new_housepi_model(model = "blank") |>
    add_parameter("theta", prior = "uniform", min = !!a, max = !!b) |>
    add_augmented("rate", prior = "exponential", rate = theta@P) |>
    add_observation("x", obs = "exponential", rate = rate)

  result <- hhmodel |>
    run_power_posteriors(progress = FALSE) |>
    suppressWarnings() # temporary warning

  expect_s3_class(result, "HousepiPowerPosteriors")
  expect_named(result$detail, c("inv_temper", "mean_log_lik", "var_log_lik"))

  # one ladder, several estimators computed on demand
  expect_equal(evaluate_log_evidence(result, "thermodynamic"), log(evidence),
               tolerance = 0.02)
  expect_equal(evaluate_log_evidence(result, "stepping-stone"), log(evidence),
               tolerance = 0.02)

  # the one-shot model path returns the same value
  hhmodel |>
    evaluate_log_evidence(method = "thermodynamic", progress = FALSE) |>
    suppressWarnings() |>
    expect_equal(log(evidence), tolerance = 0.02)

  # raw per-rung draws are posterior draws: one variable per rung
  expect_true(posterior::is_draws(result$draws))
  expect_equal(posterior::nvariables(result$draws), nrow(result$detail))
  expect_error(result$draws <- 0, "read only")
})
