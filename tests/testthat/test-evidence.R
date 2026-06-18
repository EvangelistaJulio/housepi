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
    evaluate_log_evidence(progress = FALSE) |>
    suppressWarnings() # temporary warning

  expect_equal(result$estimate, log(evidence), tolerance = 0.02)
  expect_s3_class(result, "HousepiEvidence")
  expect_named(result$detail, c("inv_temper", "mean_log_lik", "var_log_lik"))
  expect_error(result$estimate <- 0, "read only")

  # print is a compact summary, not the per-node ladder
  expect_message(print(result), "Log-evidence")
})
