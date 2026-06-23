library(testthat)

# analytical log-likelihood under power-posterior
expec_loglik_powerpost <- function(x, beta, a, b) {
  stopifnot(b > a, a >= 0.)

  ua <- a + beta * x
  ub <- b + beta * x
  const <- gamma(beta + 1) / (b - a)

  if (isTRUE(beta == 1.)) {
    A <- function(u) log(u) + x / u
    R <- function(u) 0.5 * log(u)^2 + x / u - x^2 / u^2
  } else {
    A <- function(u) u^(1 - beta) / (1 - beta) + x * u^(-beta)
    R <- function(u) {
      u^(1 - beta) / (1 - beta)^2 + x * u^(-beta) - beta * x^2 * u^(-beta - 1)
    }
  }

  Q <- const * (A(ub) - A(ua)) # normalization
  bracket <- (log(ub) * A(ub) - R(ub)) - (log(ua) * A(ua) - R(ua))

  digamma(beta + 1) - const * bracket / Q
}

test_that("Power-posterior log-evidence matches a closed-form marginal", {
  x <- 2
  a <- 1
  b <- 3
  # analytical value
  evidence <- log((x + b) / (x + a)) / (b - a) - x / (b + x) / (a + x)
  log_evidence <- log(evidence)

  hhmodel <- data.frame(house_id = "h", status = 1, x = x) |>
    new_housepi_model("blank") |>
    add_parameter("theta", prior = "uniform", min = !!a, max = !!b) |>
    add_augmented("z", prior = "exponential", rate = theta@P) |>
    add_observation("x", obs = "exponential", rate = z)

  result <- hhmodel |>
    run_power_posteriors(progress = FALSE, steps = 5e3) |>
    suppressWarnings() # temporary warning

  result |>
    expect_s3_class("HousepiPowerPosteriors")
  result$detail |>
    expect_named(c("inv_temp", "mean_log_lik", "var_log_lik"))

  # raw per-rung draws are 'posterior' draws
  expect_true(posterior::is_draws(result$draws))
  posterior::nvariables(result$draws) |>
    expect_equal(nrow(result$detail))
  expect_error(result$draws <- 0, "read only")

  # stable estimate for chosen values
  exp_log_lik <- sapply(
    result$detail$inv_temp,
    function(beta) expec_loglik_powerpost(x, beta, a, b)
  )
  expect_true(all(
    abs(result$detail$mean_log_lik - exp_log_lik) /
      sqrt(result$detail$var_log_lik) <= 0.05
  ))

  # estimators computed on demand
  evaluate_log_evidence(result, "thermodynamic") |>
    expect_equal(log(evidence), tolerance = 0.02)
  evaluate_log_evidence(result, "stepping-stone") |>
    expect_equal(log(evidence), tolerance = 0.02)

  # one-shot model path returns same value
  hhmodel |>
    evaluate_log_evidence("thermodynamic", progress = FALSE) |>
    suppressWarnings() |>
    expect_equal(log(evidence), tolerance = 0.02)

})
