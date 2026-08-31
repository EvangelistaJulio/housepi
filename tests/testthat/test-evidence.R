library(testthat)

################################################################################
###################################  HELPERS  ##################################
################################################################################

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

# Tempered-indicator regime: a hard, parameter-dependent observation support.
#   theta ~ Uniform(0, 1);  y ~ Uniform(theta - c, theta + c), observed y0.
# The likelihood is a pure indicator, so log Z = 0 for a window inside [0, 1].
# At inverse temperature 0 the tempered likelihood is 0 even where it is -Inf:
# the chain must sample the whole prior, not stay confined to [y0 - c, y0 + c].
mk_indicator_model_ <- function() {
  data.frame(house_id = "a", status = 1, y = 0.5) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform", min = 0, max = 1) |>
    add_observation(
      "y", model = "uniform", loc = theta@P, min = - 0.45, max = 0.45
    )
}

# retry seeds until a ladder starts (a random prior init may land inadmissible)
start_ladder_ <- function(make_model, ...) {
  for (s in 1:50) {
    set.seed(s)
    res <- try(suppressWarnings(
      run_power_posteriors(make_model(), progress = FALSE, ...)
    ), silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
  }
  NULL
}

################################################################################
####################################  TESTS  ###################################
################################################################################

test_that("inverse temperature 0 samples the full prior also with indicator", {
  # FIXME: time sensitive test
  fit <- mk_indicator_model_() |>
    run_inference(inv_temp = 0, chains = 1L) |>
    suppressWarnings()
  theta <- posterior::extract_variable(get_parameter_draws(fit), "theta")
  # chain must escape the observation window [0.05, 0.95]
  expect_gt(max(theta), 0.95)
  expect_lt(min(theta), 0.05)
  expect_gt(mean(abs(theta - 0.5) >= 0.45), 0.05)
})

test_that("Power-posterior log-evidence matches a closed-form marginal", {
  x <- 2
  a <- 1
  b <- 3
  # analytical value
  evidence <- log((x + b) / (x + a)) / (b - a) - x / (b + x) / (a + x)
  log_evidence <- log(evidence)

  hhmodel <- data.frame(house_id = "h", status = 1, x = x) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform", min = !!a, max = !!b) |>
    add_augmented("z", model = "exponential", rate = theta@P) |>
    add_observation("x", model = "exponential", rate = z)

  result <- hhmodel |>
    run_power_posteriors(steps = 4e3L, progress = FALSE) |>
    suppressWarnings() # once per session temporary warning

  result |>
    expect_s3_class("HousepiPowerPosteriors")
  expect_snapshot(print(result))
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
  # NOTE: 10% should be large enough to counter test stochasticity
  expect_true(all(
    abs(result$detail$mean_log_lik - exp_log_lik) /
      sqrt(result$detail$var_log_lik) < 0.1
  ))

  # estimators computed on demand
  evaluate_log_evidence(result, "thermodynamic")[["estimate"]] |>
    expect_equal(log(evidence), tolerance = 0.02)
  evaluate_log_evidence(result, "stepping-stone")[["estimate"]] |>
    expect_equal(log(evidence), tolerance = 0.02)

  # a model must be sampled into a ladder first: no one-shot path
  hhmodel |>
    evaluate_log_evidence() |>
    expect_error("run_power_posteriors")
})

test_that("Ladder sets its own inverse temperature, so it takes none", {
  mk_indicator_model_() |>
    run_power_posteriors(inv_temp = 0.5) |>
    expect_error("inv_temp")
})

test_that("stepping-stone is unbiased and thermodynamic refuses a -Inf ladder", {
  res <- start_ladder_(function() mk_indicator_model_(),
                       chains = 1L, rungs = 12L)

  # the true log-evidence is 0; the old confinement bug returned +0.105
  evaluate_log_evidence(res, "stepping-stone")[["estimate"]] |>
    abs() |>
    expect_lt(0.08)

  # the beta = 0 rung legitimately records -Inf draws
  expect_true(any(is.infinite(
    posterior::extract_variable(res$draws, "log_lik[1]")
  )))
  # thermodynamic integration is invalid there: error pointing at stepping-stone
  evaluate_log_evidence(res, "thermodynamic") |>
    expect_error("stepping-stone")
})

test_that("tempered indicator with an augmented variable runs at beta=0 without NaN", {
  x0 <- 0.5
  cc <- 0.4
  res <- start_ladder_(function() {
    data.frame(house_id = "liquid", status = 1, x = x0) |>
      new_housepi_model("test") |>
      add_parameter("theta", prior = "uniform", min = 0, max = 1) |>
      add_augmented("z", model = "uniform", min = 0, max = 1) |>
      add_observation("x", model = "uniform", min = z - !!cc, max = z + !!cc)
  }, chains = 1L, rungs = 8L)
  expect_false(is.null(res))

  # incremental totals are bypassed at beta = 0, so -Inf states do not throw
  expect_true(any(is.infinite(
    posterior::extract_variable(res$draws, "log_lik[1]")
  )))
  evaluate_log_evidence(res, "stepping-stone")[["estimate"]] |>
    is.finite() |>
    expect_true()
})
