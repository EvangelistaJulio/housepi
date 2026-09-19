library(testthat)

# Tempered-indicator regime: a hard, parameter-dependent observation support.
#   rho ~ Uniform(0, 1);  record ~ Uniform(rho - c, rho + c), observed 0.5.
# The likelihood is a pure indicator, so log Z = 0 for a window inside [0, 1].
# At inverse temperature 0 the tempered likelihood is 0 even where it is -Inf:
# the chain must sample the whole prior, not stay confined to [0.05, 0.95].
mk_indicator <- local({
  hhmodel <- data.frame(hhid = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform", min = 0, max = 1) |>
    add_observation(
      "record", model = "uniform", loc = rho@P, min = - 0.45, max = 0.45
    )
  list(model = hhmodel)
})

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

# retry until a run starts (a random prior init may land inadmissible)
retry_start <- function(run, hhmodel, ...) {
  for (attempt in 1:50) {
    results <- try(suppressWarnings(run(hhmodel, ...)), silent = TRUE)
    if (!inherits(results, "try-error")) return(results)
  }
  NULL
}
start_ladder <- function(hhmodel, ...) {
  retry_start(run_power_posteriors, hhmodel, ...)
}

###########################
## EVIDENCE : END-TO-END ##
###########################

test_that("Inverse temperature 0 samples the full prior with an indicator", {
  # FIXME: time sensitive test
  rho <- mk_indicator$model |>
    retry_start(run = run_inference, inv_temp = 0, chains = 1L) |>
    get_parameter_draws() |>
    posterior::extract_variable("rho")
  # chain must escape the observation window [0.05, 0.95]
  expect_gt(max(rho), 0.95)
  expect_lt(min(rho), 0.05)
  expect_gt(mean(abs(rho - 0.5) >= 0.45), 0.05)
})

test_that("Power-posterior log-evidence matches a closed-form marginal", {
  x <- 2
  a <- 1
  b <- 3
  # analytical value
  evidence <- log((x + b) / (x + a)) / (b - a) - x / (b + x) / (a + x)

  hhmodel <- data.frame(hhid = "h1", status = 1, record = x) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform", min = !!a, max = !!b) |>
    add_latent("z", model = "exponential", rate = theta@P) |>
    add_observation("record", model = "exponential", rate = z)

  result <- hhmodel |>
    run_power_posteriors(steps = 4e3L) |>
    suppressWarnings() # once per session temporary warning

  expect_s3_class(result, "HousepiPowerPosteriors")
  expect_snapshot(print(result))
  expect_named(result$detail, c("inv_temp", "mean_log_lik", "var_log_lik"))

  # raw per-rung draws are 'posterior' draws
  expect_true(posterior::is_draws(result$draws))
  expect_equal(posterior::nvariables(result$draws), nrow(result$detail))
  expect_error(result$draws <- 0, "read only")

  # stable estimate for chosen values
  exp_log_lik <- vapply(
    result$detail$inv_temp,
    function(beta) expec_loglik_powerpost(x, beta, a, b), numeric(1)
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
  mk_indicator$model |>
    run_power_posteriors(inv_temp = 0.5) |>
    expect_error("inv_temp")
  # a ladder holds at least two rungs
  mk_indicator$model |>
    run_power_posteriors(rungs = 1L) |>
    suppressWarnings() |>
    expect_error("[Rr]ungs was not greater than 1")
})

test_that("Stepping-stone is unbiased, thermodynamic refuses a -Inf ladder", {
  results <- start_ladder(mk_indicator$model, chains = 1L, rungs = 12L)

  # the true log-evidence is 0
  evaluate_log_evidence(results, "stepping-stone")[["estimate"]] |>
    abs() |>
    expect_lt(0.08)

  # the beta = 0 rung legitimately records -Inf draws
  expect_true(any(is.infinite(
    posterior::extract_variable(results$draws, "log_lik[1]")
  )))
  # thermodynamic integration is invalid there: error pointing at stepping-stone
  evaluate_log_evidence(results, "thermodynamic") |>
    expect_error("stepping-stone")
})

test_that("A tempered indicator with a latent runs at beta = 0 without NaN", {
  results <- data.frame(hhid = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform", min = 0, max = 1) |>
    add_latent("hidden", model = "uniform", min = 0, max = 1) |>
    add_observation("record", model = "uniform",
                    min = hidden - 0.4, max = hidden + 0.4) |>
    start_ladder(chains = 1L, rungs = 8L)
  expect_false(is.null(results))

  # incremental totals are bypassed at beta = 0, so -Inf states do not throw
  expect_true(any(is.infinite(
    posterior::extract_variable(results$draws, "log_lik[1]")
  )))
  evaluate_log_evidence(results, "stepping-stone")[["estimate"]] |>
    is.finite() |>
    expect_true()
})
