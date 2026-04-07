library(testthat)

trapezoid <- function(x, y) {
  0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))
}

test_that("Expose function has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  fun <- expose_function("sin")
  rm(fun)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in cli::cli_progress_along(seq_len(200))) {
    fun <- expose_function("sin")
    rm(fun)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})

test_that("Function exported from housepi are correct", {
  x <- seq(0, 1000, 1)

  expect_error(expose_function("pow")$func(x), regexp = "Missing non default")
  expect_error(expose_function("pow", 5), regexp = "Unnamed argument")

  expect_equal(expose_function("pow", x2 = -3)$func(x), x^(-3))
  expect_equal(expose_function("pow", x2 = 0.5)$func(x),
               expose_function("sqrt")$func(x))

  expect_equal(expose_function("sqrt")$prim(x),
               expose_function("pow", x2 = 1.5)$func(x) / 1.5)

  expect_equal(expose_function("exp")$func(x), exp(x))
  expect_equal(expose_function("exp")$prim(x), exp(x))

  expect_equal(expose_function("log")$func(x), log(x))
  expect_equal(expose_function("log")$prim(x), x * log(x) - x)

  expect_equal(expose_function("sin")$func(x), sin(x))
  expect_equal(expose_function("sin")$prim(x),
               -expose_function("cos")$func(x))

  expect_equal(expose_function("tan")$func(x), tan(x))

  args <- c(3.3, 5.5, 6.6)
  expect_equal(expose_function("sum", args)$func(0), sum(args))
  expect_equal(expose_function("prod", args)$func(1), prod(args))
  expect_equal(expose_function("min", args)$func(2), 2)
  expect_equal(expose_function("max", args)$func(2), 6.6)

  expect_error(expose_function("pow", discr = TRUE)$func(x),
               regexp = "undefined method")

  expect_error(expose_function("max", args)$prim(2),
               regexp = "undefined method")
})

test_that("Location parameter works as expected", {
  x <- seq(-10, 10, 0.1)
  sin_func <- expose_function("sin", loc = -0.5 * pi)$func(x)
  cos_func <- expose_function("cos")$func(x)

  expect_equal(sin_func, cos_func)
})

test_that("Arithmetic functions (minus, sub, div) are correct", {
  x <- seq(-5, 5, 0.5)
  expect_equal(expose_function("minus")$func(x), -x)
  expect_equal(expose_function("sub", x2 = 3)$func(x), x - 3)
  expect_equal(expose_function("div", x2 = 4)$func(x), x / 4)
  # identity-like checks
  expect_equal(expose_function("sub", x2 = 0)$func(x), x)
  expect_equal(expose_function("div", x2 = 1)$func(x), x)
})

test_that("Rounding functions (round, ceil, floor) are correct", {
  # NOTE: avoid half-values: C++ std::round uses round-half-away-from-zero
  # while R's round() uses round-half-to-even (banker's rounding)
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  expect_equal(expose_function("round")$func(x), round(x))
  expect_equal(expose_function("ceil")$func(x), ceiling(x))
  expect_equal(expose_function("floor")$func(x), floor(x))
  # integers are unchanged by floor/ceil/round
  x_int <- c(-3, -1, 0, 2, 4)
  expect_equal(expose_function("round")$func(x_int), x_int)
  expect_equal(expose_function("ceil")$func(x_int), x_int)
  expect_equal(expose_function("floor")$func(x_int), x_int)
})

test_that("cos primitive equals sin", {
  x <- seq(-pi, pi, 0.1)
  expect_equal(expose_function("cos")$prim(x), sin(x))
})

test_that("Location parameter with exp shifts argument correctly", {
  x <- seq(-3, 3, 0.5)
  loc_val <- 1.5
  expect_equal(expose_function("exp", loc = loc_val)$func(x), exp(x - loc_val))
  expect_equal(expose_function("exp", loc = loc_val)$prim(x), exp(x - loc_val))
})

test_that("Numerical recipes and distributions are correct", {
  x <- seq(1, 10, 1)
  args <- list(size = 12, prob = 0.2)
  distr <- expose_distribution("binomial", args, min = 1)

  # truncation
  shift <- stats::pbinom(0, size = args$size, prob = args$prob)

  # instantaneous
  y <- stats::dbinom(x, size = args$size, prob = args$prob)
  expect_equal(distr$pdf(x), y / (1. - shift))

  # cumulative
  y <- stats::pbinom(x, size = args$size, prob = args$prob)
  expect_equal(distr$cdf(x), (y - shift) / (1. - shift))
})

test_that("Normal distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("normal", mean = 2, sd = 3)
  expect_equal(distr$pdf(x), stats::dnorm(x, mean = 2, sd = 3))
  expect_equal(distr$cdf(x), stats::pnorm(x, mean = 2, sd = 3))
  expect_equal(expose_distribution("normal")$pdf(x), stats::dnorm(x))

  # test truncation
  low <- 0
  up <- 5
  distr <- expose_distribution("normal", mean = 2, sd = 3, min = low, max = up)
  normpdf <- function(x) stats::dnorm(x, mean = 2, sd = 3)
  normcdf <- function(x) stats::pnorm(x, mean = 2, sd = 3)
  zz <- normcdf(5) - normcdf(low)

  expec_pdf <- ifelse(x >= low & x <= up, normpdf(x), 0) / zz
  expect_equal(distr$pdf(x), expec_pdf)

  expec_cdf <- ifelse(
    x >= 0, ifelse(x <= up, (normcdf(x) - normcdf(low)) / zz, 1), 0
  )
  expect_equal(distr$cdf(x), expec_cdf)
})

test_that("Uniform distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("uniform", min = -1, max = 5)
  expect_equal(distr$pdf(x), stats::dunif(x, min = -1, max = 5))
  expect_equal(distr$cdf(x), stats::punif(x, min = -1, max = 5))
  expect_equal(expose_distribution("uniform")$pdf(x), stats::dunif(x))
})

test_that("Lognormal distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(distr$pdf(x), stats::dlnorm(x, meanlog = 1, sdlog = 0.5))
  expect_equal(distr$cdf(x), stats::plnorm(x, meanlog = 1, sdlog = 0.5))
  expect_equal(expose_distribution("lognormal")$pdf(x), stats::dlnorm(x))
})

test_that("Beta distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("beta", shape1 = 2, shape2 = 1.5)
  expect_equal(distr$pdf(x), stats::dbeta(x, shape1 = 2, shape2 = 1.5))
  expect_equal(distr$cdf(x), stats::pbeta(x, shape1 = 2, shape2 = 1.5))
})

test_that("Gamma distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("gamma", shape = 2, scale = 1.2)
  expect_equal(distr$pdf(x), stats::dgamma(x, shape = 2, scale = 1.2))
  expect_equal(distr$cdf(x), stats::pgamma(x, shape = 2, scale = 1.2))
})

test_that("Gamma (mean/var) distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("gamma_meanvar", mean = 2, var = 1.1)
  shape <- (2^2) / 1.1
  scale <- 1.1 / 2
  expect_equal(distr$pdf(x), stats::dgamma(x, shape = shape, scale = scale))
  expect_equal(distr$cdf(x), stats::pgamma(x, shape = shape, scale = scale))
})

test_that("Exponential distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("exponential", rate = 0.25)
  expect_equal(distr$pdf(x), stats::dexp(x, rate = 0.25))
  expect_equal(distr$cdf(x), stats::pexp(x, rate = 0.25))
  expect_equal(expose_distribution("exponential")$pdf(x), stats::dexp(x))
})

test_that("Poisson distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("poisson", lambda = 5)
  expect_equal(distr$pdf(x), stats::dpois(x, lambda = 5))
  expect_equal(distr$cdf(x), stats::ppois(x, lambda = 5))
})

test_that("Binomial distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("binomial", size = 10, prob = 0.5)
  expect_equal(distr$pdf(x), stats::dbinom(x, size = 10, prob = 0.5))
  expect_equal(distr$cdf(x), stats::pbinom(x, size = 10, prob = 0.5))
})

test_that("Negative binomial distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("negative_binomial", size = 10, prob = 0.5)
  expect_equal(distr$pdf(x), stats::dnbinom(x, size = 10, prob = 0.5))
  expect_equal(distr$cdf(x), stats::pnbinom(x, size = 10, prob = 0.5))
})

test_that("Geometric distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("geometric", prob = 0.5)
  expect_equal(distr$pdf(x), stats::dgeom(x, prob = 0.5))
  expect_equal(distr$cdf(x), stats::pgeom(x, prob = 0.5))
})

test_that("Weibull distribution is correct", {
  x <- seq(-100, 100, 1)
  distr <- expose_distribution("weibull", shape = 2, scale = 1)
  expect_equal(distr$pdf(x), stats::dweibull(x, shape = 2, scale = 1))
  expect_equal(distr$cdf(x), stats::pweibull(x, shape = 2, scale = 1))
})

test_that("Rayleigh distribution is correct", {
  x <- seq(0, 10, 0.5)
  scale <- 2
  rayleigh_pdf <- function(x, s) {
    ifelse(x >= 0, x / s^2 * exp(-x^2 / (2 * s^2)), 0)
  }
  rayleigh_cdf <- function(x, s) {
    ifelse(x >= 0, 1 - exp(-x^2 / (2 * s^2)), 0)
  }

  distr <- expose_distribution("rayleigh", scale = scale)
  expect_equal(distr$pdf(x), rayleigh_pdf(x, scale))
  expect_equal(distr$cdf(x), rayleigh_cdf(x, scale))
})

test_that("CDF is non-decreasing and bounded in [0, 1]", {
  x <- seq(-10, 20, 0.1)
  distrs <- list(
    list(name = "normal", mean = 3, sd = 2),
    list(name = "gamma", shape = 2, scale = 1.5),
    list(name = "exponential", rate = 0.4),
    list(name = "lognormal", meanlog = 1, sdlog = 0.5),
    list(name = "weibull", shape = 2, scale = 2),
    list(name = "rayleigh", scale = 2)
  )
  for (d in distrs) {
    distr <- do.call(expose_distribution, d)
    cdf_vals <- distr$cdf(x)
    pdf_vals <- distr$pdf(x)
    expect_true(all(diff(cdf_vals) >= -1e-12),
                label = paste(d$name, "CDF is non-decreasing"))
    expect_true(all(cdf_vals >= 0 & cdf_vals <= 1),
                label = paste(d$name, "CDF in [0, 1]"))
    expect_true(all(pdf_vals >= 0),
                label = paste(d$name, "PDF is non-negative"))
  }
})

test_that("CDF is the integral of PDF for continuous distributions", {
  # gamma: integrate pdf from 0 to 15, expect cdf(15)
  x <- seq(0.001, 15, length.out = 3000)
  distr <- expose_distribution("gamma", shape = 2, scale = 1.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
               tolerance = 1e-3)

  # lognormal: integrate pdf from 0 to 20, expect cdf(20)
  x <- seq(0.001, 20, length.out = 3000)
  distr <- expose_distribution("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
               tolerance = 1e-3)

  # rayleigh: integrate pdf from 0 to 12, expect cdf(12)
  x <- seq(0.001, 12, length.out = 3000)
  distr <- expose_distribution("rayleigh", scale = 2)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
               tolerance = 1e-3)
})

test_that("CDF equals cumsum of PDF for discrete distributions", {
  x <- seq(0, 20, 1)

  # poisson: CDF(k) = sum of PMF from 0 to k
  distr <- expose_distribution("poisson", lambda = 5)
  expect_equal(distr$cdf(x), cumsum(distr$pdf(x)), tolerance = 1e-10)

  # geometric: CDF(k) = sum of PMF from 0 to k
  distr <- expose_distribution("geometric", prob = 0.3)
  expect_equal(distr$cdf(x), cumsum(distr$pdf(x)), tolerance = 1e-10)
})

test_that("Non-integer x maps to floor for Discrete::Defined distributions", {
  lambda <- 3
  distr <- expose_distribution("poisson", lambda = lambda)
  x_floor <- seq(0, 9, 1)
  for (k in x_floor) {
    # pdf at k and any value in [k, k+1) should be identical
    expect_equal(distr$pdf(k), distr$pdf(k + 0.7))
    expect_equal(distr$cdf(k), distr$cdf(k + 0.7))
  }
})

test_that("Lognormal truncation is correct", {
  min_val <- 1
  max_val <- 5
  x <- seq(0, 8, 0.5)
  distr <- expose_distribution("lognormal", meanlog = 1, sdlog = 0.5,
                               min = min_val, max = max_val)

  lnpdf <- function(x) stats::dlnorm(x, meanlog = 1, sdlog = 0.5)
  lncdf <- function(x) stats::plnorm(x, meanlog = 1, sdlog = 0.5)
  zz <- lncdf(max_val) - lncdf(min_val)

  expec_pdf <- ifelse(x >= min_val & x <= max_val, lnpdf(x), 0) / zz
  expec_cdf <- ifelse(
    x < min_val, 0,
    ifelse(x > max_val, 1, (lncdf(x) - lncdf(min_val)) / zz)
  )

  expect_equal(distr$pdf(x), expec_pdf)
  expect_equal(distr$cdf(x), expec_cdf)
})

test_that("Exponential truncation at max is correct", {
  rate <- 0.3
  max_val <- 5
  x <- seq(0, 8, 0.5)
  distr <- expose_distribution("exponential", rate = rate, max = max_val)

  exppdf <- function(x) stats::dexp(x, rate = rate)
  expcdf <- function(x) stats::pexp(x, rate = rate)
  zz <- expcdf(max_val)

  expec_pdf <- ifelse(x >= 0 & x <= max_val, exppdf(x), 0) / zz
  expec_cdf <- ifelse(x < 0, 0, ifelse(x > max_val, 1, expcdf(x) / zz))

  expect_equal(distr$pdf(x), expec_pdf)
  expect_equal(distr$cdf(x), expec_cdf)
})

test_that("Discretized exponential uses left-discretization", {
  rate <- 0.4
  x <- seq(0, 5, 1)
  distr <- expose_distribution("exponential", rate = rate, discr = TRUE)

  # left discretization: pdf(k) = pexp(k+1) - pexp(k)
  expcdf <- function(x) stats::pexp(x, rate = rate)
  expected_pmf <- expcdf(x + 1) - expcdf(x)
  expect_equal(distr$pdf(x), expected_pmf)
  expect_equal(distr$cdf(x), cumsum(expected_pmf), tolerance = 1e-10)
})

test_that("Discretized gamma uses left-discretization", {
  shape <- 2
  scale <- 1.5
  x <- seq(0, 10, 1)
  distr <- expose_distribution("gamma", shape = shape,
                               scale = scale, discr = TRUE)

  gammcdf <- function(x) stats::pgamma(x, shape = shape, scale = scale)
  expect_equal(distr$pdf(x), gammcdf(x + 1) - gammcdf(x))
})

test_that("Discretized normal uses center-discretization", {
  mu <- 5
  sd <- 2
  x <- seq(0, 10, 1)
  distr <- expose_distribution("normal", mean = mu, sd = sd, discr = TRUE)

  # center discretization: pdf(k) = pnorm(k+0.5) - pnorm(k-0.5)
  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  expect_equal(distr$pdf(x), normcdf(x + 0.5) - normcdf(x - 0.5))
})

test_that("Ferretti20 infectivity profile behave has expected", {
  x <- seq(0, 15, 0.001)
  profile <- expose_infec_profile("ferretti20", incub_period = 5)
  expect_equal(trapezoid(x, profile$inst(x)), profile$cumul(x)[length(x)])

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  profile <- expose_infec_profile("ferretti20", incub_period = 5, discr = TRUE)
  expect_true(all(profile$inst(x) == profile$cumul(1)))
})
