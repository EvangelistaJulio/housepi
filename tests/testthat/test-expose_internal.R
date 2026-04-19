library(testthat)

trapezoid <- function(x, y) {
  0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))
}

########################################
## EXPOSE INTERNAL : MEMORY STABILITY ##
########################################

test_that("Expose function has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  fun <- expose_function("sin")
  rm(fun)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in seq_len(200)) {
    fun <- expose_function("sin")
    rm(fun)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})

#################################
## EXPOSE INTERNAL : FUNCTIONS ##
#################################

test_that("Power function is correct", {
  x <- seq(0.1, 10, 0.1)
  x2 <- 3
  fn <- expose_function("pow", x2 = x2)
  expect_equal(fn$args(), c(x2 = x2))
  expect_equal(fn$func(x), x^x2)
  expect_equal(fn$prim(x), x^(x2 + 1) / (x2 + 1))
  expect_equal(fn$log_func(x), x2 * log(x))
  expect_equal(fn$log_ratio(x, x + 0.1),
               log(fn$func(x + 0.1) / fn$func(x)))

  # x2 = -1 primitive uses log(|x|)
  expect_equal(expose_function("pow", x2 = -1)$prim(x), log(abs(x)))
  # x2 = 0.5 equals sqrt
  expect_equal(expose_function("pow", x2 = 0.5)$func(x), sqrt(x))

  # error cases
  expect_error(expose_function("pow")$func(x), regexp = "Missing non default")
  expect_error(expose_function("pow", 5), regexp = "Unnamed argument")
  expect_error(expose_function("pow", discr = TRUE)$func(x),
               regexp = "undefined method")
})

test_that("Square root function is correct", {
  x <- seq(0.1, 10, 0.1)
  fn <- expose_function("sqrt")
  expect_equal(fn$func(x), sqrt(x))
  expect_equal(fn$prim(x), x^1.5 / 1.5)
  expect_equal(fn$log_func(x), 0.5 * log(x))
  expect_equal(fn$log_ratio(x, x + 0.1),
               log(fn$func(x + 0.1) / fn$func(x)))
})

test_that("Exponential function is correct", {
  x <- seq(-5, 5, 0.1)
  fn <- expose_function("exp")
  expect_equal(fn$func(x), exp(x))
  expect_equal(fn$prim(x), exp(x))
  expect_equal(fn$log_func(x), x)
  expect_equal(fn$log_ratio(x, x + pi),
               log(fn$func(x + pi) / fn$func(x)))

  # location parameter
  loc_val <- 1.5
  expect_equal(expose_function("exp", loc = loc_val)$func(x), exp(x - loc_val))
  expect_equal(expose_function("exp", loc = loc_val)$prim(x), exp(x - loc_val))
})

test_that("Logarithm function is correct", {
  x <- seq(1.1, 20, 0.1)
  fn <- expose_function("log")
  expect_equal(fn$func(x), log(x))
  expect_equal(fn$prim(x), x * log(x) - x)
  expect_equal(fn$log_func(x), log(log(x)))
  expect_equal(fn$log_ratio(x, x + pi),
               log(fn$func(x + pi) / fn$func(x)))
})

test_that("Sine function is correct", {
  x <- seq(0.1, pi - 0.1, 0.05)
  fn <- expose_function("sin")
  expect_equal(fn$func(x), sin(x))
  expect_equal(fn$prim(x), -cos(x))
  expect_equal(fn$log_func(x), log(sin(x)))
  expect_equal(fn$log_ratio(x, x + 0.1),
               log(fn$func(x + 0.1) / fn$func(x)))

  # location parameter: sin(x + pi/2) = cos(x)
  x_full <- seq(-10, 10, 0.1)
  expect_equal(expose_function("sin", loc = -0.5 * pi)$func(x_full),
               cos(x_full))
})

test_that("Cosine function is correct", {
  x <- seq(-pi / 2 + 0.1, pi / 2 - 0.1, 0.05)
  fn <- expose_function("cos")
  expect_equal(fn$func(x), cos(x))
  expect_equal(fn$prim(x), sin(x))
  expect_equal(fn$log_func(x), log(cos(x)))
  expect_equal(fn$log_ratio(x, x + 0.1),
               log(fn$func(x + 0.1) / fn$func(x)))
})

test_that("Tangent function is correct", {
  x <- seq(0.1, pi / 2 - 0.1, 0.05)
  fn <- expose_function("tan")
  expect_equal(fn$func(x), tan(x))
  expect_equal(fn$log_func(x), log(tan(x)))
  expect_equal(fn$log_ratio(x, x + 0.1),
               log(fn$func(x + 0.1) / fn$func(x)))
  # primitive is undefined
  expect_error(fn$prim(x), regexp = "undefined method")
})

test_that("Minus function is correct", {
  x <- seq(-5, 5, 0.5)
  fn <- expose_function("minus")
  expect_equal(fn$func(x), -x)
  # default raw_log_ratio returns 0
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  expect_error(fn$prim(x), regexp = "undefined method")
  expect_error(fn$log_func(x), regexp = "undefined method")
})

test_that("Subtraction function is correct", {
  x <- seq(-5, 5, 0.5)
  fn <- expose_function("sub", x2 = 3)
  expect_equal(fn$func(x), x - 3)
  expect_equal(expose_function("sub", x2 = 0)$func(x), x)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  expect_error(fn$prim(x), regexp = "undefined method")
  expect_error(fn$log_func(x), regexp = "undefined method")
})

test_that("Division function is correct", {
  x <- seq(-5, 5, 0.5)
  fn <- expose_function("div", x2 = 4)
  expect_equal(fn$func(x), x / 4)
  expect_equal(expose_function("div", x2 = 1)$func(x), x)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  expect_error(fn$prim(x), regexp = "undefined method")
  expect_error(fn$log_func(x), regexp = "undefined method")
})

test_that("Summation function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- expose_function("sum", args)
  expect_equal(fn$func(0), sum(args))
  expect_equal(fn$log_ratio(0, pi), 0)
  expect_error(fn$prim(0), regexp = "undefined method")
  expect_error(fn$log_func(0), regexp = "undefined method")
})

test_that("Product function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- expose_function("prod", args)
  expect_equal(fn$func(1), prod(args))
  expect_equal(fn$log_ratio(1, pi), 0)
  expect_error(fn$prim(1), regexp = "undefined method")
  expect_error(fn$log_func(1), regexp = "undefined method")
})

test_that("Minimum function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- expose_function("min", args)
  expect_equal(fn$func(2), 2)
  expect_equal(fn$log_ratio(2, pi), 0)
  expect_error(fn$prim(2), regexp = "undefined method")
  expect_error(fn$log_func(1), regexp = "undefined")
})

test_that("Maximum function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- expose_function("max", args)
  expect_equal(fn$func(2), 6.6)
  expect_equal(fn$log_ratio(2, pi), 0)
  expect_error(fn$prim(2), regexp = "undefined method")
  expect_error(fn$log_func(1), regexp = "undefined")
})

test_that("Round function is correct", {
  # NOTE: avoid half-values: C++ std::round uses round-half-away-from-zero
  # while R's round() uses round-half-to-even (crazy...)
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  x_int <- c(-3, -1, 0, 2, 4)
  fn <- expose_function("round")
  expect_equal(fn$func(x), round(x))
  expect_equal(fn$func(x_int), x_int)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  expect_error(fn$prim(x), regexp = "undefined method")
  expect_error(fn$log_func(x), regexp = "undefined method")
})

test_that("Ceil function is correct", {
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  x_int <- c(-3, -1, 0, 2, 4)
  fn <- expose_function("ceil")
  expect_equal(fn$func(x), ceiling(x))
  expect_equal(fn$func(x_int), x_int)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  expect_error(fn$prim(x), regexp = "undefined method")
  expect_error(fn$log_func(x), regexp = "undefined method")
})

test_that("Floor function is correct", {
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  x_int <- c(-3, -1, 0, 2, 4)
  fn <- expose_function("floor")
  expect_equal(fn$func(x), floor(x))
  expect_equal(fn$func(x_int), x_int)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  expect_error(fn$prim(x), regexp = "undefined method")
  expect_error(fn$log_func(x), regexp = "undefined method")
})

#####################################
## EXPOSE INTERNAL : DISTRIBUTIONS ##
#####################################

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
  x <- seq(-20, 20, 1)
  expect_equal(expose_distribution("normal")$pdf(x), stats::dnorm(x))
  distr <- expose_distribution("normal", mean = 2, sd = 3)
  expect_equal(distr$pdf(x), stats::dnorm(x, mean = 2, sd = 3))
  expect_equal(distr$cdf(x), stats::pnorm(x, mean = 2, sd = 3))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))

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
  x <- seq(-20, 20, 1)
  expect_equal(expose_distribution("uniform")$pdf(x), stats::dunif(x))
  distr <- expose_distribution("uniform", min = -1, max = 5)
  expect_equal(distr$pdf(x), stats::dunif(x, min = -1, max = 5))
  expect_equal(distr$cdf(x), stats::punif(x, min = -1, max = 5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Lognormal distribution is correct", {
  x <- seq(-20, 20, 1)
  expect_equal(expose_distribution("lognormal")$pdf(x), stats::dlnorm(x))
  distr <- expose_distribution("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(distr$pdf(x), stats::dlnorm(x, meanlog = 1, sdlog = 0.5))
  expect_equal(distr$cdf(x), stats::plnorm(x, meanlog = 1, sdlog = 0.5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Beta distribution is correct", {
  x <- seq(-10, 10, 0.1)
  distr <- expose_distribution("beta", shape1 = 2, shape2 = 1.5)
  expect_equal(distr$pdf(x), stats::dbeta(x, shape1 = 2, shape2 = 1.5))
  expect_equal(distr$cdf(x), stats::pbeta(x, shape1 = 2, shape2 = 1.5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + 0.1),
               log(distr$pdf(x + 0.1) / distr$pdf(x)))
})

test_that("Gamma distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("gamma", shape = 2, scale = 1.2)
  expect_equal(distr$pdf(x), stats::dgamma(x, shape = 2, scale = 1.2))
  expect_equal(distr$cdf(x), stats::pgamma(x, shape = 2, scale = 1.2))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Gamma (mean/var) distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("gamma_meanvar", mean = 2, var = 1.1)
  shape <- (2^2) / 1.1
  scale <- 1.1 / 2
  expect_equal(distr$pdf(x), stats::dgamma(x, shape = shape, scale = scale))
  expect_equal(distr$cdf(x), stats::pgamma(x, shape = shape, scale = scale))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Exponential distribution is correct", {
  x <- seq(-20, 20, 1)
  expect_equal(expose_distribution("exponential")$pdf(x), stats::dexp(x))
  distr <- expose_distribution("exponential", rate = 0.25)
  expect_equal(distr$pdf(x), stats::dexp(x, rate = 0.25))
  expect_equal(distr$cdf(x), stats::pexp(x, rate = 0.25))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Poisson distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("poisson", lambda = 5)
  expect_equal(distr$pdf(x), stats::dpois(x, lambda = 5))
  expect_equal(distr$cdf(x), stats::ppois(x, lambda = 5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Binomial distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("binomial", size = 10, prob = 0.5)
  expect_equal(distr$pdf(x), stats::dbinom(x, size = 10, prob = 0.5))
  expect_equal(distr$cdf(x), stats::pbinom(x, size = 10, prob = 0.5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Negative binomial distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("negative_binomial", size = 10, prob = 0.5)
  expect_equal(distr$pdf(x), stats::dnbinom(x, size = 10, prob = 0.5))
  expect_equal(distr$cdf(x), stats::pnbinom(x, size = 10, prob = 0.5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Geometric distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("geometric", prob = 0.5)
  expect_equal(distr$pdf(x), stats::dgeom(x, prob = 0.5))
  expect_equal(distr$cdf(x), stats::pgeom(x, prob = 0.5))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("Weibull distribution is correct", {
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("weibull", shape = 2, scale = 1)
  expect_equal(distr$pdf(x), stats::dweibull(x, shape = 2, scale = 1))
  expect_equal(distr$cdf(x), stats::pweibull(x, shape = 2, scale = 1))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
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
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(distr$log_ratio(x, x + pi),
               log(distr$pdf(x + pi) / distr$pdf(x)))
})

test_that("CDF and quantile make sense", {
  x <- seq(-10, 20, 0.1)
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  distrs <- list(list(name = "uniform", min = -3, max = 2),
                 list(name = "normal", mean = 3, sd = 2),
                 list(name = "gamma", shape = 2, scale = 1.5),
                 list(name = "gamma_meanvar", mean = 2.2, var = 1.1),
                 list(name = "exponential", rate = 0.4),
                 list(name = "lognormal", meanlog = 1, sdlog = 0.5),
                 list(name = "weibull", shape = 2, scale = 2),
                 list(name = "rayleigh", scale = 2))
  for (d in distrs) {
    distr <- do.call(expose_distribution, d)
    cdf_vals <- distr$cdf(x)
    pdf_vals <- distr$pdf(x)
    expect_true(all(diff(cdf_vals) >= -1e-12))
    expect_true(all(cdf_vals >= 0 & cdf_vals <= 1))
    expect_true(all(pdf_vals >= 0))
    expect_equal(distr$cdf(distr$quantile(p)), p)
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

test_that("quantile matches R reference for standard distributions", {
  p <- seq(0, 1, 0.1)

  distr <- expose_distribution("normal", mean = 2, sd = 3)
  expect_equal(distr$quantile(p), stats::qnorm(p, mean = 2, sd = 3))

  distr <- expose_distribution("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(distr$quantile(p), stats::qlnorm(p, meanlog = 1, sdlog = 0.5))

  distr <- expose_distribution("gamma", shape = 2, scale = 1.5)
  expect_equal(distr$quantile(p), stats::qgamma(p, shape = 2, scale = 1.5))

  distr <- expose_distribution("exponential", rate = 0.3)
  expect_equal(distr$quantile(p), stats::qexp(p, rate = 0.3))

  distr <- expose_distribution("beta", shape1 = 2, shape2 = 3)
  expect_equal(distr$quantile(p), stats::qbeta(p, shape1 = 2, shape2 = 3))
})

test_that("quantile matches R for discrete distributions (p >= 0.5)", {
  p <- seq(0, 1, 0.1)

  distr <- expose_distribution("poisson", lambda = 5)
  expect_equal(distr$quantile(p), stats::qpois(p, lambda = 5))

  distr <- expose_distribution("binomial", size = 10, prob = 0.4)
  expect_equal(distr$quantile(p), stats::qbinom(p, size = 10, prob = 0.4))

  distr <- expose_distribution("geometric", prob = 0.3)
  expect_equal(distr$quantile(p), stats::qgeom(p, prob = 0.3))

  distr <- expose_distribution("negative_binomial", size = 5, prob = 0.4)
  expect_equal(distr$quantile(p), stats::qnbinom(p, size = 5, prob = 0.4))
})

############################################
## EXPOSE INTERNAL : INFECTIVITY PROFILES ##
############################################

test_that("Ferretti20 infectivity profile behave has expected", {
  x <- seq(0, 15, 0.001)
  profile <- expose_infec_profile("ferretti20", incub_period = 5)
  expect_equal(trapezoid(x, profile$inst(x)), profile$cumul(x)[length(x)])

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  profile <- expose_infec_profile("ferretti20", incub_period = 5, discr = TRUE)
  expect_true(all(profile$inst(x) == profile$cumul(1)))
})
