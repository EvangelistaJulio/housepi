library(testthat)

trapezoid <- function(x, y) {
  0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))
}

#######################################
## BROWSE INTERNAL : AVAILABLE NAMES ##
#######################################

# NOTE: the custom registry is populated at DLL load, so `growth` joins the
# function listing once test-model_generator.R compiles inst/examples/custom.
# This file sorts first, hence the snapshot sees a pristine registry.
test_that("Available names listing is stable", {
  expect_snapshot(browse_functions())
  expect_snapshot(browse_distributions())
  expect_snapshot(browse_infectivity_profiles())
})

########################################
## BROWSE INTERNAL : MEMORY STABILITY ##
########################################

test_that("Expose function has no leaks", {
  mem_used_mb <- function() gc()["Vcells", "(Mb)"]

  fun <- browse_functions("sin")
  rm(fun)
  invisible(gc())
  mem_before <- mem_used_mb()

  # heavy allocation loop
  for (i in seq_len(200)) {
    fun <- browse_functions("sin")
    rm(fun)
    if (i %% 20 == 0) gc()
  }
  invisible(gc())
  expect_lt(mem_used_mb() - mem_before, 5)
})

#################################
## BROWSE INTERNAL : FUNCTIONS ##
#################################

test_that("Power function is correct", {
  x <- seq(0.1, 10, 0.1)
  x2 <- 3
  fn <- browse_functions("pow", x2 = x2)
  expect_equal(fn$args(), c(x2 = x2))
  expect_equal(fn$func(x), x^x2)
  expect_equal(fn$prim(x), x^(x2 + 1) / (x2 + 1))
  expect_equal(fn$log_func(x), x2 * log(x))
  expect_equal(fn$log_ratio(x, x + 0.1),
               log(fn$func(x + 0.1) / fn$func(x)))
  fn$log_ratio(x, x[-1]) |>
    expect_error("[Ss]ize")

  # x2 = -1 primitive uses log(|x|)
  expect_equal(browse_functions("pow", x2 = -1)$prim(x), log(abs(x)))
  # x2 = 0.5 equals sqrt
  expect_equal(browse_functions("pow", x2 = 0.5)$func(x), sqrt(x))

  # error cases
  browse_functions("pow")$func(x) |>
    expect_error("Missing non default")
  browse_functions("pow", 5) |>
    expect_error("Unnamed argument")
  browse_functions("pow", x2 = 5, y2 = 3) |>
    expect_error("Unused argument.*y2.*[Dd]id you mean")
  browse_functions("pow", discr = TRUE)$func(x) |>
    expect_error("undefined method")
})

test_that("Square root function is correct", {
  x <- seq(0.1, 10, 0.1)
  fn <- browse_functions("sqrt")
  expect_equal(fn$func(x), sqrt(x))
  expect_equal(fn$prim(x), x^1.5 / 1.5)
  expect_equal(fn$log_func(x), 0.5 * log(x))
  expect_equal(
    fn$log_ratio(x, x + 0.1),
    log(fn$func(x + 0.1) / fn$func(x))
  )
})

test_that("Exponential function is correct", {
  x <- seq(-5, 5, 0.1)
  fn <- browse_functions("exp")
  expect_equal(fn$func(x), exp(x))
  expect_equal(fn$prim(x), exp(x))
  expect_equal(fn$log_func(x), x)
  expect_equal(
    fn$log_ratio(x, x + pi),
    log(fn$func(x + pi) / fn$func(x))
  )

  # location parameter
  loc_val <- 1.5
  expect_equal(browse_functions("exp", loc = loc_val)$func(x), exp(x - loc_val))
  expect_equal(browse_functions("exp", loc = loc_val)$prim(x), exp(x - loc_val))
})

test_that("Logarithm function is correct", {
  x <- seq(1.1, 20, 0.1)
  fn <- browse_functions("log")
  expect_equal(fn$func(x), log(x))
  expect_equal(fn$prim(x), x * log(x) - x)
  expect_equal(fn$log_func(x), log(log(x)))
  expect_equal(
    fn$log_ratio(x, x + pi),
    log(fn$func(x + pi) / fn$func(x))
  )
})

test_that("Sine function is correct", {
  x <- seq(0.1, pi - 0.1, 0.05)
  fn <- browse_functions("sin")
  expect_equal(fn$func(x), sin(x))
  expect_equal(fn$prim(x), -cos(x))
  expect_equal(fn$log_func(x), log(sin(x)))
  expect_equal(
    fn$log_ratio(x, x + 0.1),
    log(fn$func(x + 0.1) / fn$func(x))
  )

  # location parameter: sin(x + pi/2) = cos(x)
  x_full <- seq(-10, 10, 0.1)
  expect_equal(
    browse_functions("sin", loc = -0.5 * pi)$func(x_full),
    cos(x_full)
  )
})

test_that("Cosine function is correct", {
  x <- seq(-pi / 2 + 0.1, pi / 2 - 0.1, 0.05)
  fn <- browse_functions("cos")
  expect_equal(fn$func(x), cos(x))
  expect_equal(fn$prim(x), sin(x))
  expect_equal(fn$log_func(x), log(cos(x)))
  expect_equal(
    fn$log_ratio(x, x + 0.1),
    log(fn$func(x + 0.1) / fn$func(x))
  )
})

test_that("Tangent function is correct", {
  x <- seq(0.1, pi / 2 - 0.1, 0.05)
  fn <- browse_functions("tan")
  expect_equal(fn$func(x), tan(x))
  expect_equal(fn$log_func(x), log(tan(x)))
  expect_equal(
    fn$log_ratio(x, x + 0.1),
    log(fn$func(x + 0.1) / fn$func(x))
  )
  # primitive is undefined
  fn$prim(x) |>
    expect_error("undefined method")
})

test_that("Minus function is correct", {
  x <- seq(-5, 5, 0.5)
  fn <- browse_functions("minus")
  expect_equal(fn$func(x), -x)
  # default raw_log_ratio returns 0
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  fn$prim(x) |>
    expect_error("undefined method")
  fn$log_func(x) |>
    expect_error("undefined method")
})

test_that("Subtraction function is correct", {
  x <- seq(-5, 5, 0.5)
  fn <- browse_functions("sub", x2 = 3)
  expect_equal(fn$func(x), x - 3)
  expect_equal(browse_functions("sub", x2 = 0)$func(x), x)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  fn$prim(x) |>
    expect_error("undefined method")
  fn$log_func(x) |>
    expect_error("undefined method")
})

test_that("Division function is correct", {
  x <- seq(-5, 5, 0.5)
  fn <- browse_functions("div", x2 = 4)
  expect_equal(fn$func(x), x / 4)
  expect_equal(browse_functions("div", x2 = 1)$func(x), x)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  fn$prim(x) |>
    expect_error("undefined method")
  fn$log_func(x) |>
    expect_error("undefined method")
})

test_that("Summation function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- browse_functions("sum", args)
  expect_equal(fn$func(0), sum(args))
  expect_equal(fn$log_ratio(0, pi), 0)
  fn$prim(0) |>
    expect_error("undefined method")
  fn$log_func(0) |>
    expect_error("undefined method")
})

test_that("Product function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- browse_functions("prod", args)
  expect_equal(fn$func(1), prod(args))
  expect_equal(fn$log_ratio(1, pi), 0)
  fn$prim(1) |>
    expect_error("undefined method")
  fn$log_func(1) |>
    expect_error("undefined method")
})

test_that("Minimum function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- browse_functions("min", args)
  expect_equal(fn$func(2), 2)
  expect_equal(fn$log_ratio(2, pi), 0)
  fn$prim(2) |>
    expect_error("undefined method")
  fn$log_func(1) |>
    expect_error("undefined")
})

test_that("Maximum function is correct", {
  args <- c(3.3, 5.5, 6.6)
  fn <- browse_functions("max", args)
  expect_equal(fn$func(2), 6.6)
  expect_equal(fn$log_ratio(2, pi), 0)
  fn$prim(2) |>
    expect_error("undefined method")
  fn$log_func(1) |>
    expect_error("undefined")
})

test_that("Round function is correct", {
  # NOTE: avoid half-values: C++ std::round uses round-half-away-from-zero
  # while R's round() uses round-half-to-even (crazy...)
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  x_int <- c(-3, -1, 0, 2, 4)
  fn <- browse_functions("round")
  expect_equal(fn$func(x), round(x))
  expect_equal(fn$func(x_int), x_int)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  fn$prim(x) |>
    expect_error("undefined method")
  fn$log_func(x) |>
    expect_error("undefined method")
})

test_that("Ceil function is correct", {
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  x_int <- c(-3, -1, 0, 2, 4)
  fn <- browse_functions("ceil")
  expect_equal(fn$func(x), ceiling(x))
  expect_equal(fn$func(x_int), x_int)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  fn$prim(x) |>
    expect_error("undefined method")
  fn$log_func(x) |>
    expect_error("undefined method")
})

test_that("Floor function is correct", {
  x <- c(-2.7, -1.3, -0.3, 0.0, 0.4, 1.3, 2.7)
  x_int <- c(-3, -1, 0, 2, 4)
  fn <- browse_functions("floor")
  expect_equal(fn$func(x), floor(x))
  expect_equal(fn$func(x_int), x_int)
  expect_equal(fn$log_ratio(x, x + pi), rep(0, length(x)))
  fn$prim(x) |>
    expect_error("undefined method")
  fn$log_func(x) |>
    expect_error("undefined method")
})

#####################################
## BROWSE INTERNAL : DISTRIBUTIONS ##
#####################################

# Assert pdf, cdf, log_pdf and log_ratio of an exposed distribution against
# user-supplied R reference d/p functions.
expect_distribution_matches <- function(name, args = list(), ref_d, ref_p,
                                        x = seq(-20, 20, 1), h = pi) {
  distr <- do.call(browse_distributions, c(list(name), args))
  expect_equal(distr$pdf(x), ref_d(x))
  expect_equal(distr$cdf(x), ref_p(x))
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
  expect_equal(
    distr$log_ratio(x, x + h),
    log(distr$pdf(x + h) / distr$pdf(x))
  )
}

test_that("Truncated binomial is correct", {
  x <- seq(1, 10, 1)
  size <- 12
  prob <- 0.2
  distr <- browse_distributions("binomial", list(size = size, prob = prob),
    min = 1
  )
  shift <- stats::pbinom(0, size, prob)
  expect_equal(distr$pdf(x), stats::dbinom(x, size, prob) / (1 - shift))
  expect_equal(distr$cdf(x), (stats::pbinom(x, size, prob) - shift) / (1 - shift))
})

test_that("Standard distributions match R reference d/p/log_pdf/log_ratio", {
  rayleigh_d <- function(x, s) ifelse(x >= 0, x / s^2 * exp(-x^2 / (2 * s^2)), 0)
  rayleigh_p <- function(x, s) ifelse(x >= 0, 1 - exp(-x^2 / (2 * s^2)), 0)

  expect_distribution_matches(
    "normal", list(mean = 2, sd = 3),
    function(x) stats::dnorm(x, 2, 3), function(x) stats::pnorm(x, 2, 3)
  )
  expect_distribution_matches(
    "uniform", list(min = -1, max = 5),
    function(x) stats::dunif(x, -1, 5), function(x) stats::punif(x, -1, 5)
  )
  expect_distribution_matches(
    "lognormal", list(meanlog = 1, sdlog = 0.5),
    function(x) stats::dlnorm(x, 1, 0.5), function(x) stats::plnorm(x, 1, 0.5)
  )
  expect_distribution_matches("beta", list(shape1 = 2, shape2 = 1.5),
    function(x) stats::dbeta(x, 2, 1.5), function(x) stats::pbeta(x, 2, 1.5),
    x = seq(-10, 10, 0.1), h = 0.1
  )
  expect_distribution_matches(
    "gamma", list(shape = 2, scale = 1.2),
    function(x) stats::dgamma(x, shape = 2, scale = 1.2),
    function(x) stats::pgamma(x, shape = 2, scale = 1.2)
  )
  expect_distribution_matches(
    "exponential", list(rate = 0.25),
    function(x) stats::dexp(x, 0.25), function(x) stats::pexp(x, 0.25)
  )
  expect_distribution_matches(
    "poisson", list(lambda = 5),
    function(x) stats::dpois(x, 5), function(x) stats::ppois(x, 5)
  )
  expect_distribution_matches(
    "binomial", list(size = 10, prob = 0.5),
    function(x) stats::dbinom(x, 10, 0.5), function(x) stats::pbinom(x, 10, 0.5)
  )
  expect_distribution_matches(
    "negative_binomial",
    list(size = 10, prob = 0.5),
    function(x) stats::dnbinom(x, 10, 0.5), function(x) stats::pnbinom(x, 10, 0.5)
  )
  expect_distribution_matches(
    "geometric", list(prob = 0.5),
    function(x) stats::dgeom(x, 0.5), function(x) stats::pgeom(x, 0.5)
  )
  expect_distribution_matches(
    "weibull", list(shape = 2, scale = 1),
    function(x) stats::dweibull(x, 2, 1), function(x) stats::pweibull(x, 2, 1)
  )
  expect_distribution_matches("rayleigh", list(scale = 2),
    function(x) rayleigh_d(x, 2), function(x) rayleigh_p(x, 2),
    x = seq(0, 10, 0.5)
  )
})

test_that("Categorical distribution is correct", {
  # unnormalized weights still give a normalized pmf
  distr <- browse_distributions("categorical", fortune = 2, fatman = 1, vamp = 7)
  pmf <- c(0.2, 0.1, 0.7)
  expect_equal(distr$pdf(0:2), pmf)
  expect_equal(distr$pdf(c(-1, 3)), c(0, 0))
  expect_equal(distr$log_pdf(0:2), log(pmf))

  # cumulative pmf
  expect_equal(distr$cdf(0:2), cumsum(pmf))
  expect_equal(distr$cdf(c(-1, 5)), c(0, 1))

  # log-ratio
  distr$log_ratio(c(0, 1), c(1, 2)) |>
    expect_equal(log(pmf[c(2, 3)] / pmf[c(1, 2)]))

  browse_distributions("categorical") |>
    expect_error("no weights")
  browse_distributions("categorical", peace = 1, walker = -2) |>
    expect_error("non-strictly positive")
})

test_that("Default parameters match R reference for d", {
  x <- seq(-20, 20, 1)
  expect_equal(browse_distributions("normal")$pdf(x), stats::dnorm(x))
  expect_equal(browse_distributions("uniform")$pdf(x), stats::dunif(x))
  expect_equal(browse_distributions("lognormal")$pdf(x), stats::dlnorm(x))
  expect_equal(browse_distributions("exponential")$pdf(x), stats::dexp(x))
})

test_that("Normal both-sided truncation is correct", {
  low <- 0
  up <- 5
  x <- seq(-20, 20, 1)
  distr <- browse_distributions("normal", mean = 2, sd = 3, min = low, max = up)
  ref_d <- function(x) stats::dnorm(x, 2, 3)
  ref_p <- function(x) stats::pnorm(x, 2, 3)
  zz <- ref_p(up) - ref_p(low)

  expect_equal(distr$pdf(x), ifelse(x >= low & x <= up, ref_d(x), 0) / zz)
  expect_equal(
    distr$cdf(x),
    ifelse(x < low, 0,
      ifelse(x > up, 1, (ref_p(x) - ref_p(low)) / zz)
    )
  )
})

test_that("Gamma accepts shape/scale and mean/var parametrizations", {
  browse_distributions("gamma", mean = 2, scale = 1.1) |>
    expect_error("mix")
  browse_distributions("gamma", shape = 2, var = 1.1) |>
    expect_error("mix")

  m <- 2.5
  v <- 1.7
  shape <- m * m / v
  scale <- v / m
  x <- seq(0.1, 12, 0.3)
  d_can <- browse_distributions("gamma", shape = shape, scale = scale)
  d_alt <- browse_distributions("gamma", mean = m, var = v)
  expect_equal(d_can$pdf(x), d_alt$pdf(x))
  expect_equal(d_can$cdf(x), d_alt$cdf(x))
  expect_equal(d_can$log_pdf(x), d_alt$log_pdf(x))

  # mean/var matches stats::dgamma with derived shape/scale
  expect_equal(d_alt$pdf(x), stats::dgamma(x, shape = shape, scale = scale))
})

test_that("CDF and quantile make sense", {
  x <- seq(-10, 20, 0.1)
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  distrs <- list(
    list(name = "uniform", min = -3, max = 2),
    list(name = "normal", mean = 3, sd = 2),
    list(name = "gamma", shape = 2, scale = 1.5),
    list(name = "gamma", mean = 2.2, var = 1.1),
    list(name = "exponential", rate = 0.4),
    list(name = "lognormal", meanlog = 1, sdlog = 0.5),
    list(name = "weibull", shape = 2, scale = 2),
    list(name = "rayleigh", scale = 2)
  )
  for (d in distrs) {
    distr <- do.call(browse_distributions, d)
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
  distr <- browse_distributions("gamma", shape = 2, scale = 1.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )

  # lognormal: integrate pdf from 0 to 20, expect cdf(20)
  x <- seq(0.001, 20, length.out = 3000)
  distr <- browse_distributions("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )

  # rayleigh: integrate pdf from 0 to 12, expect cdf(12)
  x <- seq(0.001, 12, length.out = 3000)
  distr <- browse_distributions("rayleigh", scale = 2)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )
})

test_that("CDF equals cumsum of PDF for discrete distributions", {
  x <- seq(0, 20, 1)

  # poisson: CDF(k) = sum of PMF from 0 to k
  distr <- browse_distributions("poisson", lambda = 5)
  expect_equal(distr$cdf(x), cumsum(distr$pdf(x)), tolerance = 1e-10)

  # geometric: CDF(k) = sum of PMF from 0 to k
  distr <- browse_distributions("geometric", prob = 0.3)
  expect_equal(distr$cdf(x), cumsum(distr$pdf(x)), tolerance = 1e-10)
})

test_that("Non-integer x maps to floor for Discrete::Defined distributions", {
  lambda <- 3
  distr <- browse_distributions("poisson", lambda = lambda)
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
  distr <- browse_distributions("lognormal",
    meanlog = 1, sdlog = 0.5,
    min = min_val, max = max_val
  )

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
  distr <- browse_distributions("exponential", rate = rate, max = max_val)

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
  distr <- browse_distributions("exponential", rate = rate, discr = TRUE)

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
  distr <- browse_distributions("gamma",
    shape = shape,
    scale = scale, discr = TRUE
  )

  gammcdf <- function(x) stats::pgamma(x, shape = shape, scale = scale)
  expect_equal(distr$pdf(x), gammcdf(x + 1) - gammcdf(x))
})

test_that("Discretized normal uses center-discretization", {
  mu <- 5
  sd <- 2
  x <- seq(0, 10, 1)
  distr <- browse_distributions("normal", mean = mu, sd = sd, discr = TRUE)

  # center discretization: pdf(k) = pnorm(k+0.5) - pnorm(k-0.5)
  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  expect_equal(distr$pdf(x), normcdf(x + 0.5) - normcdf(x - 0.5))
})

test_that("Discretized normal stays finite far in the upper tail", {
  # far above the mean both CDFs round to 1 and the naive difference cancels
  # to 0 (-Inf log-pmf); the complementary-CDF branch keeps it finite.
  x <- seq(9, 30, 1)
  distr <- browse_distributions("normal", mean = 0, sd = 1, discr = TRUE)

  # the bug this guards against: naive CDF difference underflows to exactly 0
  naive <- stats::pnorm(x + 0.5) - stats::pnorm(x - 0.5)
  expect_true(any(naive == 0))

  # stable reference via upper-tail (complementary) CDF
  ccdf <- function(x) stats::pnorm(x, lower.tail = FALSE)
  expected_pmf <- ccdf(x - 0.5) - ccdf(x + 0.5)

  expect_true(all(is.finite(distr$log_pdf(x))))
  expect_equal(distr$pdf(x), expected_pmf)
  expect_equal(distr$log_pdf(x), log(expected_pmf))

  # symmetry: reflecting about the mean matches the lower tail
  expect_equal(distr$pdf(x), distr$pdf(-x))
})

test_that("Discretized lognormal stays finite far in the upper tail", {
  # left-discretized, heavy upper tail: exercises the gamma/erfc complement path
  meanlog <- 0
  sdlog <- 0.5
  x <- seq(60, 200, 5)
  distr <- browse_distributions("lognormal",
    meanlog = meanlog, sdlog = sdlog, discr = TRUE
  )

  lncdf <- function(x) stats::plnorm(x, meanlog = meanlog, sdlog = sdlog)
  expect_true(any(lncdf(x + 1) - lncdf(x) == 0)) # naive underflows

  ccdf <- function(x) {
    stats::plnorm(x, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE)
  }
  expected_pmf <- ccdf(x) - ccdf(x + 1)

  expect_true(all(is.finite(distr$log_pdf(x))))
  expect_equal(distr$pdf(x), expected_pmf)
})

test_that("Discretized gamma stays finite far in the upper tail", {
  # left-discretized; exercises the boost gamma_q complement override
  shape <- 2
  scale <- 1.5
  x <- seq(60, 120, 5)
  distr <- browse_distributions("gamma",
    shape = shape, scale = scale, discr = TRUE
  )

  gcdf <- function(x) stats::pgamma(x, shape = shape, scale = scale)
  expect_true(any(gcdf(x + 1) - gcdf(x) == 0)) # naive underflows

  gccdf <- function(x) {
    stats::pgamma(x, shape = shape, scale = scale, lower.tail = FALSE)
  }
  expected_pmf <- gccdf(x) - gccdf(x + 1)

  expect_true(all(is.finite(distr$log_pdf(x))))
  expect_equal(distr$pdf(x), expected_pmf)
})

test_that("rand() draws are reproducible and statistically plausible", {
  d <- browse_distributions("normal", mean = 2, sd = 3)
  set.seed(42L)
  x1 <- d$rand(50L)
  set.seed(42L)
  x2 <- d$rand(50L)
  expect_identical(x1, x2) # same R seed -> same draws
  expect_length(x1, 50L)
  # different seed -> different draws
  set.seed(43L)
  expect_false(identical(x1, d$rand(50L)))

  # mean / sd of 20k draws should be close to the truth
  set.seed(7L)
  big <- d$rand(20000L)
  expect_equal(mean(big), 2, tolerance = 0.05)
  expect_equal(sd(big), 3, tolerance = 0.05)

  # without set.seed, successive calls are non-deterministic
  expect_false(identical(d$rand(10L), d$rand(10L)))

  # invalid n
  d$rand(-1L) |>
    expect_error("non-negative")
})

test_that("rand() honours location shift (relocate_ Active branch)", {
  loc_val <- 5
  d <- browse_distributions("normal", mean = 0, sd = 1, loc = loc_val)
  expect_equal(d$pdf(loc_val), stats::dnorm(0))
  expect_equal(mean(d$rand(20000L)), loc_val, tolerance = 0.05)
})

test_that("rand() on truncated distribution lands inside the bounds", {
  d <- browse_distributions("normal", mean = 0, sd = 1, min = -1, max = 2)
  x <- d$rand(5000L)
  expect_true(all(x >= -1 & x <= 2))
})

test_that("quantile returns NA outside [0, 1]", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  expect_true(is.na(distr$quantile(-0.1)))
  expect_true(is.na(distr$quantile(1.5)))
  expect_true(all(is.na(distr$quantile(c(-2, -0.001, 1.001, 7.7)))))
})

test_that("quantile at boundaries returns +-Inf when domain is unbounded", {
  distr <- browse_distributions("normal", mean = 0, sd = 1)
  expect_equal(distr$quantile(0), -Inf)
  expect_equal(distr$quantile(1), Inf)
})

test_that("discretized normal log_pdf, cdf and log_ratio go through Center branches", {
  mu <- 5
  sd <- 2
  x <- seq(0, 10, 1)
  distr <- browse_distributions("normal", mean = mu, sd = sd, discr = TRUE)

  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  expected_pmf <- normcdf(x + 0.5) - normcdf(x - 0.5)

  # log_func_ Center branch
  expect_equal(distr$log_pdf(x), log(expected_pmf))
  # primitive_ Center branch (cdf is the primitive)
  expect_equal(distr$cdf(x), normcdf(x + 0.5))
  # log_ratio_ Center branch with distinct rounded values
  xn <- x + 1
  expect_equal(
    distr$log_ratio(x, xn),
    log(distr$pdf(xn) / distr$pdf(x))
  )
  # log_ratio when rounded values coincide should be 0
  expect_equal(distr$log_ratio(x, x + 0.1), rep(0, length(x)))
})

test_that("discretized lognormal log_pdf equals log(pdf) (Left branch)", {
  # Regression: log_func_ Left case previously returned the raw PMF
  # instead of its log. log_ratio_ Left was already correct, so this
  # asymmetry is what made the bug invisible to MCMC delta computations.
  x <- seq(0, 8, 1)
  distr <- browse_distributions("lognormal",
    meanlog = 1, sdlog = 0.5,
    discr = TRUE
  )
  expect_equal(distr$log_pdf(x), log(distr$pdf(x)))
})

test_that("upper-only truncated normal covers Trunc::Max", {
  # only max set => discr_mode_ stays Center-eligible (no min => Center
  # if discr=TRUE) and trunc_mode_ becomes Max; here we test the
  # continuous case (discr_mode_ = None) which still exercises Trunc::Max
  # through cdf_upper_ and trunc_().
  mu <- 0
  sd <- 1
  max_val <- 1.5
  x <- seq(-3, 3, 0.5)
  distr <- browse_distributions("normal", mean = mu, sd = sd, max = max_val)

  normpdf <- function(x) stats::dnorm(x, mean = mu, sd = sd)
  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  zz <- normcdf(max_val)

  expec_pdf <- ifelse(x <= max_val, normpdf(x), 0) / zz
  expec_cdf <- ifelse(x > max_val, 1, normcdf(x) / zz)

  expect_equal(distr$pdf(x), expec_pdf)
  expect_equal(distr$cdf(x), expec_cdf)
})

test_that("discretized normal with upper-only truncation exercises Center cdf_upper_", {
  mu <- 2
  sd <- 1.5
  max_val <- 4
  x <- seq(-2, 6, 1)
  distr <- browse_distributions("normal",
    mean = mu, sd = sd,
    max = max_val, discr = TRUE
  )

  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  # Center discretization with truncation at max:
  # normalizer = cdf_upper_() = pnorm(floor(max) + 0.5)
  zz <- normcdf(floor(max_val) + 0.5)
  expected_pmf <- (normcdf(x + 0.5) - normcdf(x - 0.5)) / zz
  expected_pmf[x > max_val] <- 0

  expect_equal(distr$pdf(x), expected_pmf)
})

test_that("quantile matches R reference for standard distributions", {
  p <- seq(0, 1, 0.1)
  cases <- list(
    list(
      args = list("normal", mean = 2, sd = 3),
      ref = function(p) stats::qnorm(p, 2, 3)
    ),
    list(
      args = list("lognormal", meanlog = 1, sdlog = 0.5),
      ref = function(p) stats::qlnorm(p, 1, 0.5)
    ),
    list(
      args = list("gamma", shape = 2, scale = 1.5),
      ref = function(p) stats::qgamma(p, shape = 2, scale = 1.5)
    ),
    list(
      args = list("exponential", rate = 0.3),
      ref = function(p) stats::qexp(p, 0.3)
    ),
    list(
      args = list("beta", shape1 = 2, shape2 = 3),
      ref = function(p) stats::qbeta(p, 2, 3)
    ),
    list(
      args = list("poisson", lambda = 5),
      ref = function(p) stats::qpois(p, 5)
    ),
    list(
      args = list("binomial", size = 10, prob = 0.4),
      ref = function(p) stats::qbinom(p, 10, 0.4)
    ),
    list(
      args = list("geometric", prob = 0.3),
      ref = function(p) stats::qgeom(p, 0.3)
    ),
    list(
      args = list("negative_binomial", size = 5, prob = 0.4),
      ref = function(p) stats::qnbinom(p, 5, 0.4)
    )
  )
  for (case in cases) {
    distr <- do.call(browse_distributions, case$args)
    expect_equal(distr$quantile(p), case$ref(p), info = case$args[[1]])
  }
})

test_that("quantile is truncation-aware", {
  p <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  cases <- list(
    list(
      args = list("normal", mean = 0, sd = 1, min = 0),
      ref = function(p) qnorm(pnorm(0) + p * (1 - pnorm(0)))
    ),
    list(
      args = list("normal", mean = 3, sd = 2, min = 1, max = 5),
      ref = function(p) qnorm(pnorm(1, 3, 2) + p * (pnorm(5, 3, 2) - pnorm(1, 3, 2)), 3, 2)
    ),
    list(
      args = list("exponential", rate = 1, max = 3),
      ref = function(p) qexp(p * pexp(3))
    ),
    list(
      args = list("gamma", shape = 2, scale = 1.5, min = 1, max = 8),
      ref = function(p) {
        qgamma(
          pgamma(1, 2, scale = 1.5) +
            p * (pgamma(8, 2, scale = 1.5) - pgamma(1, 2, scale = 1.5)),
          2,
          scale = 1.5
        )
      }
    )
  )
  for (case in cases) {
    distr <- do.call(browse_distributions, case$args)
    expect_equal(distr$quantile(p), case$ref(p),
      tolerance = 1e-10, info = case$args[[1]]
    )
    expect_true(all(distr$quantile(p) >= case$args$min %||% -Inf),
      info = paste(case$args[[1]], "above min")
    )
    expect_true(all(distr$quantile(p) <= case$args$max %||% Inf),
      info = paste(case$args[[1]], "below max")
    )
  }
})

test_that("Beta rejects discrete option", {
  browse_distributions("beta", shape1 = 2, shape2 = 3, discr = TRUE) |>
    expect_error("invalid.*discrete")
})

test_that("Beta rejects out-of-domain truncation bounds", {
  browse_distributions("beta", shape1 = 2, shape2 = 3, min = -0.5) |>
    expect_error("min")
  browse_distributions("beta", shape1 = 2, shape2 = 3, max = 1.5) |>
    expect_error("max")
})

test_that("Beta truncation normalises correctly", {
  low <- 0.1
  up <- 0.8
  a <- 2
  b <- 5
  x <- seq(0, 1, 0.05)
  d <- browse_distributions("beta", shape1 = a, shape2 = b, min = low, max = up)
  ref_d <- function(x) stats::dbeta(x, a, b)
  ref_p <- function(x) stats::pbeta(x, a, b)
  zz <- ref_p(up) - ref_p(low)
  expect_equal(d$pdf(x), ifelse(x >= low & x <= up, ref_d(x), 0) / zz)
})

test_that("Beta rand() produces plausible draws", {
  a <- 2
  b <- 5
  d <- browse_distributions("beta", shape1 = a, shape2 = b)
  set.seed(99L)
  x1 <- d$rand(50L)
  set.seed(99L)
  x2 <- d$rand(50L)
  expect_identical(x1, x2)
  big <- d$rand(2000L)
  expect_true(all(big >= 0 & big <= 1))
  expect_equal(mean(big), a / (a + b), tolerance = 0.02)
})

test_that("Gamma pdf and log_pdf handle x == 0 across shape regimes", {
  scale <- 2
  # test case: shape < 1
  d1 <- browse_distributions("gamma", shape = 0.5, scale = scale)
  expect_equal(d1$pdf(0), Inf)
  expect_equal(d1$log_pdf(0), Inf)

  # test case: shape == 1
  d2 <- browse_distributions("gamma", shape = 1, scale = scale)
  expect_equal(d2$pdf(0), 1 / scale)
  expect_equal(d2$log_pdf(0), -log(scale))

  # test case: shape > 1
  d3 <- browse_distributions("gamma", shape = 2, scale = scale)
  expect_equal(d3$pdf(0), 0)
  expect_equal(d3$log_pdf(0), -Inf)
})

test_that("Binomial rejects out-of-domain truncation bounds", {
  browse_distributions("binomial", size = 10, prob = 0.5, min = -1) |>
    expect_error("min")
  browse_distributions("binomial", size = 10, prob = 0.5, min = 11) |>
    expect_error("max")
})

test_that("rand() is reproducible and plausible for all distributions", {
  cases <- list(
    list(
      args = list("rayleigh", scale = 2),
      mean = 2 * sqrt(pi / 2), int = FALSE
    ),
    list(
      args = list("exponential", rate = 0.5),
      mean = 2, int = FALSE
    ),
    list(
      args = list("poisson", lambda = 5),
      mean = 5, int = TRUE
    ),
    list(
      args = list("binomial", size = 20, prob = 0.3),
      mean = 6, int = TRUE, max = 20
    ),
    list(
      args = list("negative_binomial", size = 5, prob = 0.4),
      mean = 5 * 0.6 / 0.4, int = TRUE
    ),
    list(
      args = list("geometric", prob = 0.3),
      mean = 0.7 / 0.3, int = TRUE
    )
  )
  for (cc in cases) {
    d <- do.call(browse_distributions, cc$args)
    nm <- cc$args[[1]]
    set.seed(42L)
    x1 <- d$rand(50L)
    set.seed(42L)
    expect_identical(d$rand(50L), x1, info = nm)
    set.seed(1L)
    big <- d$rand(2000L)
    expect_true(all(big >= 0), info = nm)
    if (cc$int)
      expect_true(all(big == floor(big)), info = nm)
    if (!is.null(cc$max))
      expect_true(all(big <= cc$max), info = nm)
    expect_equal(
      mean(big), cc$mean,
      tolerance = 0.15, info = nm
    )
  }
})

############################################
## BROWSE INTERNAL : INFECTIVITY PROFILES ##
############################################

test_that("Ferretti20 infectivity profile behave has expected", {
  x <- seq(0, 15, 0.001)
  profile <- browse_infectivity_profiles("ferretti20", incub_period = 5)
  expect_equal(trapezoid(x, profile$inst(x)), profile$cumul(x)[length(x)])

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  profile <- browse_infectivity_profiles("ferretti20", incub_period = 5, discr = TRUE)
  expect_true(all(profile$inst(x) == profile$cumul(1)))
})

####################################
## BROWSE INTERNAL : PRINT METHOD ##
####################################

test_that("Exposed object print behavior is stable", {
  expect_snapshot(print(browse_functions("pow", x2 = 3)))
  expect_snapshot(print(browse_distributions("normal", mean = 2, sd = 3)))
  expect_snapshot(print(browse_infectivity_profiles("gamma", shape = 2, scale = 1)))
  # no-argument object prints only the header
  expect_snapshot(print(browse_functions("sin")))
})
