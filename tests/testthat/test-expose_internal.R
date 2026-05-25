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
  expect_equal(
    fn$log_ratio(x, x + 0.1),
    log(fn$func(x + 0.1) / fn$func(x))
  )
  expect_error(fn$log_ratio(x, x[-1]), regexp = "[Ss]ize")

  # x2 = -1 primitive uses log(|x|)
  expect_equal(expose_function("pow", x2 = -1)$prim(x), log(abs(x)))
  # x2 = 0.5 equals sqrt
  expect_equal(expose_function("pow", x2 = 0.5)$func(x), sqrt(x))

  # error cases
  expect_error(expose_function("pow")$func(x), regexp = "Missing non default")
  expect_error(expose_function("pow", 5), regexp = "Unnamed argument")
  expect_error(expose_function("pow", discr = TRUE)$func(x),
    regexp = "undefined method"
  )
})

test_that("Square root function is correct", {
  x <- seq(0.1, 10, 0.1)
  fn <- expose_function("sqrt")
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
  fn <- expose_function("exp")
  expect_equal(fn$func(x), exp(x))
  expect_equal(fn$prim(x), exp(x))
  expect_equal(fn$log_func(x), x)
  expect_equal(
    fn$log_ratio(x, x + pi),
    log(fn$func(x + pi) / fn$func(x))
  )

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
  expect_equal(
    fn$log_ratio(x, x + pi),
    log(fn$func(x + pi) / fn$func(x))
  )
})

test_that("Sine function is correct", {
  x <- seq(0.1, pi - 0.1, 0.05)
  fn <- expose_function("sin")
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
    expose_function("sin", loc = -0.5 * pi)$func(x_full),
    cos(x_full)
  )
})

test_that("Cosine function is correct", {
  x <- seq(-pi / 2 + 0.1, pi / 2 - 0.1, 0.05)
  fn <- expose_function("cos")
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
  fn <- expose_function("tan")
  expect_equal(fn$func(x), tan(x))
  expect_equal(fn$log_func(x), log(tan(x)))
  expect_equal(
    fn$log_ratio(x, x + 0.1),
    log(fn$func(x + 0.1) / fn$func(x))
  )
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

# Assert pdf, cdf, log_pdf and log_ratio of an exposed distribution against
# user-supplied R reference d/p functions.
expect_distribution_matches <- function(name, args = list(), ref_d, ref_p,
                                        x = seq(-20, 20, 1), h = pi) {
  distr <- do.call(expose_distribution, c(list(name), args))
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
  distr <- expose_distribution("binomial", list(size = size, prob = prob),
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

test_that("Default parameters match R reference for d", {
  x <- seq(-20, 20, 1)
  expect_equal(expose_distribution("normal")$pdf(x), stats::dnorm(x))
  expect_equal(expose_distribution("uniform")$pdf(x), stats::dunif(x))
  expect_equal(expose_distribution("lognormal")$pdf(x), stats::dlnorm(x))
  expect_equal(expose_distribution("exponential")$pdf(x), stats::dexp(x))
})

test_that("Normal both-sided truncation is correct", {
  low <- 0
  up <- 5
  x <- seq(-20, 20, 1)
  distr <- expose_distribution("normal", mean = 2, sd = 3, min = low, max = up)
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
  expect_error(expose_distribution("gamma", mean = 2, scale = 1.1), "mix")
  expect_error(expose_distribution("gamma", shape = 2, var = 1.1), "mix")

  m <- 2.5
  v <- 1.7
  shape <- m * m / v
  scale <- v / m
  x <- seq(0.1, 12, 0.3)
  d_can <- expose_distribution("gamma", shape = shape, scale = scale)
  d_alt <- expose_distribution("gamma", mean = m, var = v)
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
    tolerance = 1e-3
  )

  # lognormal: integrate pdf from 0 to 20, expect cdf(20)
  x <- seq(0.001, 20, length.out = 3000)
  distr <- expose_distribution("lognormal", meanlog = 1, sdlog = 0.5)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )

  # rayleigh: integrate pdf from 0 to 12, expect cdf(12)
  x <- seq(0.001, 12, length.out = 3000)
  distr <- expose_distribution("rayleigh", scale = 2)
  expect_equal(trapezoid(x, distr$pdf(x)), distr$cdf(tail(x, 1)),
    tolerance = 1e-3
  )
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
  distr <- expose_distribution("lognormal",
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
  distr <- expose_distribution("gamma",
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
  distr <- expose_distribution("normal", mean = mu, sd = sd, discr = TRUE)

  # center discretization: pdf(k) = pnorm(k+0.5) - pnorm(k-0.5)
  normcdf <- function(x) stats::pnorm(x, mean = mu, sd = sd)
  expect_equal(distr$pdf(x), normcdf(x + 0.5) - normcdf(x - 0.5))
})

test_that("rand() draws are reproducible and statistically plausible", {
  d <- expose_distribution("normal", mean = 2, sd = 3)
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
  expect_error(d$rand(-1L), regexp = "non-negative")
})

test_that("rand() honours location shift (relocate_ Active branch)", {
  loc_val <- 5
  d <- expose_distribution("normal", mean = 0, sd = 1, loc = loc_val)
  expect_equal(d$pdf(loc_val), stats::dnorm(0))
  set.seed(11L)
  expect_equal(mean(d$rand(20000L)), loc_val, tolerance = 0.05)
})

test_that("rand() on truncated distribution lands inside the bounds", {
  d <- expose_distribution("normal", mean = 0, sd = 1, min = -1, max = 2)
  set.seed(5L)
  x <- d$rand(5000L)
  expect_true(all(x >= -1 & x <= 2))
})

test_that("quantile returns NA outside [0, 1]", {
  distr <- expose_distribution("normal", mean = 0, sd = 1)
  expect_true(is.na(distr$quantile(-0.1)))
  expect_true(is.na(distr$quantile(1.5)))
  expect_true(all(is.na(distr$quantile(c(-2, -0.001, 1.001, 7.7)))))
})

test_that("quantile at boundaries returns +-Inf when domain is unbounded", {
  distr <- expose_distribution("normal", mean = 0, sd = 1)
  expect_equal(distr$quantile(0), -Inf)
  expect_equal(distr$quantile(1), Inf)
})

test_that("discretized normal log_pdf, cdf and log_ratio go through Center branches", {
  mu <- 5
  sd <- 2
  x <- seq(0, 10, 1)
  distr <- expose_distribution("normal", mean = mu, sd = sd, discr = TRUE)

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
  distr <- expose_distribution("lognormal",
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
  distr <- expose_distribution("normal", mean = mu, sd = sd, max = max_val)

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
  distr <- expose_distribution("normal",
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
    distr <- do.call(expose_distribution, case$args)
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
    distr <- do.call(expose_distribution, case$args)
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
