library(testthat)

test_that("Function exported from housepi are correct", {
  x <- seq(0, 1000, 1)
  expect_equal(expose_function("pow")$func(x), x)
  expect_equal(expose_function("pow", x2 = 2)$func(x), x^2)
  expect_equal(expose_function("pow", x2 = -3)$func(x), x^(-3))
  expect_equal(expose_function("exp")$func(x), exp(x))
  expect_equal(expose_function("log")$func(x), log(x))
  expect_equal(expose_function("sin")$func(x), sin(x))
  expect_equal(expose_function("cos")$func(x), cos(x))
  expect_equal(expose_function("sin")$func(x), sin(x))
  expect_equal(expose_function("tan")$func(x), tan(x))

  args <- c(3.3, 5.5, 6.6)
  expect_equal(expose_function("sum", args)$func(0), sum(args))
  expect_equal(expose_function("prod", args)$func(1), prod(args))
  expect_equal(expose_function("min", args)$func(2), 2)
  expect_equal(expose_function("max", args)$func(2), 6.6)
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

test_that("Ferretti20 infectivity profile behave has expected", {
  x <- seq(0, 15, 0.001)
  profile <- expose_infec_profile("ferretti20")

  trapezoid <- function(x, y) {
    0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))
  }
  expect_equal(
    trapezoid(x, profile$inst(x)),
    profile$cumul(x)[length(x)]
  )

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  profile <- expose_infec_profile("ferretti20", discr = TRUE)
  expect_true(all(profile$inst(x) == profile$cumul(1)))
})
