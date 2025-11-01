library(testthat)

test_that("Function exported from housepi are correct", {

  x <- seq(0, 1000, 1)

  expect_equal(test_function("pow", x), x)
  expect_equal(test_function("pow", x, args = c(x2 = 2)), x^2)
  expect_equal(test_function("pow", x, args = c(x2 = -3)), x^(-3))
  expect_equal(test_function("exp", x), exp(x))
  expect_equal(test_function("log", x), log(x))
  expect_equal(test_function("sin", x), sin(x))
  expect_equal(test_function("cos", x), cos(x))
  expect_equal(test_function("sin", x), sin(x))
  expect_equal(test_function("tan", x), tan(x))
})

test_that("Distributions exported from housepi are correct", {
  x <- seq(0, 1000, 1)
  r <- unique(c(-x, x))

  # Test normal distribution
  y <- test_distribution("normal", r)
  expect_equal(y, dnorm(r))

  # test R-like prefix
  expect_equal(test_distribution("normal", r), y)

  # Test uniform distribution
  y <- test_distribution("uniform", x, args = c(min = 0, max = 1))
  expect_equal(y, dunif(x, min = 0, max = 1))

  # Test lognormal distribution
  y <- test_distribution("lognormal", x, args = c(meanlog = 0, sdlog = 1))
  expect_equal(y, dlnorm(x, meanlog = 0, sdlog = 1))

  # Test beta distribution
  y <- test_distribution("beta", x, args = c(shape1 = 2, shape2 = 5))
  expect_equal(y, dbeta(x, shape1 = 2, shape2 = 5))

  # Test gamma distribution
  y <- test_distribution("gamma", x, args = c(shape = 2, scale = 2))
  expect_equal(y, dgamma(x, shape = 2, scale = 2))

  # Test exponential distribution
  y <- test_distribution("exponential", x, args = c(rate = 1))
  expect_equal(y, dexp(x, rate = 1))

  # Test poisson distribution
  y <- test_distribution("poisson", x, args = c(lambda = 5))
  expect_equal(y, dpois(x, lambda = 5))

  # Test binomial distribution
  y <- test_distribution("binomial", x, args = c(size = 10, prob = 0.5))
  expect_equal(y, dbinom(x, size = 10, prob = 0.5))

  # Test negative binomial distribution
  y <- test_distribution("negative_binomial", x,
                         args = c(size = 10, prob = 0.5))
  expect_equal(y, dnbinom(x, size = 10, prob = 0.5))

  # Test geometric distribution
  y <- test_distribution("geometric", x, args = c(proba = 0.5))
  expect_equal(y, dgeom(x, prob = 0.5))

  # Test weibull distribution
  y <- test_distribution("weibull", x, args = c(shape = 2, scale = 1))
  expect_equal(y, dweibull(x, shape = 2, scale = 1))
})

test_that("Infectivity profile behave has expected", {
  x <- seq(0, 15, 0.001)
  y <- test_infec_profile("ferretti20", x)

  trapezoid <- function(x, y) {
    0.5 * sum(diff(x) * (y[-length(x)] + y[-1]))
  }
  z <- test_infec_profile("ferretti20", x, cumul = TRUE)
  expect_equal(trapezoid(x, y), z[length(x)])

  # left discretization works as expected
  x <- seq(1.0, 1.9, 0.1)
  y <- test_infec_profile("ferretti20", x, discr = TRUE)
  z <- test_infec_profile("ferretti20", 1, cumul = TRUE, discr = TRUE)
  expect_true(all(z == z))
})

test_that("Numerical recipes and distributions are correct", {
  x <- seq(1, 10, 1)
  args <- list(size = 12, prob = 0.2)
  prof_inst <- test_distribution("binomial", x,
                                 args = args, min = 1)
  prof_cumul <- test_distribution("binomial", x,
                                  args = args, cumul = TRUE, min = 1)

  # truncation
  shift <- stats::pbinom(0, size = args$size, prob = args$prob)
  trunc <- 1. / (1. - shift)
  # instantaneous
  y <- stats::dbinom(x, size = args$size, prob = args$prob)
  expect_equal(prof_inst, y * trunc)

  # cumulative
  y <- stats::pbinom(x, size = args$size, prob = args$prob)
  expect_equal(prof_cumul, (y - shift) * trunc)
})
