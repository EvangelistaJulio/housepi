library(testthat)

test_that("Function exported from housepi are correct", {

  x <- seq(0, 1000, 1)
  
  expect_equal(housepi_function("pow", x), x)
  expect_equal(housepi_function("pow", x, args = c(y = 2)), x^2)
  expect_equal(housepi_function("pow", x, args = c(y = -3)), x^(-3))
  expect_equal(housepi_function("exp", x), exp(x))
  expect_equal(housepi_function("log", x), log(x))
  expect_equal(housepi_function("sin", x), sin(x))
  expect_equal(housepi_function("cos", x), cos(x))
  expect_equal(housepi_function("sin", x), sin(x))
  expect_equal(housepi_function("tan", x), tan(x))
})
  
test_that("Distributions exported from housepi are correct", {
  x <- seq(0, 1000, 1)
  xR <- unique(c(-x, x))
  
  # Test normal distribution
  y <- housepi_distribution("norm", xR)
  expect_equal(y, dnorm(xR))

  # test R-like prefix
  expect_equal(housepi_distribution("dnorm", xR), y)
  
  # Test uniform distribution
  y <- housepi_distribution("unif", x, args = c(min = 0, max = 1))
  expect_equal(y, dunif(x, min = 0, max = 1))
  
  # Test lognormal distribution
  y <- housepi_distribution("lnorm", x, args = c(meanlog = 0, sdlog = 1))
  expect_equal(y, dlnorm(x, meanlog = 0, sdlog = 1))
  
  # Test beta distribution
  y <- housepi_distribution("beta", x, args = c(shape1 = 2, shape2 = 5))
  expect_equal(y, dbeta(x, shape1 = 2, shape2 = 5))
  
  # Test gamma distribution
  y <- housepi_distribution("gamma", x, args = c(shape = 2, scale = 2))
  expect_equal(y, dgamma(x, shape = 2, scale = 2))
  
  # Test exponential distribution
  y <- housepi_distribution("exp", x, args = c(rate = 1))
  expect_equal(y, dexp(x, rate = 1))
  
  # Test poisson distribution
  y <- housepi_distribution("pois", x, args = c(lambda = 5))
  expect_equal(y, dpois(x, lambda = 5))
  
  # Test binomial distribution
  y <- housepi_distribution("binom", x, args = c(size = 10, prob = 0.5))
  expect_equal(y, dbinom(x, size = 10, prob = 0.5))
  
  # Test negative binomial distribution
  y <- housepi_distribution("nbinom", x, args = c(size = 10, prob = 0.5))
  expect_equal(y, dnbinom(x, size = 10, prob = 0.5))
  
  # Test geometric distribution
  y <- housepi_distribution("geom", x, args = c(proba = 0.5))
  expect_equal(y, dgeom(x, prob = 0.5))
  
  # Test weibull distribution
  y <- housepi_distribution("weibull", x, args = c(shape = 2, scale = 1))
  expect_equal(y, dweibull(x, shape = 2, scale = 1))
})

test_that("Ferretti infectivity profile is correctly normalized", {
  x <- seq(0, 15, 0.001)
  y <- housepi_infec_profile("ferretti20", x)
  
  trapezoid <- function(x, y) {
    if (length(x) < 2) stop("Need at least 2 points")
    heights <- 0.5 * (y[-length(x)] + y[-1])
    return(sum(diff(x) * heights))
  }
  z <- housepi_infec_profile("pferretti20", x)
  expect_equal(trapezoid(x, y), z[length(x)])
})

test_that("String function extraction returns expected results", {
  x <- housepi:::extract_func_("=alpha@P")
  expect_equal(x$func, "eq")
  expect_equal(x$key, "x")
  expect_equal(x$val, "alpha")
  expect_equal(x$at, "P")

  expect_equal(housepi:::extract_func_("alpha@P+beta@P")$val,
               housepi:::type_join_(c("alpha", "beta")))

  x <- housepi:::extract_func_("log ( alpha @ A ) ")
  expect_equal(x$func, "log")
  expect_equal(x$key, "x")
  expect_equal(x$val, "alpha")
  expect_equal(x$at, "A")

  expect_error(housepi:::extract_augm_func_("=alpha@P+"))
  expect_error(housepi:::extract_func_("log ( + alpha@P  )"))
  expect_error(housepi:::extract_func_("arc ( alpha@P  )"))
  expect_error(housepi:::extract_func_("log ( alpha@P  "))
  expect_error(housepi:::extract_func_("log alpha@P )"))
  expect_error(housepi:::extract_func_("log ( alpha @ B )  "))
})


test_that("Numerical recipes and distributions are correct", {
  x <- seq(1, 10, 1)
  args <- list(size = 12, prob = 0.2)
  prof_inst <- housepi_distribution("binom", x, args = args, min = 1)
  prof_cumul <- housepi_distribution("pbinom", x, args = args, min = 1)

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
