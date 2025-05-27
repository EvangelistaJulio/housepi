library(testthat)

test_that("Model class inconsistent calls", {
  # blank household model object
  hhmodel <- HouseholdModel$new("blank")
  # cleanup at end of test
  withr::defer(rm(hhmodel))

  # data
  expect_error(hhmodel$dataset)
  # traces
  expect_error(hhmodel$traces)
  expect_error(hhmodel$traces_results)

  # household traces
  expect_error(hhmodel$loo())
  expect_error(hhmodel$waic())

  # inference without definitions
  expect_error(hhmodel$infer())
})

test_that("Model class blank calls", {

  expect_warning(HouseholdModel$new())

  # blank household model object
  hhmodel <- HouseholdModel$new("blank")
  # cleanup at end of test
  withr::defer(rm(hhmodel))

  # templates
  expect_warning(hhmodel$load_param_template())
  expect_warning(hhmodel$load_augm_template())

  # zero loglikelihood
  x <- hhmodel$loglikelihood()
  expect_equal(x, 0)

  # blank simulation
  expect_error(hhmodel$simulate())
})

test_that("Numerical recipes and distributions are correct", {
  hhmodel <- housepi::HouseholdModel$new("blank")
  withr::defer(rm(hhmodel))

  # test regularised incomplete beta function #
  p <- list(shape1 = 0.02, shape2 = 105.6)
  hhmodel$add_infec_profile(name = "beta",
                        args = c(shape1 = p$shape1, shape2 = p$shape2))
  profile <- hhmodel$showcase_infec_profile()$profile
  # instantaneous
  y <- stats::dbeta(profile$gentime, shape1 = p$shape1, shape2 = p$shape2)
  expect_equal(profile$inst, y)
  # cumulative
  y <- stats::pbeta(profile$gentime, shape1 = p$shape1, shape2 = p$shape2)
  expect_equal(profile$cumul, y)

  # test regularised incomplete gamma function #
  p <- list(shape = 0.5, scale = 10.3)
  hhmodel$add_infec_profile(name = "gamma",
                            args = c(shape = p$shape, scale = p$scale))
  profile <- hhmodel$showcase_infec_profile()$profile
  # instantaneous
  y <- stats::dgamma(profile$gentime, shape = p$shape, scale = p$scale)
  expect_equal(profile$inst, y)
  # cumulative
  y <- stats::pgamma(profile$gentime, shape = p$shape, scale = p$scale)
  expect_equal(profile$cumul, y)

  # test binomial cdf #
  # discrete infectivity starts from 1: test truncation
  p <- list(size = 12, prob = 0.2)
  hhmodel$add_infec_profile(name = "binom",
                              args = c(size = p$size, prob = p$prob))
  profile <- hhmodel$showcase_infec_profile()$profile
  # truncation
  shift <- stats::pbinom(0, size = p$size, prob = p$prob)
  trunc <- 1. / (1. - shift)
  # instantaneous
  y <- stats::dbinom(profile$gentime, size = p$size, prob = p$prob)
  expect_equal(profile$inst, y * trunc)

  # cumulative
  y <- stats::pbinom(profile$gentime, size = p$size, prob = p$prob)
  expect_equal(profile$cumul, (y - shift) * trunc)
})

test_that("Chain binomial models are coherent", {
  # loglikelihood from chain binomial model
  hhmodel1 <- HouseholdModel$new("chain-binomial")
  withr::defer(rm(hhmodel1))
  hhmodel1$load_param_list(alpha = 0.15, beta = 0.22)
  ll1 <- hhmodel1$loglikelihood(housepi::survey)

  # loglikelihood from heterogeneous chain binomial model
  # with homogenous rates
  hhmodel2 <- HouseholdModel$new("chain-binomial-heterogeneous")
  withr::defer(rm(hhmodel2))

  hhmodel2$load_param_list(alpha = 0.15, beta = 0.22)
  ll2 <- hhmodel2$loglikelihood(housepi::survey)

  expect_equal(ll1, ll2)
})

test_that("Augmentation example on incubation from delays is coherent", {
  # loglikelihood from simulated delays
  shape <- 2.0
  scale <- 1.1
  dataset <- data.frame(status = 1,
                        ip = rgamma(223, shape = shape, scale = scale))
  exact <- sum(log(stats::dgamma(dataset$ip, shape = shape, scale = scale)))

  # generate an instance of the household model class
  hhmodel <- HouseholdModel$new("augm-example")
  withr::defer(rm(hhmodel))
  hhmodel$load_augm_list(
    proba = list(func = "gamma",
                 func_args = c(x = "ip@D", shape = shape, scale = scale)),
    log_proba = "log(proba@A)"
    )

  expect_equal(hhmodel$loglikelihood(dataset), exact)
})

