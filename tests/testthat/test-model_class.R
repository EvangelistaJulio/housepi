library(testthat)

test_that("Model class inconsistent calls", {
  # blank household model object
  hmc <- HouseholdModel$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))

  # data
  expect_error(hmc$dataset)
  # traces
  expect_error(hmc$traces)
  expect_error(hmc$traces_results)

  # household traces
  expect_error(hmc$loo())
  expect_error(hmc$waic())

  # inference without definitions
  expect_error(hmc$infer())
})

test_that("Model class blank calls", {

  expect_warning(HouseholdModel$new())

  # blank household model object
  hmc <- HouseholdModel$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))

  # templates
  expect_warning(hmc$load_param_template())
  expect_warning(hmc$load_augm_template())

  # zero loglikelihood
  expect_warning(x <- hmc$loglikelihood())
  expect_equal(x, 0)

  # blank simulation
  expect_error(hmc$simulate())
})

test_that("Numerical recipes and distributions are correct", {

  hmc <- housepi::HouseholdModel$new("blank")
  withr::defer(rm(hmc))

  # test regularised incomplete beta function #
  p <- list(shape1 = 0.02, shape2 = 105.6)
  hmc$add_infec_profile(name = "beta",
                        args = c(shape1 = p$shape1, shape2 = p$shape2))
  profile <- hmc$showcase_infec_profile()$profile
  # instantaneous
  y <- stats::dbeta(profile$gentime, shape1 = p$shape1, shape2 = p$shape2)
  expect_equal(profile$inst, y)
  # cumulative
  y <- stats::pbeta(profile$gentime, shape1 = p$shape1, shape2 = p$shape2)
  expect_equal(profile$cumul, y)

  # test regularised incomplete gamma function #
  p <- list(shape = 0.5, scale = 10.3)
  hmc$add_infec_profile(name = "gamma",
                        args = c(shape = p$shape, scale = p$scale))
  profile <- hmc$showcase_infec_profile()$profile
  # instantaneous
  y <- stats::dgamma(profile$gentime, shape = p$shape, scale = p$scale)
  expect_equal(profile$inst, y)
  # cumulative
  y <- stats::pgamma(profile$gentime, shape = p$shape, scale = p$scale)
  expect_equal(profile$cumul, y)

  # test binomial cdf #
  # discrete infectivity starts from 1: test truncation
  p <- list(size = 12, prob = 0.2)
  hmc$add_infec_profile(name = "binom",
                              args = c(size = p$size, prob = p$prob))
  profile <- hmc$showcase_infec_profile()$profile
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
  hmc1 <- HouseholdModel$new("chain-binomial")
  withr::defer(rm(hmc1))
  hmc1$load_param_list(alpha = 0.15, beta = 0.22)
  ll1 <- hmc1$loglikelihood(housepi::survey)

  # loglikelihood from heterogeneous chain binomial model
  # with homogenous rates
  hmc2 <- HouseholdModel$new("chain-binomial-heterogeneous")
  withr::defer(rm(hmc2))

  hmc2$load_param_list(alpha = 0.15, beta = 0.22)
  ll2 <- hmc2$loglikelihood(housepi::survey)

  expect_equal(ll1, ll2)
})

test_that("Augmentation example on incubation from delays is coherent", {
  # loglikelihood from simulated delays
  shape <- 2.0
  scale <- 1.1
  dataset <- data.frame(t1 = 0,
                        t2 = rgamma(223, shape = shape, scale = scale))
  exact <- sum(log(stats::dgamma(dataset$t2, shape = shape, scale = scale)))

  # generate an instance of the household model class
  hmc <- HouseholdModel$new("augm-incubation-from-delays")
  withr::defer(rm(hmc))
  hmc$load_param_list(shape = shape, scale = scale)
  hmc$load_augm_list(exposure = "t1@D", onset = "t2@D")

  expect_warning(
    from_model <- hmc$loglikelihood(dataset),
    "Unknown household structure: missing data column \"house_id\"."
  )

  expect_equal(from_model, exact)
})
