library(testthat)

test_that("Model class inconsistent calls", {
  # blank household model object
  hmc <- HouseholdModelClass$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))

  # data
  expect_error(hmc$dataset,
               "Impossible to export data before model build.")
  # traces
  expect_error(hmc$traces,
               "Impossible to export MCMC traces before model inference.")
  expect_error(hmc$traces_results,
               "Impossible to export MCMC traces before model inference.")

  # household traces
  expect_error(hmc$loo(),
               "Impossible to run method 'loo' before model inference.")
  expect_error(hmc$waic(),
               "Impossible to run method 'waic' before model inference.")

  # inference without definitions
  expect_error(
    hmc$infer(),
    "Nothing to sample. Check your parameter and augmented variable definitions."
  )
})

test_that("Model class blank calls", {

  expect_warning(HouseholdModelClass$new(),
                 "No model specified. Falling back to \"blank\" model.")

  # blank household model object
  hmc <- HouseholdModelClass$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))

  # templates
  expect_warning(hmc$load_param_template(),
                 "No template parameters available for blank model.")
  expect_warning(hmc$load_augm_template(),
                 "No template augmented variables available for blank model.")

  # zero loglikelihood
  expect_warning(hmc$compute_loglikelihood(),
                 "Computing the likelihood with \"blank\" model.")
  expect_equal(hmc$compute_loglikelihood(), 0)

  # blank simulation
  expect_error(hmc$simulate(),
               "Impossible to perform simulation with \"blank\" model.")
})

test_that("Numerical recipes are correct", {

  hmc <- housepi::HouseholdModelClass$new("blank")
  withr::defer(rm(hmc))

  # test binomial cdf #
  p <- list(size = 12, prob = 0.2)
  hmc$add_infectivity_profile(name = "binom",
                              args = list(size = p$size, prob = p$prob))
  out <- hmc$export_infectivity_profile()
  # instantaneous
  y <- stats::dbinom(out[[1]], size = p$size, prob = p$prob)
  expect_equal(out[[2]], y)
  # cumulative
  y <- stats::pbinom(out[[1]], size = p$size, prob = p$prob)
  expect_equal(out[[3]], y)

  # test regularised incomplete beta function #
  p <- list(shape1 = 0.02, shape2 = 105.6)
  hmc$add_infectivity_profile(name = "beta",
                              args = list(shape1 = p$shape1, shape2 = p$shape2))
  out <- hmc$export_infectivity_profile()
  # instantaneous
  y <- stats::dbeta(out[[1]], shape1 = p$shape1, shape2 = p$shape2)
  expect_equal(out[[2]], y)
  # cumulative
  y <- stats::pbeta(out[[1]], shape1 = p$shape1, shape2 = p$shape2)
  expect_equal(out[[3]], y)

  # test regularised incomplete gamma function #
  p <- list(shape = 0.5, scale = 10.3)
  hmc$add_infectivity_profile(name = "gamma",
                              args = list(shape = p$shape, scale = p$scale))
  out <- hmc$export_infectivity_profile()
  # instantaneous
  y <- stats::dgamma(out[[1]], shape = p$shape, scale = p$scale)
  expect_equal(out[[2]], y)
  # cumulative
  y <- stats::pgamma(out[[1]], shape = p$shape, scale = p$scale)
  expect_equal(out[[3]], y)
})

test_that("Chain binomial models are coherent", {
  # loglikelihood from chain binomial model
  hmc1 <- HouseholdModelClass$new("chain-binomial")
  withr::defer(rm(hmc1))
  hmc1$load_param_list(alpha = 0.15, beta = 0.22)
  ll1 <- hmc1$compute_loglikelihood(housepi::survey)

  # loglikelihood from heterogeneous chain binomial model
  # with homogenous rates
  hmc2 <- HouseholdModelClass$new("chain-binomial-heterogeneous")
  withr::defer(rm(hmc2))

  hmc2$load_param_list(alpha = 0.15, beta = 0.22,
                       inf_child = 1., com_sus_child = 1.,
                       com_sus_woman = 1.0, sus_child = 1.0,
                       sus_woman = 1., child_age = 15)
  ll2 <- hmc2$compute_loglikelihood(housepi::survey)

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
  hmc <- HouseholdModelClass$new("eg-augm-incubation-from-delays")
  withr::defer(rm(hmc))
  hmc$load_param_list(shape = shape, scale = scale)
  hmc$load_augm_list(exposure = "t1@D", onset = "t2@D")

  expect_warning(
    from_model <- hmc$compute_loglikelihood(dataset),
    "Unknown household structure: missing data column \"house_id\"."
  )

  expect_equal(from_model, exact)
})
