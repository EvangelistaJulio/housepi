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
  expect_error(hmc$infer(), 
               "Nothing to sample. Check your parameter and augmented variable definitions.")
})

test_that("Model class blank calls", {

  expect_warning(HouseholdModelClass$new(), 
                 "No model specified. Falling back to \"blank\" model.")

  # blank household model object
  hmc <- HouseholdModelClass$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))

  # templates
  expect_warning(hmc$add_template_param(), 
                 "No template parameters available for blank model.")
  expect_warning(hmc$add_template_augm(), 
                 "No template augmented variables available for blank model.")
  
  # zero loglikelihood
  expect_warning(hmc$compute_loglikelihood(), 
                 "Computing the likelihood with \"blank\" model.")
  expect_equal(hmc$compute_loglikelihood(), 0)
  
  # blank simulation
  expect_error(hmc$simulate(), 
               "Impossible to perform simulation with \"blank\" model.")
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

  hmc2$load_param_list(
      alpha = 0.15, beta = 0.22, 
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
  exact <- sum(log(dgamma(dataset$t2, shape = shape, scale = scale)))
  
  # generate an instance of the household model class
  hmc <- HouseholdModelClass$new("eg-augm-incubation-from-delays")
  hmc$load_param_list(shape = shape, scale = scale)
  hmc$load_augm_list(exposure = "t1@D", onset = "t2@D")
  
  expect_warning(from_model <- hmc$compute_loglikelihood(dataset), 
                 "Unknown household structure: missing data column \"house_id\".")

  expect_equal(exact, from_model)
})
