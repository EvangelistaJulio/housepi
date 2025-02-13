library(testthat)

test_that("Model class inconsistent calls", {
  # blank household model object
  hmc <- HouseholdModelClass$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))
  
  # data
  expect_error(hmc$export_data(), 
               "Impossible to export data before model build.")
  expect_error(hmc$summarise_data(), 
               "Impossible to show data summary before model build.")
  # traces
  expect_error(hmc$export_traces(), 
               "Impossible to export MCMC traces before model inference.")
  expect_error(hmc$summarise_traces(), 
               "Impossible to export MCMC traces before model inference.")

  # household traces
  expect_error(hmc$loo(), 
               "Impossible to run method 'loo' before model inference.")
  expect_error(hmc$waic(), 
               "Impossible to run method 'waic' before model inference.")
  
  # inference without definitions
  expect_error(hmc$infer(), 
               "Nothing to infer. Check your parameter and augmented variable definitions.")
})

test_that("Model class blank calls", {

  expect_warning(HouseholdModelClass$new(), 
                 "No model specified. Falling back to \"blank\" model.")

  # blank household model object
  hmc <- HouseholdModelClass$new("blank")
  # cleanup at end of test
  withr::defer(rm(hmc))

  # templates
  expect_warning(hmc$template_param(), 
                 "No template parameters available for blank model.")
  expect_warning(hmc$template_augm(), 
                 "No template augmented variables available for blank model.")
  
  # zero loglikelihood
  expect_warning(hmc$loglikelihood(), 
                 "Computing the likelihood with \"blank\" model.")
  expect_equal(hmc$loglikelihood(), 0)
  
  # blank simulation
  expect_error(hmc$simulate(), 
               "Impossible to perform simulation with \"blank\" model.")
})

test_that("Chain binomial models are coherent", {
  # loglikelihood from chain binomial model
  hmc1 <- HouseholdModelClass$new("chain-binomial")
  withr::defer(rm(hmc1))
  ll1 <- hmc1$load_param_list(list(
    alpha = 0.15, beta = 0.22
  ))$loglikelihood(survey)

  # loglikelihood from heterogeneous chain binomial model 
  # with homogenous rates
  hmc2 <- HouseholdModelClass$new("chain-binomial-heterogeneous")
  withr::defer(rm(hmc2))

  ll2 <- hmc2$load_param_list(list(
      alpha = 0.15, beta = 0.22, 
      inf_child = 1., com_sus_child = 1.,
      com_sus_woman = 1.0, sus_child = 1.0, 
      sus_woman = 1., child_age = 15
  ))$loglikelihood(survey)

  expect_equal(ll1, ll2)
})