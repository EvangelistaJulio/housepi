library(testthat)

test_that("Able to extract results from coda", {
  expect_error(retrieve_traces_results(data.frame()),
               "Input should be a 'coda::mcmc' or a 'coda::mcmc.list' object.")
})