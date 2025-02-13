library(testthat)

test_that("Able to extract results from coda", {
  
  expect_error(extract_results_from_coda(data.frame()), 
               "Input should be a coda 'mcmc' or 'mcmc.list' object.")
})