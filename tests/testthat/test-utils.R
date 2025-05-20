library(testthat)

test_that("String function extraction", {
  s <- "=alpha@P"
  extract_augm_func_(s)

  s <- "=alpha@P+"
  expect_error(extract_augm_func_(s))

  s <- "alpha@P+beta@P"
  expect_equal(extract_augm_func_(s)$var, type_join_(c("alpha", "beta")))

  s <- "alpha@P+beta@P+"
  expect_error(extract_augm_func_(s))
})
