library(testthat)

test_that("Rlang standard functions work as expected", {
  expect_error(housepi::extract_func_(~log ( alpha @ B )),
               regexp = "Invalid type \"B\"")
  expect_equal(housepi::extract_func_(~log(alpha)),
               list(func = "log", key = "x", val = "alpha", at = "A"))
  expect_equal(housepi::extract_func_(~log(alpha @ P)),
               list(func = "log", key = "x", val = "alpha", at = "P"))
})

test_that("Rlang basic identity assignments work as expected", {
  expect_equal(housepi::extract_func_(+3),
               list(func = "..identity..", key = "x", val = "3", at = "S"))
  expect_equal(housepi::extract_func_(-1),
               list(func = "..identity..", key = "x", val = "-1", at = "S"))
  expect_equal(housepi::extract_func_(~ alpha),
               list(func = "..identity..", key = "x", val = "alpha", at = "A"))
  expect_equal(housepi::extract_func_(~ alpha @ P),
               list(func = "..identity..", key = "x", val = "alpha", at = "P"))
  expect_equal(housepi::extract_func_(~ +alpha),
               list(func = "..identity..", key = "x", val = "alpha", at = "A"))
  expect_equal(housepi::extract_func_(~ +alpha@P),
               list(func = "..identity..", key = "x", val = "alpha", at = "P"))
  expect_equal(housepi::extract_func_(~ -alpha),
               list(func = "minus", key = "x", val = "alpha", at = "A"))
  expect_equal(housepi::extract_func_(~ -alpha@D),
               list(func = "minus", key = "x", val = "alpha", at = "D"))
})

test_that("Rlang arithmetic operations work as expected", {
  expect_equal(housepi::extract_func_(~ alpha + beta),
               list(func = "add", key = "x/x2", val = "alpha/beta", at = "A/A"))
  expect_equal(housepi::extract_func_(~ alpha @ P+ beta),
               list(func = "add", key = "x/x2", val = "alpha/beta", at = "P/A"))
  expect_equal(housepi::extract_func_(~ alpha + beta @ P),
               list(func = "add", key = "x/x2", val = "alpha/beta", at = "A/P"))
  expect_equal(housepi::extract_func_(~ alpha @ P + beta @ P),
               list(func = "add", key = "x/x2", val = "alpha/beta", at = "P/P"))
  expect_equal(housepi::extract_func_(~ (alpha) -(beta)),
               list(func = "sub", key = "x/x2", val = "alpha/beta", at = "A/A"))
  expect_equal(housepi::extract_func_(~ alpha * (beta)),
               list(func = "mult", key = "x/x2", val = "alpha/beta", at = "A/A"))
  expect_equal(housepi::extract_func_(~ alpha * (beta @ P)),
               list(func = "mult", key = "x/x2", val = "alpha/beta", at = "A/P"))
  expect_equal(extract_func_(~ alpha @P * beta @ P),
               list(func = "mult", key = "x/x2", val = "alpha/beta", at = "P/P"))
})

test_that("Rlang nested expression raise errors", {
  # UPDATE: nested expression to be developed!
  expect_error(housepi::extract_func_(~a + 3 * b),
               regexp = "not supported")
  expect_error(housepi::extract_func_(~a + b - c),
               regexp = "not supported")
  expect_error(housepi::extract_func_(~log ( + alpha@P  )),
               regexp = "not supported")
})
