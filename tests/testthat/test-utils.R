library(testthat)

test_that("String function extraction", {
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
