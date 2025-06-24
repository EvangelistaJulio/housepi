library(testthat)

test_that("Function exported from housepi are correct", {

  x <- c(0, 1, 5, 10, 50, 100, 500, 1000)

  expect_true(all(housepi:::housepi_function("pow", x) == x))
  expect_true(all(housepi:::housepi_function("pow", x, args = c(y = 2)) == x^2))
  expect_true(all(housepi:::housepi_function("pow", x, args = c(y = -3)) == x^(-3)))
  expect_true(all(housepi:::housepi_function("exp", x) == exp(x)))
  expect_true(all(housepi:::housepi_function("log", x) == log(x)))
  expect_true(all(housepi:::housepi_function("sin", x) == sin(x)))
  expect_true(all(housepi:::housepi_function("cos", x) == cos(x)))
  expect_true(all(housepi:::housepi_function("sin", x) == sin(x)))
  expect_true(all(housepi:::housepi_function("tan", x) == tan(x)))
})

test_that("String function extraction returns expected results", {
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
