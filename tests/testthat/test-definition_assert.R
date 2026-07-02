library(testthat)

test_that("Helper format errors keeps full multi-line text of the pick", {
  # single error: returns the full (multi-line) text for the cli split
  expect_equal(helper_format_error(list("first line\nlast line")),
               "first line\nlast line")

  # multiple errors: "Unused" or "Missing" are prioritized
  err <- list(a = "some error", b = "Unused argument x")
  expect_equal(helper_format_error(err), "Unused argument x")
  err <- list(a = "some error", b = "Missing argument x")
  expect_equal(helper_format_error(err), "Missing argument x")

  # "func" value in err triggers func branch (checks value, not key)
  err <- list(dist = "dist error", func = "func")
  expect_equal(helper_format_error(err), list(func = "func"))

  # fallback to first
  err <- list(a = "first", b = "second")
  expect_equal(helper_format_error(err), "first")
})

test_that("Assert valid definition name rejects invalid names", {
  # empty string
  expect_error(assert_valid_definition_name("", "parameter"),
               regexp = "non-empty")
  # NA
  expect_error(assert_valid_definition_name(NA, "parameter"),
               regexp = "non-empty")
  # non-character
  expect_error(assert_valid_definition_name(42, "parameter"),
               regexp = "non-empty")
  # special characters
  expect_error(assert_valid_definition_name("a b", "parameter"),
               regexp = "invalid characters")
  expect_error(assert_valid_definition_name("a@b", "parameter"),
               regexp = "invalid characters")
  # reserved name pi
  expect_error(assert_valid_definition_name("pi", "parameter"),
               regexp = "reserved")
  # valid names pass
  expect_no_error(assert_valid_definition_name("alpha", "parameter"))
  expect_no_error(assert_valid_definition_name("beta_1", "parameter"))
  expect_no_error(assert_valid_definition_name("ga.2", "parameter"))
})

test_that("Assert valid status rejects non-positive integers", {
  # valid
  expect_no_error(assert_valid_status(c(NA, 1, 2, 3)))
  # zero
  expect_error(assert_valid_status(c(0, 1)),
               regexp = "positive integer")
  # negative
  expect_error(assert_valid_status(c(-1, 2)),
               regexp = "positive integer")
  # decimal
  expect_error(assert_valid_status(c(1.5)),
               regexp = "positive integer")
})
