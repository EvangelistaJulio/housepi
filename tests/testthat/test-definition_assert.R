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
  assert_valid_definition_name("", "parameter") |>
    expect_error("non-empty")
  # NA
  assert_valid_definition_name(NA, "parameter") |>
    expect_error("non-empty")
  # non-character
  assert_valid_definition_name(42, "parameter") |>
    expect_error("non-empty")
  # special characters
  assert_valid_definition_name("a b", "parameter") |>
    expect_error("invalid characters")
  assert_valid_definition_name("a@b", "parameter") |>
    expect_error("invalid characters")
  # reserved name pi
  assert_valid_definition_name("pi", "parameter") |>
    expect_error("reserved")
  # a leading dot marks a hidden definition, forbidden by default
  assert_valid_definition_name(".ration", "observation") |>
    expect_error("cannot start with a dot")
  # valid names pass
  expect_no_error(assert_valid_definition_name("alpha", "parameter"))
  expect_no_error(assert_valid_definition_name("beta_1", "parameter"))
})

test_that("Assert valid status rejects non-positive integers", {
  # valid
  expect_no_error(assert_valid_status(c(NA, 1, 2, 3)))
  # zero
  assert_valid_status(c(0, 1)) |>
    expect_error("positive integer")
  # negative
  assert_valid_status(c(-1, 2)) |>
    expect_error("positive integer")
  # decimal
  assert_valid_status(c(1.5)) |>
    expect_error("positive integer")
})
