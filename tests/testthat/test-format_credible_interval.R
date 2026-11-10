library(testthat)

####  Definitions  ####

test_that("Credible intervals are formatted on vectors and frames", {
  # basic formatting
  format_credible_interval(1.23, 0.9, 1.5) |> expect_match("\\[")

  # as_integer
  format_credible_interval(2.1, 1.7, 2.3, as_integer = TRUE,
                           fmt = "%s-%s-%s") |>
    expect_identical("2-2-2")

  # NA mean returns NA_character_
  format_credible_interval(NA, NA, NA) |> expect_identical(NA_character_)

  # equal bounds
  format_credible_interval(2, 2, 2) |> expect_type("character")

  # two significant digits of the estimate, one of the width, none above units
  format_credible_interval(c(1.97, 9.9, 1234, 0.52013),
                           c(0.673, 2.6, 1180, 0.51904),
                           c(3.24, 23.7, 1290, 0.52097)) |>
    expect_identical(c("2.0 [0.7, 3.2]", "9.9 [2.6, 23.7]",
                       "1234 [1180, 1290]", "0.520 [0.519, 0.521]"))

  # interval covering zero -> two significant digits of the width
  format_credible_interval(c(0.0012, 0), c(-0.5, 0), c(0.5, 0.3)) |>
    expect_identical(c("0.0 [-0.5, 0.5]", "0.00 [0.00, 0.30]"))

  # inconsistent interval warns
  format_credible_interval(5, 10, 20) |>
    expect_warning("inconsistent")

  # mismatched lengths error
  format_credible_interval(c(1, 2), c(1), c(2)) |>
    expect_error("same length")

  # wrong str_format error
  format_credible_interval(1, 0, 2, fmt = "%s only one") |>
    expect_error("three.*placeholder")

  # data.frame input
  mean_frame <- data.frame(a = c(1, 2), b = c(3, 4))
  upper_frame <- data.frame(a = c(2, 3), b = c(4, 5))
  result <- format_credible_interval(
    mean_frame, data.frame(a = c(0, 1), b = c(2, 3)), upper_frame
  )
  result |> expect_s3_class("data.frame")
  result |> expect_dim(dim(mean_frame))

  # mismatched data.frame dimensions
  format_credible_interval(mean_frame, data.frame(a = c(0)), upper_frame) |>
    expect_error("same dimensions")
})
