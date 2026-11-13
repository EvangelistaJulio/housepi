library(testthat)

####  Definitions  ####

test_that("Overwrite is a no-op on the same path and raises on failure", {
  # no-op on identical paths
  path <- tempfile()
  file.create(path)
  overwrite(path, path) |> expect_invisible()

  # destination directory does not exist -> file.copy returns FALSE
  overwrite(path, file.path(tempfile(), "missing", "x")) |>
    expect_error("user permissions")
})

test_that("Draws are thinned evenly, first and last kept", {
  helper_thin_even_(seq(10), 4) |> expect_identical(c(1L, 4L, 7L, 10L))
  hh_data <- data.frame(datum = seq(10))
  helper_thin_even_(hh_data, 4)$datum |> expect_identical(c(1L, 4L, 7L, 10L))
  # fewer draws than asked, or no target: every draw kept
  helper_thin_even_(seq(3), 4) |> expect_identical(seq(3))
  helper_thin_even_(hh_data, NA) |> expect_identical(hh_data)
})
