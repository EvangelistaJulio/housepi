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

test_that("Deep clone helper work as expected on default types", {
  # scalar passthrough
  helper_clone(42) |> expect_equal(42)
  helper_clone("hello") |> expect_identical("hello")

  # data.frame independence
  original <- data.frame(x = 1:3)
  cloned <- helper_clone(original)
  cloned$x[1] <- 99
  original$x[1] |> expect_equal(1)
  cloned$x[1] |> expect_equal(99)

  # list of data.frames recursively cloned
  frames <- list(a = data.frame(x = 1), b = data.frame(y = 2))
  cloned <- helper_clone(frames)
  cloned$a$x[1] <- 99
  frames$a$x[1] |> expect_equal(1)
  cloned$a$x[1] |> expect_equal(99)
})

test_that("Draws are thinned evenly, first and last kept", {
  helper_thin_even_(seq(10), 4) |> expect_identical(c(1L, 4L, 7L, 10L))
  hh_data <- data.frame(datum = seq(10))
  helper_thin_even_(hh_data, 4)$datum |> expect_identical(c(1L, 4L, 7L, 10L))
  # fewer draws than asked, or no target: every draw kept
  helper_thin_even_(seq(3), 4) |> expect_identical(seq(3))
  helper_thin_even_(hh_data, NA) |> expect_identical(hh_data)
})
