library(testthat)

test_that("Default infectivity profile definitions are correct", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Add a profile with default arguments
  infec$add(name = "dummy", discr = FALSE, args = list(a = 1, b = 2), max = 10)

  expect_equal(infec$defs$name, "dummy")
  expect_equal(infec$defs$discr, FALSE)
  expect_equal(infec$defs$arg_key, "a/b")
  expect_equal(infec$defs$arg_val, "1/2")
  expect_equal(infec$defs$arg_at, "S/S")
  expect_equal(infec$defs$max, "10")
  expect_equal(infec$defs$max_at, "S")
})

test_that("Infectivity profile with discrete values is correct", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Add a discrete profile
  infec$add(name = "discrete_profile", discr = TRUE, 
            args = list(x = 5), max = 20)

  expect_equal(infec$defs$name, "discrete_profile")
  expect_equal(infec$defs$discr, TRUE)
  expect_equal(infec$defs$arg_key, "x")
  expect_equal(infec$defs$arg_val, "5")
  expect_equal(infec$defs$arg_at, "S")
  expect_equal(infec$defs$max, "20")
  expect_equal(infec$defs$max_at, "S")
})

test_that("Error handling in infectivity profile definitions", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Missing name
  expect_error(infec$add(name = NA, discr = FALSE, args = list(a = 1)),
               regexp = "Incomplete definition: missing name")

  # Invalid function in arguments
  expect_error(
    infec$add(name = "invalid_args", args = list(a = quote(func()))),
    regexp = "not supported"
  )

  # Invalid function in max
  expect_error(infec$add(name = "invalid_max", max = quote(func())),
               regexp = "not supported")
})
