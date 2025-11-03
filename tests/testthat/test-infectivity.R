library(testthat)

test_that("Default infectivity profile definitions are correct", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Add a profile with default arguments
  infec$add(name = "dummy", discr = FALSE, args = c(a = 1, b = 2), max = 10)

  expect_equal(infec$defs$name, "dummy|a:b|{1|N}:{2|N}")
  expect_equal(infec$defs$discr, FALSE)
  expect_equal(infec$defs$max, "10|N")
})

test_that("Infectivity profile with discrete values is correct", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Add a discrete profile
  infec$add(name = "discrete_profile", discr = TRUE,
            args = list(x = 5), max = 20)

  expect_equal(infec$defs$name, "discrete_profile|x|{5|N}")
  expect_equal(infec$defs$discr, TRUE)
  expect_equal(infec$defs$max, "20|N")
})

test_that("Error handling in infectivity profile definitions", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Missing name
  expect_error(infec$add(name = NA, discr = FALSE, args = list(a = 1)),
               regexp = "Incomplete definition: missing name")
})
