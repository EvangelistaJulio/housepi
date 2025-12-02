library(testthat)

add_add <- function(self, name, status = NA_integer_,
                    discr = FALSE, ..., max = NULL) {
  self$add_quo(name = name, status = status, discr = discr,
               rlang::enquos(...), max = rlang::enquo(max))
  invisible(self)
}

test_that("Default infectivity profile definitions are correct", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Add a profile with default arguments
  add_add(infec, name = "dummy", discr = FALSE, a = 1, b = 2, max = 10)

  expect_equal(infec$defs$name, "dummy|a:b|{1|N}:{2|N}")
  expect_equal(infec$defs$status, NA_integer_)
  expect_equal(infec$defs$discr, FALSE)
  expect_equal(infec$defs$max, "10|N")
})

test_that("Infectivity profile with discrete values is correct", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Add a discrete profile
  add_add(infec, name = "discrete_profile", discr = TRUE,
          x = 5, max = 20)

  expect_equal(infec$defs$name, "discrete_profile|x|{5|N}")
  expect_equal(infec$defs$discr, TRUE)
  expect_equal(infec$defs$max, "20|N")
})

test_that("Error handling in infectivity profile definitions", {
  infec <- housepi:::HousepiInfecProfile$new()

  # Missing name
  expect_error(add_add(infec, name = NA, discr = FALSE, a = 1),
               regexp = "Incomplete definition: missing name")
})
