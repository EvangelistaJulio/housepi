library(testthat)

test_that("Augmented function definition is correct", {
  augm <- HousepiAugmented$new()

  # empty definition
  expect_error(augm$add())

  # e. g. infection date definition
  augm$add(name = "infec", value = " onset@A -   incub@A ")
  mask <- augm$defs$name == "infec"
  expect_equal(augm$defs[mask, "func"], "sub")
  expect_equal(augm$defs[mask, "farg_key"], "x/y")
  expect_equal(augm$defs[mask, "farg_val"], "onset/incub")
  expect_equal(augm$defs[mask, "farg_at"], "A/A")
  expect_equal(augm$defs[mask, "value"], NA)
})
