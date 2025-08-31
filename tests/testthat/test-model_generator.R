library(testthat)

test_that("Model generator inconsistent calls", {
  # blank household model object
  hhgen <- housepi:::HousepiGenerator$new("blank")

  # helper with wrong template name specification
  expect_error(housepi:::HousepiGenerator$new("augment example"),
               regexp = "augmented example")
  
  # inference without definitions
  expect_error(hhgen$infer())
})
