library(testthat)

test_that("Model generator inconsistent calls", {
  # blank household model object
  hhgen <- housepi:::HousepiGenerator$new("blank")

  # traces
  expect_error(hhgen$traces)

  # household traces
  expect_error(hhgen$traces, regexp = "No traces")
  expect_error(hhgen$house_loglike_traces, regexp = "No traces")

  # inference without definitions
  expect_error(hhgen$infer())
})
