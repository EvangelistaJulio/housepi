library(testthat)

test_that("IndividualSampler initializes with levels encoding", {
  sampler <- housepi:::IndividualSampler$new(
    sex = list(levels = c("M", "F"))
  )
  expect_s3_class(sampler, "IndividualSampler")
})

test_that("IndividualSampler initializes with breaks encoding", {
  sampler <- housepi:::IndividualSampler$new(
    age = list(breaks = c(18, 65))
  )
  expect_s3_class(sampler, "IndividualSampler")
})

test_that("IndividualSampler errors on invalid input", {
  # no arguments
  expect_error(housepi:::IndividualSampler$new(),
               regexp = "No definition provided")
  # non-list argument
  expect_error(housepi:::IndividualSampler$new(sex = c("M", "F")),
               regexp = "only lists")
  # missing breaks/levels
  expect_error(housepi:::IndividualSampler$new(sex = list(foo = "bar")),
               regexp = "Incomplete encoding")
})

test_that("IndividualSampler encode/decode round-trips correctly", {
  sampler <- housepi:::IndividualSampler$new(
    sex = list(levels = c("M", "F"))
  )

  # create test data with levels only
  input <- data.frame(sex = c("M", "F", "M"), stringsAsFactors = FALSE)

  # round-trip
  encoded <- sampler$encode(input)
  expect_true(is.numeric(encoded))
  expect_equal(length(encoded), 3)

  decoded <- sampler$decode(encoded)
  re_encoded <- sampler$encode(decoded)
  expect_equal(encoded, re_encoded)
  expect_equal(decoded$sex, input$sex)
})

test_that("IndividualSampler multi-feature encode/decode works", {
  sampler <- housepi:::IndividualSampler$new(
    sex = list(levels = c("M", "F")),
    vaccinated = list(levels = c("yes", "no"))
  )

  input <- data.frame(
    sex = c("M", "F", "M", "F"),
    vaccinated = c("yes", "no", "no", "yes"),
    stringsAsFactors = FALSE
  )

  encoded <- sampler$encode(input)
  decoded <- sampler$decode(encoded)
  re_encoded <- sampler$encode(decoded)
  expect_equal(encoded, re_encoded)
})

test_that("HousepiSampler initializes and sets ascending", {
  sampler <- housepi:::HousepiSampler$new(
    sex = list(levels = c("M", "F")),
    ascending = TRUE
  )
  expect_s3_class(sampler, "HousepiSampler")
  expect_s3_class(sampler, "IndividualSampler")
})

test_that("HousepiSampler hierarchy works", {
  sampler <- new_household_sampler(sex = list(levels = c("M", "F")))
  # hierarchy uses expand.grid on the levels themselves (no "unknown" row)
  h <- sampler$hierarchy
  expect_true(is.data.frame(h))
  expect_true("sex" %in% colnames(h))
  expect_true(".features" %in% colnames(h))

  # read only
  expect_error(sampler$hierarchy <- "foo", regexp = "read only")
})

test_that("HousepiSampler compute_household_marginals works", {
  sampler <- new_household_sampler(sex = list(levels = c("M", "F")))
  data_wide <- data.frame(
    id = paste0("p", 1:6),
    house_id = c("h1", "h1", "h2", "h2", "h2", "h3"),
    status = c(1, 0, 1, 0, 0, 1),
    sex = c("M", "F", "M", "M", "F", "F")
  )

  expect_no_error(sampler$compute_household_marginals(data_wide))

  # missing house_id column errors
  bad_data <- data.frame(id = "p1", sex = "M")
  expect_error(sampler$compute_household_marginals(bad_data),
               regexp = "house_id")
})

test_that("new_household_sampler exported constructor works", {
  sampler <- new_household_sampler(sex = list(levels = c("M", "F")))
  expect_s3_class(sampler, "HousepiSampler")
})

test_that("compute_household_marginals exported function clones", {
  sampler <- new_household_sampler(sex = list(levels = c("M", "F")))
  data_wide <- data.frame(
    id = paste0("p", 1:4),
    house_id = c("h1", "h1", "h2", "h2"),
    status = c(1, 0, 1, 0),
    sex = c("M", "F", "M", "F")
  )

  result <- compute_household_marginals(sampler, data_wide)
  expect_s3_class(result, "HousepiSampler")
  # returned object is a different instance (deep clone)
  expect_false(identical(sampler, result))
})

test_that("sample_households exported function validates input", {
  expect_error(sample_households("not a sampler", 10))
})

test_that("HousepiSampler sample_households size validation", {
  sampler <- housepi:::HousepiSampler$new(
    sex = list(levels = c("M", "F"))
  )
  expect_error(sampler$sample_households(size = 0), regexp = "greater than 0")
})
