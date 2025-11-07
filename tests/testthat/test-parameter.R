library(testthat)

test_that("Default parameter definitions are correct", {
  param <- housepi:::HousepiParameter$new()

  expect_error(
    param$add(name = "dummy", value = 3, const = FALSE),
    regexp = "No prior"
  )

  # default sample framework
  param$load_list(dummy = list(value = 3, prior = "uniform",
                               min = 0, max = 5))

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|N")
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  expect_equal(param$defs$proposal, "sample")
  expect_equal(param$defs$proposal_scale, NA_real_)
  expect_equal(param$defs$prior, "uniform")
  expect_equal(param$defs$min, "0|N")
  expect_equal(param$defs$max, "5|N")

  # load_list and single add should be coherent
  tmp <- housepi:::HousepiParameter$new()
  tmp$add(name = "dummy", value = 3, prior = "uniform",
          min = 0, max = 5)
  expect_equal(param$defs, tmp$defs)

  # proposal
  param <- housepi:::HousepiParameter$new()
  param$add(name = "dummy", prior = "uniform", min = 0, max = 5,
            proposal = "normal")

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, NA_real_)
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  expect_equal(param$defs$proposal, "normal")
  expect_equal(param$defs$proposal_scale, 1)
  expect_equal(param$defs$prior, "uniform")
  expect_equal(param$defs$min, "0|N")
  expect_equal(param$defs$max, "5|N")
})

test_that("Constant parameter definitions are correct", {
  param <- housepi:::HousepiParameter$new()
  expect_error(param$add(name = "dummy"), regexp = "Incomplete definition")
  expect_error(param$add(name = "dummy", prior = "uniform",
                         min = 0, max = 5, const = TRUE))

  # load_list and single add should be coherent
  param <- param$load_list(dummy = 3)
  tmp <- housepi:::HousepiParameter$new()$add(name = "dummy", value = 3)
  expect_equal(param$defs, tmp$defs)

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|N")
  expect_equal(param$defs$const, TRUE)
  expect_equal(param$defs$proposal, NA_character_)
  expect_equal(param$defs$proposal_scale, NA_real_)

  expect_warning(
    param <- param$add(name = "dummy", value = 3, prior = "uniform",
                       min = 0, max = 5, const = TRUE),
    regexp = "Overwriting"
  )
  expect_equal(param$defs, tmp$defs)
})
