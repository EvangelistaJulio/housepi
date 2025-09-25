library(testthat)

test_that("Default parameter definitions are correct", {
  param <- housepi:::HousepiParameter$new()
  
  expect_error(
    param$add(name = "dummy", value = 3, const = FALSE),
    regexp = "No prior"
  )
  
  # default sample framework
  param$load_list(dummy = list(value = 3, prior = "unif",
                               args = c(min = 0, max = 5)))
  
  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|S")
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  expect_equal(param$defs$rule, "sample")
  expect_equal(param$defs$scale, NA)
  expect_equal(param$defs$prior, "unif|min:max|{0|S}:{5|S}")
  
  # load_list and single add should be coherent
  tmp <- housepi:::HousepiParameter$new()
  tmp$add(name = "dummy", value = 3, prior = "unif",
          args = c(min = 0, max = 5))
  expect_equal(param$defs, tmp$defs)
  
  # given proposer rule framework
  param <- housepi:::HousepiParameter$new()
  param$add(name = "dummy", prior = "unif", args = c(min = 0, max = 5),
            rule = "norm")
  
  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, NA)
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  expect_equal(param$defs$rule, "norm")
  expect_equal(param$defs$scale, 1)
  expect_equal(param$defs$prior, "unif|min:max|{0|S}:{5|S}")
})

test_that("Constant parameter definitions are correct", {
  param <- housepi:::HousepiParameter$new()
  expect_error(param$add(name = "dummy"), regexp = "Missing value")
  expect_error(param$add(name = "dummy", prior = "unif",
                         args = c(min = 0, max = 5), const = TRUE))
  
  # load_list and single add should be coherent
  param <- param$load_list(dummy = 3)
  tmp <- housepi:::HousepiParameter$new()$add(name = "dummy", value = 3)
  expect_equal(param$defs, tmp$defs)
  
  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|S")
  expect_equal(param$defs$const, TRUE)
  expect_equal(param$defs$rule, NA)
  expect_equal(param$defs$scale, NA)
  
  expect_warning(param <- param$add(name = "dummy", value = 3, prior = "unif",
                                    args = c(min = 0, max = 5), const = TRUE),
                 regexp = "Overwriting")
  expect_equal(param$defs, tmp$defs)
})
