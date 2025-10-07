library(testthat)

test_that("Augmented definition error are recognized", {
  augm <- housepi:::HousepiAugmented$new()
  expect_error(augm$add(name = "infec", args = c(min = 0)))
})

test_that("Constant value definition are coherent", {
  augm <- housepi:::HousepiAugmented$new()$load_list(dummy = 3)
  tmp <- housepi:::HousepiAugmented$new()$add(name = "dummy", value = 3)
  expect_equal(augm$defs, tmp$defs)

  augm <- housepi:::HousepiAugmented$new()$load_list(dummy = ~3@D)
  tmp <- housepi:::HousepiAugmented$new()$add(name = "dummy", expr = ~3@D)
  expect_equal(augm$defs, tmp$defs)
})

test_that("Undefined definitions are controoled", {
  augm <- housepi:::HousepiAugmented$new()
  expect_error(augm$add(name = "test", rule = "sample", value = NA),
               regexp = "No model")

  expect_error(augm$add(name = "test", rule = "norm", value = NA),
               regexp = "undefined")

  expect_no_error(augm$add(name = "test", value = NA))

  augm <- housepi:::HousepiAugmented$new()
  augm$add(name = "test", value = 3, const = TRUE, rule = "sample")
  expect_true(is.na(augm$defs$rule))
})

test_that("Augmented variable definitions are correct", {
  # basic constant definition
  augm <- housepi:::HousepiAugmented$new()
  augm$add(name = "dummy", value = 3)
  expect_equal(augm$defs$name, "dummy")
  expect_equal(augm$defs$status, NA_integer_)
  expect_equal(augm$defs$const, TRUE)
  expect_equal(augm$defs$discr, FALSE)
  expect_equal(augm$defs$value, "3|N")
  expect_equal(augm$defs$rule, NA)
  expect_equal(augm$defs$scale, NA)
  
  # variable definition
  augm <- housepi:::HousepiAugmented$new()
  augm$add(name = "dummy", value = ~expo@D, model = "uniform",
           args = c(min = ~expo_min@A, max = ~expo_max@P),
           max = 3)
  expect_equal(augm$defs$name, "dummy")
  expect_equal(augm$defs$status, as.integer(NA))
  expect_equal(augm$defs$value, "expo|D")
  expect_equal(augm$defs$model, "uniform|min:max|{expo_min|A}:{expo_max|P}")
  expect_equal(augm$defs$max, "3|N")

  # consequent required definitions
  private <- augm$.__enclos_env__$private
  expect_equal(names(private$.find_def_required("P")), "augmented")
  expect_equal(unname(private$.find_def_required("P")), "expo_max")
  expect_equal(unname(private$.find_def_required("D")), "expo")
  expect_equal(unname(private$.find_def_required("A")), "expo_min")
})

test_that("Augmented sampling default is correct", {
  augm <- housepi:::HousepiAugmented$new()
  augm$add(name = "dummy", model = "norm")

  expect_equal(augm$defs$const, FALSE)
  expect_equal(augm$defs$rule, "sample")
})

test_that("Augmented function definition is correct", {
  augm <- housepi:::HousepiAugmented$new()
  expect_error(augm$add(), regexp = "Incomplete definition")

  # difference definition
  augm$add(name = "dummy", expr = ~onset - incub)
  expect_equal(augm$defs$func, "sub|x:x2|{onset|A}:{incub|A}")
  expect_equal(augm$defs$value, NA)

  # overwriting warnings
  expect_warning(
    augm$add(name = "dummy", func = "log", func_args = c(x = ~infec@A)),
    "Overwriting existing"
  )
  augm$add(name = "rate", value = 3)
  expect_warning(
    augm$add(name = "dummy", expr = ~exp(-rate)),
    "Overwriting existing"
  )
})

test_that("Augmented load enable special definitions", {
  # alternative function call definitions
  augm <- housepi:::HousepiAugmented$new()
  augm$add(name = "dummy", func = "log", func_args = c(x = ~infec@A))

  tmp <- housepi:::HousepiAugmented$new()
  tmp$load_list(dummy = ~log(infec@A))
  expect_equal(tmp$defs, augm$defs)

  expect_warning(tmp$add(name = "dummy", expr = ~log(infec@A)),
                 regexp = "Overwriting existing")
  expect_equal(tmp$defs, augm$defs)

  # multi-status definition
  tmp$defs <- NULL
  tmp$load_list(dummy = list(expr = ~log(infec@A), status = c(1, 3)))
  tmp$defs$status <- NULL

  expect_equal(tmp$defs[1, , drop = TRUE], tmp$defs[2, , drop = TRUE])
})
