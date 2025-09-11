library(testthat)

test_that("Augmented definition error are recognized", {
  augm <- housepi:::HousepiAugmented$new()
  expect_error(augm$add(name = "infec", args = c(min = 0)))
})

test_that("Constant value definition are correct", {
  augm <- housepi:::HousepiAugmented$new()$load_list(dummy = 3)
  tmp <- housepi:::HousepiAugmented$new()$add(name = "dummy", value = 3)
  expect_equal(augm$defs, tmp$defs)

  augm <- housepi:::HousepiAugmented$new()$load_list(dummy = ~3@D)
  tmp <- housepi:::HousepiAugmented$new()$add(name = "dummy", expr = ~3@D)
  testthat::expect_equal(augm$defs, tmp$defs)
})

test_that("Augmented variable definition is correct", {
  augm <- housepi:::HousepiAugmented$new()

  augm$add(name = "infec", value = ~expo@D, model = "unif",
           args = c(min = ~expo_min@A, max = ~expo_max@P),
           max = 3)

  expect_equal(augm$defs$name, "infec")
  expect_equal(augm$defs$status, as.integer(NA))
  expect_equal(augm$defs$value, "expo")
  expect_equal(augm$defs$value_at, "D")
  expect_equal(augm$defs$model, "unif")
  expect_equal(augm$defs$arg_key, "min/max")
  expect_equal(augm$defs$arg_val, "expo_min/expo_max")
  expect_equal(augm$defs$arg_at, "A/P")
  expect_equal(augm$defs$max, "3")
  expect_equal(augm$defs$max_at, "S")

  # consequent required definitions
  private <- augm$.__enclos_env__$private
  expect_equal(names(private$.find_def_required("P")), "augmented")
  expect_equal(unname(private$.find_def_required("P")), "expo_max")
  expect_equal(unname(private$.find_def_required("D")), "expo")
  expect_equal(unname(private$.find_def_required("A")), "expo_min")
})


test_that("Augmented sampling default stays coherent", {
  augm <- housepi:::HousepiAugmented$new()
  augm$add(name = "dummy", model = "norm")

  expect_equal(augm$defs$const, FALSE)
  expect_equal(augm$defs$init_rand, TRUE)
  expect_equal(augm$defs$rule, "sample")
})

test_that("Augmented function definition is correct", {
  augm <- housepi:::HousepiAugmented$new()

  # empty definition
  expect_error(augm$add(), regexp = "Incomplete definition")

  # difference definition
  augm$add(name = "dummy", expr = ~onset - incub)

  expect_equal(augm$defs$func, "sub")
  expect_equal(augm$defs$farg_key, "x/x2")
  expect_equal(augm$defs$farg_val, "onset/incub")
  expect_equal(augm$defs$farg_at, "A/A")
  expect_equal(augm$defs$value, NA)

  # overwriting warning
  expect_warning(
    augm$add(name = "dummy", func = "log", func_args = c(x = ~infec@A))
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

test_that("Status-specific definitions are safe", {
  augm <- housepi:::HousepiAugmented$new()

  required <- c(model = "log_proba",
                augmented = c("..proba", "..layer"))
  infectious_status_set <- c(4, 2, 3, 1)

  # no error
  expect_no_error(
    augm$load_list(
      ..proba = list(model = "norm", status = 1, min = 0,
                     args = c(mean = ~param@D)),
      ..proba = list(expr = ~..layer, status = c(2, 3, 4)),
      ..layer = list(model = "lnorm",
                     args = c(meanlog = ~param@D)),
      log_proba = ~log(..proba)
    )$assert(required, infectious_status_set)
  )

  # first generation of definitions
  err_msg <- as.character(expect_error(
    augm$load_list(
      ..proba = list(model = "norm", status = 1, min = 0,
                     args = c(mean = ~param@D)),
      ..proba = list(expr = ~..layer, status = 2),
      ..layer = list(model = "lnorm", status = 2,
                     args = c(meanlog = ~param@D)),
      #> NO ..proba provided for status 3 and 4
      log_proba = ~log(..proba)
    )$assert(required, infectious_status_set)
  ))
  expect_match(err_msg, "Missing definition")
  expect_match(err_msg, "status 3, 4")

  # second generation of definitions
  err_msg <- as.character(expect_error(
    augm$load_list(
      ..proba = list(model = "norm", status = 1, min = 0,
                     args = c(mean = ~param@D)),
      ..proba = list(expr = ~..layer, status = c(2, 3, 4)),
      ..layer = list(model = "lnorm", status = c(2, 3),
                     args = c(meanlog = ~param@D)),
      #> NO ..layer provided for status 3, 4
      log_proba = ~log(..proba)
    )$assert(required, infectious_status_set)
  ))
  expect_match(err_msg, "Missing definition")
  expect_match(err_msg, "status 4")
})
