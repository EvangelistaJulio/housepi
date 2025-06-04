library(testthat)

test_that("Augmented definition error are recognized", {
  augm <- housepi:::HousepiAugmented$new()
  expect_error(augm$add(name = "infec", args = c(min = 0)))
})

test_that("Augmented variable definition is correct", {
  augm <- housepi:::HousepiAugmented$new()

  augm$add(name = "infec", value = "expo@D", model = "unif",
           args = c(min = "expo_min@A", max = "expo_max@P"),
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
  expect_equal(private$.find_def_required("P"), "expo_max")
  expect_equal(private$.find_def_required("D"), "expo")
  expect_equal(private$.find_def_required("A"), "expo_min")
})

test_that("Augmented function definition is correct", {
  augm <- housepi:::HousepiAugmented$new()

  # empty definition
  expect_error(augm$add())

  # difference definition
  augm$add(name = "dummy", value = " onset -   incub ")
  expect_equal(augm$defs$func, "sub")
  expect_equal(augm$defs$farg_key, "x/y")
  expect_equal(augm$defs$farg_val, "onset/incub")
  expect_equal(augm$defs$farg_at, "A/A")
  expect_equal(augm$defs$value, NA)

  # overwriting warning
  expect_warning(
    augm$add(name = "dummy", func = "log", func_args = c(x = "infec@A"))
  )

  # alternative function call definitions
  tmp <- housepi:::HousepiAugmented$new()
  tmp$load_list(dummy = "log(infec@A)")
  expect_true(all(tmp$defs == augm$defs, na.rm = TRUE))

  tmp <- housepi:::HousepiAugmented$new()
  tmp$add(name = "dummy", func = "log(infec@A)")
  expect_true(all(tmp$defs == augm$defs, na.rm = TRUE))
})
