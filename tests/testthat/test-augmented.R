library(testthat)

aux_add <- function(self, name = NA_character_,  value = NULL,
                    status = NA_integer_, const = NA, discr = FALSE,
                    func = NULL, func_args = NULL,
                    model = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_,
                    proposal_scale = NA_real_) {
  self$add_quo(name = name, value = rlang::enquo(value),
               status = status, const = const, discr = discr,
               func = rlang::enquo(func),
               func_args = rlang::enquos(func_args),
               model = model, rlang::enquos(...),
               min = rlang::enquo(min), max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale)
  invisible(self)
}

test_that("Augmented definition error are recognized", {
  augm <- HousepiAugmented$new()
  expect_error(aux_add(augm, name = "infec@"),
               regexp = "invalid characters")
  expect_error(aux_add(augm, name = "infec"),
               regexp = "Incomplete definition")
})

test_that("Constant value definition are coherent", {
  augm <- helper_test_aux_load(HousepiAugmented$new(), dummy = 3)
  tmp <- aux_add(HousepiAugmented$new(), name = "dummy", value = 3)
  expect_equal(augm$defs, tmp$defs)

  augm <- helper_test_aux_load(HousepiAugmented$new(), dummy = 3@D)
  tmp <- aux_add(HousepiAugmented$new(), name = "dummy", value = 3@D)
  expect_equal(augm$defs, tmp$defs)
})

test_that("Undefined definitions are controoled", {
  augm <- HousepiAugmented$new()
  expect_error(aux_add(augm, name = "test", proposal = "sample", value = NA),
               regexp = "No model")
  expect_error(aux_add(augm, name = "test", value = NA),
               regexp = "Missing value")

  augm <- HousepiAugmented$new()
  aux_add(augm, name = "test", value = 3, const = TRUE, proposal = "sample")
  expect_true(is.na(augm$defs$proposal))
})

test_that("Augmented variable definitions are correct", {
  # basic constant definition
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", value = 3)
  expect_equal(augm$defs$name, "dummy")
  expect_equal(augm$defs$status, NA_integer_)
  expect_equal(augm$defs$const, TRUE)
  expect_equal(augm$defs$discr, FALSE)
  expect_equal(augm$defs$value, "3|N")
  expect_equal(augm$defs$proposal, NA_character_)
  expect_equal(augm$defs$proposal_scale, NA_real_)

  # variable definition
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", value = expo@D, model = "uniform",
          min = expo_min@A, max = expo_max@P)
  expect_equal(augm$defs$name, "dummy")
  expect_equal(augm$defs$status, NA_integer_)
  expect_equal(augm$defs$value, "expo|D")
  expect_equal(augm$defs$model, "uniform")
  expect_equal(augm$defs$min, "expo_min|A")
  expect_equal(augm$defs$max, "expo_max|P")

  # consequent required definitions
  expect_equal(names(augm$param_required), "augmented")
  expect_equal(unname(augm$param_required), "expo_max")
  expect_equal(unname(augm$data_required), "expo")
  expect_equal(unname(augm$augm_required), "expo_min")
})

test_that("Augmented sampling default is correct", {
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", model = "norm")

  expect_equal(augm$defs$const, FALSE)
  expect_equal(augm$defs$proposal, "sample")
})

test_that("Augmented function definition is correct", {
  augm <- HousepiAugmented$new()
  expect_error(aux_add(augm), regexp = "non-empty")

  # difference definition
  aux_add(augm, name = "dummy", func = onset - incub)
  expect_equal(augm$defs$func, "sub|x:x2|{onset|A}:{incub|A}")
  expect_equal(augm$defs$value, NA_real_)

  # overwriting warnings
  expect_warning(
    aux_add(augm, name = "dummy", func = "log", func_args = c(x = infec@A)),
    regexp = "Overwriting"
  )
  aux_add(augm, name = "rate", value = 3)
  expect_warning(aux_add(augm, name = "dummy", func = exp(-rate)),
                 regexp = "Overwriting")
})

test_that("Augmented load enable special definitions", {
  # alternative function call definitions
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", func = "log", func_args = c(x = infec@A))

  tmp <- HousepiAugmented$new()
  helper_test_aux_load(tmp, dummy = log(infec@A))
  expect_equal(tmp$defs, augm$defs)

  expect_warning(aux_add(tmp, name = "dummy", func = log(infec@A)),
                 regexp = "Overwriting")
  expect_equal(tmp$defs, augm$defs)

  # multi-status definition
  tmp$defs <- NULL
  helper_test_aux_load(tmp, dummy = list(func = log(infec@A), status = c(1, 3)))
  tmp$defs$status <- NULL

  expect_equal(tmp$defs[1, , drop = TRUE], tmp$defs[2, , drop = TRUE])
})

test_that("Augmented support rlang double bang", {
  augm <- HousepiAugmented$new()
  external_value <- 33
  str_value <- sprintf("%s|N", external_value)
  helper_test_aux_load(
    augm,
    A = !!external_value,
    B = log(infec@A + !!external_value),
    C = c(model = "uniform",
          min = !!external_value, max = log(infec@A + !!external_value)),
    D = list(model = "uniform",
             min = !!external_value, max = log(infec@A + !!external_value))
  )
  expect_equal(subset(augm$defs, name == "A")$value, str_value)
  expect_equal(subset(augm$defs, name == "B")$func,
               sprintf("log|x|{sum|x:x2|{infec|A}:{%s}}", str_value))
  subset(augm$defs, name == "C")
})

test_that("Symbolic definitions are equivalent", {
  augm <- HousepiAugmented$new()
  val <- 1
  helper_test_aux_load(
    augm,
    var1 = list(model = "normal", mean = indiv_shape + 1),
    var2 = list(model = "normal", mean = indiv_shape + !!val),
    var3 = list(model = "normal", mean = ~ indiv_shape + 1),
    var4 = list(model = "normal", mean = ~ indiv_shape + !!val)
  )
  aux <- augm$defs
  aux$name <- NULL
  # check that all rows are identical
  expect_equal(nrow(unique(aux)), 1)

  # check injection from external
  external <- list(
    some_func = list(func = "gamma",
                     func_args = c(shape = ~ shape, scale = ~ scale))
  )
  expect_error(helper_test_aux_load(augm, external),
               regexp = "rlang injection")

  helper_test_aux_load(augm, !!!external)
  expect_equal(augm$defs$func, "gamma|shape:scale|{shape|A}:{scale|A}")
})
