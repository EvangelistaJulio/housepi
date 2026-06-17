library(testthat)

aux_add <- function(self, name = NA_character_,  value = NULL,
                    status = NA_integer_, const = NA, discr = FALSE,
                    func = NULL, func_args = NULL,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_,
                    proposal_scale = NA_real_) {
  self$add_quo(name = name, value = rlang::enquo(value),
               status = status, const = const, discr = discr,
               func = rlang::enquo(func),
               func_args = rlang::enquos(func_args),
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min), max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale)
  invisible(self)
}

test_that("Augmented definition error are recognized", {
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "infec@") |>
    expect_error("invalid characters")
  aux_add(augm, name = "infec") |>
    expect_error("Incomplete definition")
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
  aux_add(augm, name = "test", proposal = "sample", value = NA) |>
    expect_error("No prior")
  aux_add(augm, name = "test", value = NA) |>
    expect_error("Missing value")

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
  aux_add(augm, name = "dummy", value = expo@D, prior = "uniform",
          min = expo_min@A, max = expo_max@P)
  expect_equal(augm$defs$name, "dummy")
  expect_equal(augm$defs$status, NA_integer_)
  expect_equal(augm$defs$value, "expo|D")
  expect_equal(augm$defs$prior, "uniform")
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
  aux_add(augm, name = "dummy", prior = "norm")

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
  aux_add(augm, name = "dummy", func = "log", func_args = c(x = infec@A)) |>
    expect_warning("Overwriting")
  aux_add(augm, name = "rate", value = 3)
  aux_add(augm, name = "dummy", func = exp(-rate)) |>
    expect_warning("Overwriting")
})

test_that("Augmented load enable special definitions", {
  # alternative function call definitions
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", func = "log", func_args = c(x = infec@A))

  tmp <- HousepiAugmented$new()
  helper_test_aux_load(tmp, dummy = log(infec@A))
  expect_equal(tmp$defs, augm$defs)

  aux_add(tmp, name = "dummy", func = log(infec@A)) |>
    expect_warning("Overwriting")
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
    C = c(prior = "uniform",
          min = !!external_value, max = log(infec@A + !!external_value)),
    D = list(prior = "uniform",
             min = !!external_value, max = log(infec@A + !!external_value))
  )
  expect_equal(subset(augm$defs, name == "A")$value, str_value)
  expect_equal(subset(augm$defs, name == "B")$func,
               sprintf("log|x|{sum|x:x2|{infec|A}:{%s}}", str_value))
  subset(augm$defs, name == "C")
})

###############################
## AUGMENTED : ACTIVE FIELDS ##
###############################

test_that("Active fields are read-only", {
  HousepiAugmented$new() |>
    helper_test_read_only(
      "augm_required", "data_required",
      "group_rates_required", "param_required"
    )
})

test_that("Symbolic definitions are equivalent", {
  augm <- HousepiAugmented$new()
  val <- 1
  helper_test_aux_load(
    augm,
    var1 = list(prior = "normal", mean = indiv_shape + 1),
    var2 = list(prior = "normal", mean = indiv_shape + !!val),
    var3 = list(prior = "normal", mean = ~ indiv_shape + 1),
    var4 = list(prior = "normal", mean = ~ indiv_shape + !!val)
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
  helper_test_aux_load(augm, external) |>
    expect_error("rlang injection")

  helper_test_aux_load(augm, !!!external)
  expect_equal(augm$defs$func, "gamma|shape:scale|{shape|A}:{scale|A}")
})

test_that("Observation node definitions are encoded", {
  augm <- HousepiAugmented$new()
  helper_test_aux_load(
    augm,
    onset_date = list(obs = "normal", mean = onset_time@A, sd = 1),
    .kind = "obs"
  )
  d <- augm$defs
  # marked as observation, data-bound and constant (never sampled)
  expect_identical(d[["kind"]], "obs")
  expect_true(d[["const"]])
  expect_identical(d[["value"]], "onset_date|D")
  expect_true(is.na(d[["proposal"]]))
  # distribution string carries the observation arguments
  expect_identical(d[["prior"]], "normal|mean:sd|{onset_time|A}:{1|N}")
  # the data column is required, the latent is an augmented dependency
  expect_identical(unname(augm$data_required), "onset_date")
  expect_identical(unname(augm$augm_required), "onset_time")
  # reported as a generated data column (excluded from sim requirement)
  expect_identical(augm$obs_data_names, "onset_date")
})

test_that("Observation node rejects conflicting options", {
  augm <- HousepiAugmented$new()
  # prior and obs are mutually exclusive
  helper_test_aux_load(augm, x = list(prior = "normal", obs = "normal"),
                       .kind = "obs") |>
    expect_error("only one between 'prior' and 'obs'")
  # observation nodes are data-bound, never sampled
  helper_test_aux_load(augm, x = list(obs = "normal", proposal = "sample"),
                       .kind = "obs") |>
    expect_error("data-bound")
})

test_that("Deterministic observation node is encoded as a function obs", {
  augm <- HousepiAugmented$new()
  helper_test_aux_load(
    augm,
    incub_period = list(prior = "gamma", shape = 2, scale = 1),
    infec_time = -incub_period
  )
  # a bare expression (no 'obs =' distribution) is a deterministic observation
  helper_test_aux_load(augm, onset_date = floor(infec_time + incub_period),
                       .kind = "obs")
  d <- augm$defs[augm$defs[["name"]] == "onset_date", ]
  # tagged as observation, carrying a deterministic function and no distribution
  expect_identical(d[["kind"]], "obs")
  expect_true(is.na(d[["prior"]]))
  expect_true(is.na(d[["value"]]))
  expect_match(d[["func"]], "floor")
  # still reported as a generated data column
  expect_identical(augm$obs_data_names, "onset_date")
})

test_that("Observation and augmented loads each own only their definitions", {
  augm <- HousepiAugmented$new()
  helper_test_aux_load(
    augm,
    incub_period = list(prior = "gamma", shape = 2, scale = 1),
    infec_time = -incub_period
  )
  # a fresh (non-append) observation load keeps the augmented variables it needs
  helper_test_aux_load(
    augm,
    onset_date = list(obs = "normal", mean = infec_time@A, sd = 1),
    .kind = "obs"
  )
  expect_setequal(augm$defs[["name"]],
                  c("incub_period", "infec_time", "onset_date"))

  # symmetrically, a fresh augmented load keeps the observation node
  helper_test_aux_load(
    augm, incub_period = list(prior = "gamma", shape = 2, scale = 1)
  )
  expect_true("onset_date" %in% augm$defs[["name"]])
  expect_identical(augm$obs_data_names, "onset_date")
})

test_that("Household observations coexist with augmented and observation defs", {
  augm <- HousepiAugmented$new()
  helper_test_aux_load(
    augm, incub_period = list(prior = "gamma", shape = 2, scale = 1)
  )
  helper_test_aux_load(
    augm,
    onset_date = list(obs = "normal", mean = incub_period@A, sd = 1),
    .kind = "obs"
  )
  # a fresh household-observation load owns only the "house-obs" rows, leaving the
  # augmented variable and the observation node it depends on untouched
  helper_test_aux_load(augm, first_test = min(test_date@D), .kind = "house-obs")
  expect_setequal(augm$defs[["name"]],
                  c("incub_period", "onset_date", "first_test"))
  expect_identical(augm$defs[["kind"]][augm$defs[["name"]] == "first_test"],
                   "house-obs")
  # the reduction over a data source is encoded as a variadic call (a call lands
  # in 'func', a bare terminal would land in 'value')
  expect_identical(augm$defs[["func"]][augm$defs[["name"]] == "first_test"],
                   "min|x|{test_date|D}")

  # a fresh (non-append) household load drops only the house rows
  helper_test_aux_load(augm, last_test = max(test_date@D), .kind = "house-obs")
  expect_false("first_test" %in% augm$defs[["name"]])
  expect_setequal(augm$defs[["name"]],
                  c("incub_period", "onset_date", "last_test"))
})
