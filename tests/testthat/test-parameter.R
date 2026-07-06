library(testthat)

aux_add <- function(self, name = NA_character_, value = NULL,
                    const = NA, discr = FALSE,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_,
                    proposal_scale = NA_real_) {
  self$add_quo(name = name, value = rlang::enquo(value),
               const = const, discr = discr,
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min), max = rlang::enquo(max),
               proposal = proposal,
               proposal_scale = proposal_scale)
  invisible(self)
}

test_that("Default parameter definitions are correct", {
  param <- HousepiParameter$new()

  aux_add(param, name = "dummy", value = 3, const = FALSE) |>
    expect_error("prior")

  # active fields are read-only
  helper_test_read_only(param, "param_required")

  # default sample framework
  helper_test_aux_load(
    param, dummy = list(value = 3, prior = "uniform", min = 0, max = 5)
  )

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
  tmp <- HousepiParameter$new()
  aux_add(tmp, name = "dummy", value = 3, prior = "uniform",
          min = 0, max = 5)
  expect_equal(param$defs, tmp$defs)

  # proposal
  param <- HousepiParameter$new()
  aux_add(param, name = "dummy", prior = "uniform", min = 0, max = 5,
          proposal = "normal")

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, NA_real_)
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  expect_equal(param$defs$proposal, "normal")
  # NA_real_ marks the variable as adaptive (C++ side uses NaN as sentinel).
  # A finite proposal_scale would instead pin the scale (no adaptation).
  expect_equal(param$defs$proposal_scale, NA_real_)
  expect_equal(param$defs$prior, "uniform")
  expect_equal(param$defs$min, "0|N")
  expect_equal(param$defs$max, "5|N")
})

test_that("Constant parameter definitions are correct", {
  param <- HousepiParameter$new()
  aux_add(param, name = "dummy") |>
    expect_error("Incomplete")
  aux_add(param, name = "dummy", prior = "uniform",
          min = 0, max = 5, const = TRUE) |>
    expect_error("Missing value")

  # load_list and single add should be coherent
  helper_test_aux_load(param, dummy = 3)
  tmp <- aux_add(HousepiParameter$new(), name = "dummy", value = 3)
  expect_equal(param$defs, tmp$defs)

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|N")
  expect_equal(param$defs$const, TRUE)
  expect_equal(param$defs$proposal, NA_character_)
  expect_equal(param$defs$proposal_scale, NA_real_)

  aux_add(param, name = "dummy", value = 3, prior = "uniform",
          min = 0, max = 5, const = TRUE) |>
    expect_warning("Overwriting")
  expect_equal(param$defs, tmp$defs)
})

test_that("Non default arguments are loaded correctly", {
  # simple add
  param <- HousepiParameter$new()
  aux_add(param, name = "dummy", prior = "gamma",
          shape = 2, scale = 1)
  expect_equal(param$defs$prior, "gamma|shape:scale|{2|N}:{1|N}")

  # load list
  expect_no_error(helper_test_aux_load(param,
    normal = list(prior = "gamma", shape = 2, scale = 1,
                  proposal = "normal"),
    lognormal = list(prior = "gamma", shape = 2, scale = 1,
                     proposal = "lognormal")
  ))
})

test_that("Symbolic definitions are equivalent", {
  param <- HousepiParameter$new()
  val <- 1
  helper_test_aux_load(
    param,
    var1 = list(prior = "normal", mean = indiv_shape + 1),
    var2 = list(prior = "normal", mean = indiv_shape + !!val),
    var3 = list(prior = "normal", mean = ~ indiv_shape + 1),
    var4 = list(prior = "normal", mean = ~ indiv_shape + !!val)
  )
  aux <- param$defs
  aux$name <- NULL

  # check that all rows are identical
  expect_equal(nrow(unique(aux)), 1)
})

test_that("Distribution weight labels are not swallowed by formal matching", {
  # list/bulk path: the reported reproducer, where single-letter weight 'c'
  # used to be partial-matched to the 'const' formal of add_quo
  param <- HousepiParameter$new() |>
    helper_test_aux_load(
      big_boss = list(prior = "categorical", a = 0.2, b = 0.1, c = 0.7)
    )
  expect_equal(param$defs$name, "big_boss")
  expect_equal(param$defs$prior,
               "categorical|a:b:c|{0.2|N}:{0.1|N}:{0.7|N}")

  # direct single add_quo form with c/d/p labels colliding with
  # const/discr/prior formals
  param2 <- HousepiParameter$new()
  param2$add_quo(name = "ocelot", prior = "categorical",
                 a = 0.2, b = 0.1, c = 0.7, d = 0.4, p = 0.6)
  expect_equal(param2$defs$name, "ocelot")
  expect_equal(param2$defs$prior,
               "categorical|a:b:c:d:p|{0.2|N}:{0.1|N}:{0.7|N}:{0.4|N}:{0.6|N}")
})
