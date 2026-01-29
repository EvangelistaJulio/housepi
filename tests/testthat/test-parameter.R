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

aux_load <- function(self, ..., append = FALSE) {
  self$load_quosures(rlang::enquos(...), append = append)
  invisible(self)
}

test_that("Default parameter definitions are correct", {
  param <- housepi:::HousepiParameter$new()

  expect_error(
    aux_add(param, name = "dummy", value = 3, const = FALSE),
    regexp = "No prior"
  )

  # default sample framework
  aux_load(param, dummy = list(value = 3, prior = "uniform",
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
  aux_add(tmp, name = "dummy", value = 3, prior = "uniform",
          min = 0, max = 5)
  expect_equal(param$defs, tmp$defs)

  # proposal
  param <- housepi:::HousepiParameter$new()
  aux_add(param, name = "dummy", prior = "uniform", min = 0, max = 5,
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
  expect_error(aux_add(param, name = "dummy"), regexp = "Incomplete definition")
  expect_error(aux_add(param, name = "dummy", prior = "uniform",
                       min = 0, max = 5, const = TRUE))

  # load_list and single add should be coherent
  aux_load(param, dummy = 3)
  tmp <- aux_add(housepi:::HousepiParameter$new(), name = "dummy", value = 3)
  expect_equal(param$defs, tmp$defs)

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|N")
  expect_equal(param$defs$const, TRUE)
  expect_equal(param$defs$proposal, NA_character_)
  expect_equal(param$defs$proposal_scale, NA_real_)

  expect_warning(
    aux_add(param, name = "dummy", value = 3, prior = "uniform",
            min = 0, max = 5, const = TRUE),
    regexp = "Overwriting"
  )
  expect_equal(param$defs, tmp$defs)
})

test_that("Non default arguments are loaded correctly", {
  # simple add
  param <- housepi:::HousepiParameter$new()
  aux_add(param, name = "dummy", prior = "gamma",
          shape = 2, scale = 1)
  expect_equal(param$defs$prior, "gamma|shape:scale|{2|N}:{1|N}")

  # load list
  expect_no_error(aux_load(param,
    normal = list(prior = "gamma", shape = 2, scale = 1,
                  proposal = "normal"),
    lognormal = list(prior = "gamma", shape = 2, scale = 1,
                     proposal = "lognormal")
  ))
})

test_that("Symbolic definitions are equivalent", {
  param <- housepi:::HousepiParameter$new()
  val <- 1
  aux_load(param,
           var1 = list(prior = "normal", mean = indiv_shape + 1),
           var2 = list(prior = "normal", mean = indiv_shape + !!val),
           var3 = list(prior = "normal", mean = ~ indiv_shape + 1),
           var4 = list(prior = "normal", mean = ~ indiv_shape + !!val))
  aux <- param$defs
  aux$name <- NULL

  # check that all rows are identical
  expect_equal(nrow(unique(aux)), 1)
})
