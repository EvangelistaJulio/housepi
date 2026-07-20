library(testthat)

aux_add <- function(self, name = NA_character_, value = NULL,
                    const = NA, discr = FALSE,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_) {
  self$add_quo(
    name = name, value = rlang::enquo(value),
    const = const, discr = discr,
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal
  )
  invisible(self)
}

##############################
## PARAMETER : DEFINITIONS  ##
##############################

test_that("Default parameter definitions are correct", {
  param <- HousepiParameter$new()

  aux_add(param, name = "dummy", value = 3, const = FALSE) |>
    expect_error("prior")

  # active fields are read-only
  helper_test_read_only(param, "param_required")

  # default sample framework
  param |>
    helper_test_aux_load(
      dummy = list(value = 3, prior = "uniform", min = 0, max = 5)
    )

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|N")
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  expect_equal(param$defs$proposal, "rw")
  expect_equal(param$defs$prior, "uniform")
  expect_equal(param$defs$min, "0|N")
  expect_equal(param$defs$max, "5|N")

  # load_list and single add should be coherent
  tmp <- HousepiParameter$new() |>
    aux_add(
      name = "dummy", value = 3, prior = "uniform",
      min = 0, max = 5
    )
  expect_equal(param$defs, tmp$defs)

  # proposal
  param <- HousepiParameter$new() |>
    aux_add(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "rw"
    )

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, NA_real_)
  expect_equal(param$defs$const, FALSE)
  expect_equal(param$defs$discr, FALSE)
  # a bare proposal (no folded scale) stays adaptive
  expect_equal(param$defs$proposal, "rw")
  expect_equal(param$defs$prior, "uniform")
  expect_equal(param$defs$min, "0|N")
  expect_equal(param$defs$max, "5|N")

  # a user-defined proposal_scale is folded into the proposal call as an
  # argument of the DSL, pinning the scale (no adaptation)
  param <- HousepiParameter$new() |>
    aux_add(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "rw", proposal_scale = 0.5
    )
  expect_equal(param$defs$proposal, "rw|scale|{0.5|N}")
})

test_that("Proposal families route correctly", {
  # a distribution name is an independence sampler; its arguments fold in
  param <- HousepiParameter$new() |>
    aux_add(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_mean = 0, proposal_sd = 1
    )
  expect_equal(param$defs$proposal, "normal|mean:sd|{0|N}:{1|N}")

  # sampling from the prior is the default keyword
  param <- HousepiParameter$new() |>
    aux_add(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "prior"
    )
  expect_equal(param$defs$proposal, "prior")

  # proposal_scale only tunes a random walk
  param <- HousepiParameter$new() |>
    aux_add(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_scale = 0.5
    ) |>
    expect_error("random walk")

  # distribution arguments are not random-walk knobs
  param <- HousepiParameter$new() |>
    aux_add(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "rw", proposal_sd = 1
    ) |>
    expect_error("only")

  # advanced options need an explicit proposal
  param <- HousepiParameter$new() |>
    aux_add(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal_mean = 0
    ) |>
    expect_error("explicit")
})

test_that("Constant parameter definitions are correct", {
  param <- HousepiParameter$new()
  aux_add(param, name = "dummy") |>
    expect_error("Incomplete")
  param |>
    aux_add(
      name = "dummy", prior = "uniform",
      min = 0, max = 5, const = TRUE
    ) |>
    expect_error("Missing value")

  # load_list and single add should be coherent
  helper_test_aux_load(param, dummy = 3)
  tmp <- aux_add(HousepiParameter$new(), name = "dummy", value = 3)
  expect_equal(param$defs, tmp$defs)

  expect_equal(param$defs$name, "dummy")
  expect_equal(param$defs$value, "3|N")
  expect_equal(param$defs$const, TRUE)
  expect_equal(param$defs$proposal, NA_character_)

  param |>
    aux_add(
      name = "dummy", value = 3, prior = "uniform",
      min = 0, max = 5, const = TRUE
    ) |>
    expect_warning("Overwriting")
  expect_equal(param$defs, tmp$defs)
})

test_that("Non default arguments are loaded correctly", {
  # simple add
  param <- HousepiParameter$new() |>
    aux_add(
      name = "dummy", prior = "gamma",
      shape = 2, scale = 1
    )
  expect_equal(param$defs$prior, "gamma|shape:scale|{2|N}:{1|N}")
})

test_that("Variadic labels are not swallowed by formal matching", {
  # single definition
  param <- HousepiParameter$new()
  param$add_quo(
    name = "vulcan_raven", prior = "categorical",
    a = 0.2, b = 0.1, c = 0.7, d = 0.4, p = 0.6
  )
  expect_equal(param$defs$name, "vulcan_raven")
  expect_equal(
    param$defs$prior,
    "categorical|a:b:c:d:p|{0.2|N}:{0.1|N}:{0.7|N}:{0.4|N}:{0.6|N}"
  )

  # multiple definitions
  param <- HousepiParameter$new() |>
    helper_test_aux_load(
      vulcan = list(prior = "categorical", a = 0.2, b = 0.1, c = 0.7),
      raven = list(prior = "categorical", d = 0.4, p = 0.6)
    )
  expect_equal(param$defs$name, c("vulcan", "raven"))
  expect_equal(
    param$defs$prior,
    c(
      "categorical|a:b:c|{0.2|N}:{0.1|N}:{0.7|N}",
      "categorical|d:p|{0.4|N}:{0.6|N}"
    )
  )
})

############################
## PARAMETER : END-TO-END ##
############################

test_that("Unsupported parameter definitions throw errors", {
  hhmodel <- new_housepi_model(model = "test")

  # eager rebuild: errors surface at definition, not at run_inference
  hhmodel |>
    add_parameter("a", prior = "uniform", min = min_val@D) |>
    expect_error("not supported")
  # hint is thrown
  hhmodel |>
    add_parameter("b", prior = "uniform", min = min_val) |>
    expect_error("rlang injection")
})

test_that("Degenerate range raises an error for parameters", {
  new_housepi_model(model = "test") |>
    add_parameter("venom", prior = "uniform", min = 3, max = 3) |>
    expect_error("degenerate range")
})

test_that("Invalid proposal name is rejected", {
  new_housepi_model(model = "test") |>
    add_parameter("snake", prior = "uniform", proposal = "venom") |>
    run_inference(steps = 1, chains = 1, progress = FALSE) |>
    expect_error("Unknown proposal.*venom")
})

test_that("Lognormal proposal non-positive domain", {
  # prior with unbounded lower bound
  new_housepi_model(model = "test") |>
    add_parameter("huey", prior = "normal", proposal = "rw_log") |>
    expect_error("non-positive.*huey")

  # prior with negative lower bound
  new_housepi_model(model = "test") |>
    add_parameter(
      name = "miller", prior = "uniform",
      min = -1, proposal = "rw_log"
    ) |>
    expect_error("non-positive")

  # location parameter
  suppressWarnings(
    param_draws <- new_housepi_model(model = "test") |>
      add_parameter(
        name = "amanda", prior = "gamma", shape = 2,
        max = 5, loc = -1, proposal = "rw_log"
      ) |>
      run_inference(steps = 20, chains = 1, adapt = 0, thin = 1) |>
      get_parameter_draws()
  )
  expect_true(all(param_draws[["amanda"]] > -1))
  expect_true(all(param_draws[["amanda"]] < 4))
})

test_that("Rlang calls work as expected", {
  # rlang calls and injections are evaluated for parameters
  # here testing for equivalent symbolic definitions
  val <- 1
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(
      param_shape = 5,
      var1 = list(prior = "normal", mean = param_shape@p + 1),
      var2 = list(prior = "normal", mean = param_shape@p + !!val),
      var3 = list(prior = "normal", mean = ~ param_shape@P + 1),
      var4 = list(prior = "normal", mean = ~ param_shape@P + !!val)
    )
  aux <- hhmodel$parameters
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 2L)
})

test_that("Discretization is enabled as expected", {
  discr <- TRUE
  val <- list(min = 1, max = 5)
  hhmodel <- new_housepi_model(model = "test")
  hhmodel <- hhmodel |>
    add_parameter(
      sample = list(
        prior = "gamma", shape = 2, scale = 1, discr = TRUE,
        proposal = "prior", min = !!val$min, max = !!val$max
      ),
      uniform = list(
        prior = "gamma", shape = 2, scale = 1, discr = TRUE,
        proposal = "rw_unif", min = !!val$min, max = !!val$max
      ),
      normal = list(
        prior = "gamma", shape = 2, scale = 1, discr = TRUE,
        proposal = "rw", min = !!val$min, max = !!val$max
      ),
      lognormal = list(
        prior = "gamma", shape = 2, scale = 1, discr = TRUE,
        proposal = "rw_log", min = !!val$min, max = !!val$max
      )
    )

  results <- suppressWarnings(run_inference(hhmodel,
    steps = 50, burn = 0, thin = 1, adapt = 50, chains = 1
  ))
  y <- posterior::as_draws_matrix(get_parameter_draws(results))
  expect_true(all(y == floor(y)))
  expect_true(all(y >= val$min))
  expect_true(all(y <= val$max))
})

test_that("Parameter defined as a function tracks its dependencies", {
  set.seed(1)
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(
      alpha = 1e-4, beta = hhsize_ref@P^-gamma@P,
      gamma = list(prior = "uniform"), hhsize_ref = 4
    )
  res <- suppressWarnings(
    hhmodel |>
      run_inference(steps = 100, chains = 1, adapt = 0, thin = 1)
  )
  vars <- posterior::variables(res$parameter_draws)
  # both sampled and function-derived parameters appear in the trace
  expect_true(all(c("gamma", "beta") %in% vars))
  # the function-derived parameter has no prior contribution
  # UPDATE: this will not require private access in the future
  expect_false("lprior_beta" %in% posterior::variables(
    res$.__enclos_env__$private$param_prior_draws_
  ))

  df <- as.data.frame(posterior::as_draws_df(res$parameter_draws))
  # every draw of beta equals hhsize_ref^(-gamma) exactly
  expect_equal(df$beta, 4^(-df$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hhmodel <- new_housepi_model(model = "test")

  # parameters do not accept individual-specifc definitons
  hhmodel |>
    add_augmented("cigar", prior = "normal") |>
    add_parameter("kaz", prior = "normal", mean = cigar + 1) |>
    expect_error("[Ii]ndividual.*specific")

  # bare names default to augmented source! here @P is missing
  hhmodel |>
    add_parameter(beta = hhsize_ref^-gamma) |>
    expect_error("beta.*cannot.*depend.*augmented")

  # explicit reference to a data column is also rejected
  hhmodel |>
    add_data_wide(dplyr::mutate(dummy_household, hush_puppy = 1)) |> 
    add_parameter(beta = hush_puppy@D + 1) |>
    expect_error("beta.*cannot.*depend.*data")
})
