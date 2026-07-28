library(testthat)

# model carrying no parameter of its own.
mk_param <- local({
  list(model = new_housepi_model(model = "test"))
})

##############################
## PARAMETER : DEFINITIONS  ##
##############################

test_that("Default parameter definitions are correct", {
  mk_param$model |>
    add_parameter() |>
    expect_error("[Nn]o.*definitions.*provided")

  mk_param$model |>
    add_parameter(name = "dummy", value = 3, const = FALSE) |>
    expect_error("[Nn]o .prior. provided for parameter")

  # default sample framework
  param <- mk_param$model |>
    add_parameter(
      dummy = list(value = 3, prior = "uniform", min = 0, max = 5)
    ) |>
    get_parameters()

  expect_equal(param$name, "dummy")
  expect_equal(param$value, 3)
  expect_equal(param$const, FALSE)
  expect_equal(param$discr, FALSE)
  expect_equal(param$proposal, "rw")
  expect_equal(param$prior, "uniform")
  expect_equal(param$min, 0)
  expect_equal(param$max, 5)

  # bulk and single add should be coherent
  tmp <- mk_param$model |>
    add_parameter(
      name = "dummy", value = 3, prior = "uniform",
      min = 0, max = 5
    ) |>
    get_parameters()
  expect_equal(param, tmp)

  # proposal
  param <- mk_param$model |>
    add_parameter(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "rw_log"
    ) |>
    get_parameters()

  expect_equal(param$name, "dummy")
  expect_equal(param$value, NA_real_)
  expect_equal(param$const, FALSE)
  expect_equal(param$discr, FALSE)
  # a bare proposal (no folded scale) stays adaptive
  expect_equal(param$proposal, "rw_log")
  expect_equal(param$prior, "uniform")
  expect_equal(param$min, 0)
  expect_equal(param$max, 5)

  # a user-defined proposal_scale is folded into the proposal call as an
  # argument of the DSL, pinning the scale (no adaptation)
  param <- mk_param$model |>
    add_parameter(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "rw", proposal_scale = 0.5
    ) |>
    get_parameters()
  expect_equal(param$proposal, "rw")
  expect_equal(param$proposal_args[[1]], c(scale = 0.5))
})

test_that("Proposal families route correctly", {
  # a distribution name is an independence sampler; its arguments fold in
  param <- mk_param$model |>
    add_parameter(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_mean = 0, proposal_sd = 1
    ) |>
    get_parameters()
  expect_equal(param$proposal, "normal")
  expect_equal(param$proposal_args[[1]], c(mean = 0, sd = 1))

  # sampling from the variable's own distribution uses the "forward" keyword
  param <- mk_param$model |>
    add_parameter(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "forward"
    ) |>
    get_parameters()
  expect_equal(param$proposal, "forward")

  # proposal_scale only tunes a random walk
  mk_param$model |>
    add_parameter(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_scale = 0.5
    ) |>
    expect_error("random.?walk")

  # distribution arguments are not random-walk knobs
  mk_param$model |>
    add_parameter(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal = "rw", proposal_sd = 1
    ) |>
    expect_error("[Rr]andom.?walk proposal .*accepts only")

  # advanced options need an explicit proposal
  mk_param$model |>
    add_parameter(
      name = "otacon", prior = "uniform", min = 0, max = 5,
      proposal_mean = 0
    ) |>
    expect_error("[Aa]dvanced .proposal.* require an explicit")
})

test_that("Constant parameter definitions are correct", {
  mk_param$model |>
    add_parameter(name = "dummy") |>
    expect_error("[Ii]ncomplete")
  mk_param$model |>
    add_parameter(
      name = "dummy", prior = "uniform",
      min = 0, max = 5, const = TRUE
    ) |>
    expect_error("[Mm]issing value")

  # bulk and single add should be coherent
  hhmodel <- add_parameter(mk_param$model, dummy = 3)
  tmp <- add_parameter(mk_param$model, name = "dummy", value = 3)
  expect_equal(get_parameters(hhmodel), get_parameters(tmp))

  param <- get_parameters(hhmodel)
  expect_equal(param$name, "dummy")
  expect_equal(param$value, 3)
  expect_equal(param$const, TRUE)
  expect_equal(param$proposal, NA_real_)

  expect_warning(
    overwritten <- add_parameter(
      hhmodel, name = "dummy", value = 3, prior = "uniform",
      min = 0, max = 5, const = TRUE
    ),
    "[Oo]verwriting"
  )
  expect_equal(get_parameters(overwritten), get_parameters(tmp))
})

test_that("Non default arguments are loaded correctly", {
  # simple add
  param <- mk_param$model |>
    add_parameter(
      name = "dummy", prior = "gamma",
      shape = 2, scale = 1
    ) |>
    get_parameters()
  expect_equal(param$prior, "gamma")
  expect_equal(param$args[[1]], c(shape = 2, scale = 1))
})

test_that("Variadic labels are not swallowed by formal matching", {
  param <- mk_param$model |>
    add_parameter(
      vulcan = list(prior = "categorical", a = 0.2, b = 0.1, c = 0.7),
      raven = list(prior = "categorical", d = 0.4, p = 0.6)
    ) |>
    get_parameters()
  expect_equal(param$name, c("vulcan", "raven"))
  expect_equal(param$prior, c("categorical", "categorical"))
  expect_equal(param$args[[1]], c(a = 0.2, b = 0.1, c = 0.7))
  expect_equal(param$args[[2]], c(d = 0.4, p = 0.6))
})

############################
## PARAMETER : END-TO-END ##
############################

test_that("Unsupported parameter definitions throw errors", {
  hhmodel <- new_housepi_model(model = "test")

  # eager rebuild: errors surface at definition, not at run_inference
  hhmodel |>
    add_parameter("a", prior = "uniform", min = min_val@D) |>
    expect_error("[Nn]ot supported")
  # hint is thrown
  hhmodel |>
    add_parameter("b", prior = "uniform", min = min_val) |>
    expect_error("[Rr]lang injection")
})

test_that("Degenerate range raises an error for parameters", {
  new_housepi_model(model = "test") |>
    add_parameter("venom", prior = "uniform", min = 3, max = 3) |>
    expect_error("[Dd]egenerate range")
})

test_that("Invalid proposal name is rejected", {
  new_housepi_model(model = "test") |>
    add_parameter("snake", prior = "uniform", proposal = "venom") |>
    run_inference(steps = 1, chains = 1, progress = FALSE) |>
    expect_error("[Uu]nknown proposal.*venom")
})

test_that("Lognormal proposal non-positive domain", {
  # prior with unbounded lower bound
  new_housepi_model(model = "test") |>
    add_parameter("huey", prior = "normal", proposal = "rw_log") |>
    expect_error("non.?positive.*huey")

  # prior with negative lower bound
  new_housepi_model(model = "test") |>
    add_parameter(
      name = "miller", prior = "uniform",
      min = -1, proposal = "rw_log"
    ) |>
    expect_error("non.?positive")

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
        proposal = "forward", min = !!val$min, max = !!val$max
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

  df <- as.data.frame(posterior::as_draws_df(res$parameter_draws))
  # every draw of beta equals hhsize_ref^(-gamma) exactly
  expect_equal(df$beta, 4^(-df$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hhmodel <- new_housepi_model(model = "test")

  # parameters do not accept individual-specifc definitons
  hhmodel |>
    add_augmented("cigar", model = "normal") |>
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

test_that("Parameter template import loads the shipped definitions", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template")

  defs <- get_parameters(hhmodel)
  expect_setequal(defs$name, c("alpha", "beta", "gamma", "hhsize_ref"))

  # a model shipping no template definition file says so
  flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_augmented(import = "template") |>
    expect_error("file.*not found")

  # the test model carries no template at all: warn and add nothing
  hhmodel <- new_housepi_model(model = "test")
  expect_warning(
    result <- add_parameter(hhmodel, import = "template"),
    regexp = "test model"
  )
  expect_identical(nrow(get_parameters(result)), 0L)
})
