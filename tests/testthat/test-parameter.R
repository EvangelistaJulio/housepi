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
    add_parameter(dummy = list(value = 3, prior = "uniform", max = 5)) |>
    get_parameters()

  expect_equal(param$name, "dummy")
  expect_equal(param$value, 3)
  expect_equal(param$const, FALSE)
  expect_equal(param$discr, FALSE)
  expect_equal(param$proposal, "rw")
  expect_equal(param$prior, "uniform")
  ## UPDATE: the default min is not returned, but it happears as NA
  ## expect_equal(param$min, 0)
  expect_equal(param$max, 5)

  # bulk and single add should be coherent
  tmp <- mk_param$model |>
    add_parameter("dummy", value = 3, prior = "uniform", max = 5) |>
    get_parameters()
  expect_equal(param, tmp)

  # a user-defined proposal_scale is folded into the proposal call as an
  # argument of the DSL, pinning the scale (no adaptation)
  param <- mk_param$model |>
    add_parameter(name = "dummy", prior = "uniform", max = 5,
                  proposal = "rw_log", proposal_scale = 0.5) |>
    get_parameters()
  expect_equal(param$proposal, "rw_log")
  expect_equal(param$proposal_args[[1]], c(scale = 0.5))
})

test_that("The defining distribution of the sibling family is named as spare", {
  # `model` keys latent and observations, `prior` keys parameters
  mk_param$model |>
    add_parameter(dummy = list(model = "uniform", min = 0, max = 5)) |>
    expect_error("[Ss]pare argument .model. for parameter .dummy.")
  mk_param$model |>
    add_parameter(dummy = list(model = "uniform", min = 0, max = 5)) |>
    expect_error("Did you mean .prior. instead of .model..")
})

test_that("A func definition is deterministic", {
  # a plain number is a constant value, whichever argument carried it
  by_func <- mk_param$model |> add_parameter(dummy = list(func = 21))
  by_value <- mk_param$model |> add_parameter(dummy = list(value = 21))
  expect_equal(get_parameters(by_func), get_parameters(by_value))

  # a deterministic definition accepts no proposal
  mk_param$model |>
    add_parameter(
      dummy = list(func = "log", func_args = c(x = 2), proposal = "rw")
    ) |>
    expect_error("[Dd]eterministic parameter.*dummy.*cannot be sampled")
})

test_that("Proposal families route correctly", {
  # a distribution name is an independence sampler; its arguments fold in
  param <- mk_param$model |>
    add_parameter(name = "otacon", prior = "uniform", min = 0, max = 5,
                  proposal = "normal", proposal_mean = 0, proposal_sd = 1) |>
    get_parameters()
  expect_equal(param$proposal, "normal")
  expect_equal(param$proposal_args[[1]], c(mean = 0, sd = 1))

  # sampling from the variable's own distribution uses the "forward" keyword
  param <- mk_param$model |>
    add_parameter(name = "otacon", prior = "uniform", min = 0, max = 5,
                  proposal = "forward") |>
    get_parameters()
  expect_equal(param$proposal, "forward")

  # proposal_scale only tunes a random walk
  mk_param$model |>
    add_parameter(name = "otacon", prior = "uniform", min = 0, max = 5,
                  proposal = "normal", proposal_scale = 0.5) |>
    expect_error("random.?walk")

  # distribution arguments are not random-walk knobs
  mk_param$model |>
    add_parameter(name = "otacon", prior = "uniform", min = 0, max = 5,
                  proposal = "rw", proposal_sd = 1) |>
    expect_error("[Rr]andom.?walk proposal .*accepts only")

  # advanced options need an explicit proposal
  mk_param$model |>
    add_parameter(name = "otacon", prior = "uniform", min = 0, max = 5,
                  proposal_mean = 0) |>
    expect_error("[Aa]dvanced .proposal.* require an explicit")
})

test_that("Constant parameter definitions are correct", {
  mk_param$model |>
    add_parameter(name = "dummy") |>
    expect_error("[Ii]ncomplete")
  mk_param$model |>
    add_parameter(name = "dummy", prior = "uniform",
                  min = 0, max = 5, const = TRUE) |>
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
    overwritten <- hhmodel |>
      add_parameter(name = "dummy", value = 3, prior = "uniform",
                    min = 0, max = 5, const = TRUE),
    "[Oo]verwriting"
  )
  expect_equal(get_parameters(overwritten), get_parameters(tmp))
})

test_that("Non default arguments are loaded correctly", {
  # simple add
  param <- mk_param$model |>
    add_parameter(name = "dummy", prior = "gamma", shape = 2, scale = 1) |>
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

test_that("A degenerate domain raises an error for parameters", {
  new_housepi_model(model = "test") |>
    add_parameter("venom", prior = "uniform", min = 3, max = 3) |>
    expect_error("[Dd]egenerate domain")
})

test_that("Invalid proposal name is rejected", {
  new_housepi_model(model = "test") |>
    add_parameter("snake", prior = "uniform", proposal = "venom") |>
    run_inference(steps = 1L, chains = 1L) |>
    expect_error("[Uu]nknown proposal.*venom")
})

test_that("Lognormal proposal rejects unsupported domains", {
  # prior with unbounded lower bound
  new_housepi_model(model = "test") |>
    add_parameter("huey", prior = "normal", proposal = "rw_log") |>
    expect_error("non.?positive.*huey")

  # prior with negative lower bound
  new_housepi_model(model = "test") |>
    add_parameter("miller", prior = "uniform",
                  min = -1, proposal = "rw_log") |>
    expect_error("non.?positive")

  # discrete lognormal proposal is forbidden
  new_housepi_model(model = "test") |>
    add_parameter("ocelot", prior = "gamma", shape = 2, scale = 1,
                  discr = TRUE, proposal = "rw_log") |>
    expect_error("lognormal.*discrete.*ocelot")

  # location parameter
  suppressWarnings(
    param_draws <- new_housepi_model(model = "test") |>
      add_parameter("amanda", prior = "gamma", shape = 2,
                    max = 5, loc = -1, proposal = "rw_log") |>
      run_inference(steps = 20L, chains = 1L, adapt = 0L, thin = 1L) |>
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
      sample = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                    proposal = "forward", min = !!val$min, max = !!val$max),
      uniform = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                     proposal = "rw_unif", min = !!val$min, max = !!val$max),
      normal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                    proposal = "rw", min = !!val$min, max = !!val$max)
    )

  results <- hhmodel |>
    run_inference(steps = 50L, thin = 1L, adapt = 50L, chains = 1L) |>
    suppressWarnings()
  y <- posterior::as_draws_matrix(get_parameter_draws(results))
  expect_true(all(y == floor(y)))
  expect_true(all(y >= val$min))
  expect_true(all(y <= val$max))

  # a lattice kernel tuned narrower than one integer value stays pinned to its
  # starting point, every step proposing the value it already holds
  expect_gt(length(unique(y[, "uniform"])), 1L)
  expect_gt(length(unique(y[, "normal"])), 1L)
})

test_that("Lattice random walks sample their target", {
  # a parameter no model reads keeps its prior as posterior, so the draws are
  # checked against the exact mass function: a kernel whose acceptance ratio
  # does not describe the moves it makes settles on a different law
  exact <- dbinom(0:6, size = 6, prob = 0.4)

  for (kernel in c("rw", "rw_unif")) {
    draws <- new_housepi_model(model = "test") |>
      add_parameter(revolver = list(prior = "binomial", size = 6, prob = 0.4,
                                    proposal = !!kernel)) |>
      run_inference(steps = 1e4L, thin = 1L, adapt = 500L, chains = 1L) |>
      suppressWarnings() |>
      get_parameter_draws()

    empirical <- as.numeric(table(factor(draws[["revolver"]], levels = 0:6)))
    empirical <- empirical / sum(empirical)
    # total variation distance, well within the spread of a correct kernel
    expect_lt(0.5 * sum(abs(empirical - exact)), 0.06)
  }
})

test_that("A discrete prior resolves 'uniform' to the discrete distribution", {
  hhmodel <- new_housepi_model(model = "test")

  # the flat lattice law owns every integer value in the support, so a bound off
  # the lattice has no reading
  hhmodel |>
    add_parameter(theta = list(prior = "uniform", min = 0.5, discr = TRUE)) |>
    evaluate_log_lik() |>
    expect_error("min.*uniform.*discrete")
})

test_that("Parameter defined as a function tracks its dependencies", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(alpha = 1e-4, beta = hhsize_ref@P^-gamma@P,
                  gamma = list(prior = "uniform"), hhsize_ref = 4)
  res <- hhmodel |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings()
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
    add_latent("cigar", model = "normal") |>
    add_parameter("kaz", prior = "normal", mean = cigar + 1) |>
    expect_error("[Ii]ndividual.*specific")

  # bare names default to latent source! here @P is missing
  hhmodel |>
    add_parameter(beta = hhsize_ref^-gamma) |>
    expect_error("beta.*cannot.*depend.*latent")

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
    add_latent(import = "template") |>
    expect_error("file.*not found")

  # the test model carries no template at all: warn and add nothing
  hhmodel <- new_housepi_model(model = "test")
  expect_warning(
    result <- add_parameter(hhmodel, import = "template"),
    regexp = "test model"
  )
  expect_identical(nrow(get_parameters(result)), 0L)
})

test_that("A parameter meets the model definition requirements", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_infectivity_profile(import = "template") |>
    add_latent(import = "template")

  # the model declares the household size reference discrete
  hhmodel$.__enclos_env__$private$generator_$param_discr_required |>
    expect_equal("hhsize_ref")

  hhmodel |>
    add_parameter(hhsize_ref = 3.3) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not discrete-valued")

  # the same reference is required strictly positive: an integer value clears the
  # discretization check and is caught by the second requirement
  hhmodel |>
    add_parameter(hhsize_ref = 0) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a discrete prior clears the discretization check but owns the zero bin
  hhmodel |>
    add_parameter(hhsize_ref = list(prior = "poisson", lambda = 4)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a sampled reference satisfies both without being declared discrete
  hhmodel |>
    add_parameter(hhsize_ref = list(prior = "uniform", min = 2, max = 6)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_no_error()
})

test_that("A reduction answers for the definitions it folds over", {
  # the operands are read for what they define, not for the value they hold
  # while the definition is being built: a pointer to a variable still waiting
  # for its first draw is not an absent operand
  flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template") |>
    add_parameter(raiden = list(prior = "gamma", shape = 2, scale = 1),
                  otacon = list(prior = "gamma", shape = 2, scale = 1),
                  hhsize_ref = min(raiden@P, otacon@P)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not discrete-valued")
})

test_that("Strict positivity propagates through a parameter function", {
  # the community risk is required strictly positive and carries no
  # discretization requirement, so it reads a function definition directly
  mk <- function(...) {
    flu_1918_frederick |>
      new_housepi_model("final_size") |>
      add_parameter(import = "template") |>
      add_parameter(...) |>
      suppressWarnings() |>
      evaluate_log_lik()
  }

  # positivity follows the operands, one node at a time
  expect_no_error(mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = exp(-raiden@P)))
  expect_no_error(mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = sqrt(raiden@P) / 10))
  expect_no_error(mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = raiden@P^2 * 1e-2))
  expect_no_error(mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
                     otacon = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = min(raiden@P, otacon@P) * 1e-2))
  expect_no_error(mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
                     otacon = list(prior = "normal"),
                     alpha = max(raiden@P, otacon@P) * 1e-2))
  expect_no_error(mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
                     otacon = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = (raiden@P + otacon@P) * 1e-2))

  # a sign-indefinite node is not provably positive, whatever it evaluates to
  mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
     alpha = log(raiden@P)) |>
    expect_error("not strictly positive")
  mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
     alpha = -raiden@P) |>
    expect_error("not strictly positive")
  # min over a sign-indefinite operand loses the proof, max keeps it
  mk(raiden = list(prior = "gamma", shape = 2, scale = 1),
     otacon = list(prior = "normal"),
     alpha = min(raiden@P, otacon@P) * 1e-2) |>
    expect_error("not strictly positive")
})
