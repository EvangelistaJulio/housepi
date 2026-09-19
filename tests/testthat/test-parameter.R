library(testthat)

# model carrying no parameter of its own
mk_param <- local({
  list(model = new_housepi_model(model = "test"))
})

# transmission model reading its parameters, bare and with its template
mk_final_size <- local({
  model <- flu_1918_frederick |> new_housepi_model("final-size")
  list(model = model, template = model |> add_parameter(import = "template"))
})

##############################
## PARAMETER : DEFINITIONS  ##
##############################

test_that("Default parameter definitions are correct", {
  mk_param$model |>
    add_parameter() |>
    expect_error("[Mm]issing.*empty.*definitions")

  mk_param$model |>
    add_parameter("rho", value = 3, const = FALSE) |>
    expect_error("[Nn]o .prior. provided for parameter")

  # default sample framework
  param <- mk_param$model |>
    add_parameter("rho", value = 3, prior = "uniform", max = 5) |>
    get_parameters()

  expect_equal(param$name, "rho")
  expect_equal(param$value, 3)
  expect_equal(param$const, FALSE)
  expect_equal(param$discr, FALSE)
  expect_equal(param$proposal, "rw")
  expect_equal(param$prior, "uniform")
  ## UPDATE: the default min is not returned, but it happears as NA
  ## expect_equal(param$min, 0)
  expect_equal(param$max, 5)

  # bulk and single modes definitions should be identical
  mk_param$model |>
    add_parameter(rho = list(value = 3, prior = "uniform", max = 5)) |>
    get_parameters() |>
    expect_identical(param)

  # non-default prior arguments are folded into the prior call
  param <- mk_param$model |>
    add_parameter("rho", prior = "gamma", shape = 2, scale = 1) |>
    get_parameters()
  expect_equal(param$prior, "gamma")
  expect_equal(param$args[[1]], c(shape = 2, scale = 1))

  # a user-defined proposal_scale is folded into the proposal call as an
  # argument of the DSL, pinning the scale (no adaptation)
  param <- mk_param$model |>
    add_parameter("rho", prior = "uniform", max = 5,
                  proposal = "rw_log", proposal_scale = 0.5) |>
    get_parameters()
  expect_equal(param$proposal, "rw_log")
  expect_equal(param$proposal_args[[1]], c(scale = 0.5))
})

test_that("The defining distribution of the sibling family is named as spare", {
  # `model` keys latent and observations, `prior` keys parameters
  mk_param$model |>
    add_parameter(rho = list(model = "uniform", min = 0, max = 5)) |>
    expect_error("[Ss]pare argument .model. for parameter .rho.")
  mk_param$model |>
    add_parameter(rho = list(model = "uniform", min = 0, max = 5)) |>
    expect_error("[Dd]id you mean .prior. instead of .model.")
})

test_that("A func definition is deterministic", {
  # a plain number is a constant value, whichever argument carried it
  by_func <- mk_param$model |> add_parameter(rho = list(func = 21))
  by_value <- mk_param$model |> add_parameter(rho = list(value = 21))
  expect_equal(get_parameters(by_func), get_parameters(by_value))

  # a deterministic definition accepts no proposal
  mk_param$model |>
    add_parameter(rho = list(func = "log", x = 2, proposal = "rw")) |>
    expect_error("[Dd]eterministic parameter.*rho.*cannot be sampled")

  # the spare arguments reach the function, whichever spelling declares it
  hhmodel <- mk_param$model |> add_parameter("rho_2", value = 2)
  by_args <- hhmodel |> add_parameter("rho", func = "log", x = rho_2@P)
  by_expr <- hhmodel |> add_parameter(rho = log(rho_2@P))
  expect_equal(get_parameters(by_args), get_parameters(by_expr))
  expect_equal(get_parameters(by_args)$func_args[[2]], c(x = "rho_2|P"))
})

test_that("A parameter takes a single defining argument", {
  mk_param$model |>
    add_parameter("rho", func = "log", x = 2, prior = "uniform") |>
    expect_error("[Oo]nly one between .func. and .prior.")
  mk_param$model |>
    add_parameter("rho", value = 2, func = "log", x = 1) |>
    expect_error("[Oo]nly one between .value. and .func.")
})

test_that("Spare arguments reach the function they belong to", {
  # neither a prior nor a func takes them
  mk_param$model |>
    add_parameter("rho", value = 3, mean = 1, sd = 2) |>
    expect_error("2 spare arguments for parameter .rho.: .mean. and .sd.")
  mk_param$model |>
    add_parameter("rho", value = 3, mean = 1) |>
    expect_error("[Pp]rovide the .prior. or the .func. they belong to")
  # a func folded into a constant takes none either
  mk_param$model |>
    add_parameter("rho", func = log(3), mean = 1) |>
    expect_error("[Ee]xpression already carries its arguments")

  # a parameter function is checked against the registry at declaration
  mk_param$model |>
    add_parameter("rho", func = "log", base = 2) |>
    expect_error("Unused argument .base.")
  # a distribution is not a parameter function
  mk_param$model |>
    add_parameter("rho", func = "gamma", shape = 2, scale = 1) |>
    expect_error("Unknown function .gamma.")
})

test_that("Proposal families route correctly", {
  # a distribution name is an independence sampler; its arguments fold in
  param <- mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "normal", proposal_mean = 0, proposal_sd = 1) |>
    get_parameters()
  expect_equal(param$proposal, "normal")
  expect_equal(param$proposal_args[[1]], c(mean = 0, sd = 1))

  # sampling from the variable's own distribution uses the "forward" keyword
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "forward") |>
    get_parameters() |>
    getElement("proposal") |>
    expect_equal("forward")

  # proposal_scale only tunes a random walk
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "normal", proposal_scale = 0.5) |>
    expect_error("random.?walk")

  # distribution arguments are not random-walk knobs
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "rw", proposal_sd = 1) |>
    expect_error("[Rr]andom.?walk proposal .*accepts only")

  # advanced options need an explicit proposal
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal_mean = 0) |>
    expect_error("[Aa]dvanced .proposal.* require an explicit")

  # sampling from the variable's own distribution takes no tuning
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "forward", proposal_mean = 0) |>
    expect_error("[Ii]ndependence sampling from .forward. takes no")
})

test_that("Constant parameter definitions are correct", {
  mk_param$model |>
    add_parameter("rho") |>
    expect_error("[Ii]ncomplete")
  mk_param$model |>
    add_parameter("rho", prior = "uniform", const = TRUE) |>
    expect_error("[Mm]issing value")

  # bulk and single add should be coherent
  hhmodel <- mk_param$model |> add_parameter(rho = 3)
  hhmodel_2 <- mk_param$model |> add_parameter("rho", value = 3)
  expect_identical(get_parameters(hhmodel), get_parameters(hhmodel_2))

  param <- get_parameters(hhmodel)
  expect_equal(param$name, "rho")
  expect_equal(param$value, 3)
  expect_equal(param$const, TRUE)
  expect_equal(param$proposal, NA_real_)

  expect_warning(
    overwritten <- hhmodel |>
      add_parameter("rho", value = 3, prior = "uniform",
                    min = 0, max = 5, const = TRUE),
    "[Oo]verwriting"
  )
  expect_equal(get_parameters(overwritten), get_parameters(hhmodel_2))
})

test_that("Parameter options typed as values refuse other types", {
  mk_param$model |>
    add_parameter("rho", prior = sum) |>
    expect_error("[Aa]rgument .prior. must be character, numeric, or logical")
  mk_param$model |>
    add_parameter("rho", value = 1, const = NULL) |>
    expect_error("[Aa]rgument .const. must be character, numeric, or logical")
  mk_param$model |>
    add_parameter(rho = list(prior = sum)) |>
    expect_error("[Aa]rgument .prior. must be character, numeric, or logical")
})

############################
## PARAMETER : END-TO-END ##
############################

test_that("Unsupported parameter definitions throw errors", {
  # eager rebuild: errors surface at definition, not at run_inference
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = datum@D) |>
    expect_error("[Nn]ot supported")
  # a bare name nothing declares carries no source yet: it waits for the wall
  expect_no_error(
    hhmodel <- add_parameter(mk_param$model, "rho", prior = "uniform",
                             min = datum)
  )
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("[Cc]ould not resolve the source of .datum.")

  # a constant argument outside the function domain
  mk_param$model |>
    add_parameter(rho = list(func = "log", x = -1)) |>
    evaluate_log_lik() |>
    expect_error("[Vv]alue -1 for argument .x. is outside of domain")
})

test_that("A degenerate domain raises an error for parameters", {
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 3, max = 3) |>
    expect_error("[Dd]egenerate domain")
})

test_that("Invalid proposal name is rejected", {
  mk_param$model |>
    add_parameter("rho", prior = "uniform", proposal = "unknown") |>
    run_inference(steps = 1L, chains = 1L) |>
    expect_error("[Uu]nknown proposal.*unknown")
})

test_that("Lognormal proposal rejects unsupported domains", {
  # prior with unbounded lower bound
  mk_param$model |>
    add_parameter("rho", prior = "normal", proposal = "rw_log") |>
    expect_error("non.?positive.*rho")

  # prior with negative lower bound
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = -1, proposal = "rw_log") |>
    expect_error("non.?positive")

  # discrete lognormal proposal is forbidden
  mk_param$model |>
    add_parameter("rho", prior = "gamma", shape = 2, scale = 1,
                  discr = TRUE, proposal = "rw_log") |>
    expect_error("lognormal.*discrete.*rho")

  # location parameter
  draws <- mk_param$model |>
    add_parameter("rho", prior = "gamma", shape = 2,
                  max = 5, loc = -1, proposal = "rw_log") |>
    run_inference(steps = 20L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings() |>
    get_parameter_draws()
  expect_true(all(draws[["rho"]] > -1))
  expect_true(all(draws[["rho"]] < 4))
})

test_that("Discretization is enabled as expected", {
  val <- list(min = 1, max = 5)
  # one parameter per proposal: forward, lattice uniform, lattice normal walk
  draws <- mk_param$model |>
    add_parameter(
      rho = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                 proposal = "forward", min = !!val$min, max = !!val$max),
      rho_2 = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                   proposal = "rw_unif", min = !!val$min, max = !!val$max),
      rho_3 = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                   proposal = "rw", min = !!val$min, max = !!val$max)
    ) |>
    run_inference(steps = 50L, thin = 1L, adapt = 50L, chains = 1L) |>
    suppressWarnings() |>
    get_parameter_draws() |>
    posterior::as_draws_matrix()
  expect_true(all(draws == floor(draws)))
  expect_true(all(draws >= val$min))
  expect_true(all(draws <= val$max))

  # a lattice kernel tuned narrower than one integer value stays pinned to its
  # starting point, every step proposing the value it already holds
  expect_gt(length(unique(draws[, "rho_2"])), 1L)
  expect_gt(length(unique(draws[, "rho_3"])), 1L)
})

test_that("Lattice random walks sample their target", {
  # a parameter no model reads keeps its prior as posterior, so the draws are
  # checked against the exact mass function: a kernel whose acceptance ratio
  # does not describe the moves it makes settles on a different law
  exact <- dbinom(0:6, size = 6, prob = 0.4)

  for (kernel in c("rw", "rw_unif")) {
    draws <- mk_param$model |>
      add_parameter(rho = list(prior = "binomial", size = 6, prob = 0.4,
                               proposal = !!kernel)) |>
      run_inference(steps = 1e4L, thin = 1L, adapt = 500L, chains = 1L) |>
      suppressWarnings() |>
      get_parameter_draws()

    empirical <- as.numeric(table(factor(draws[["rho"]], levels = 0:6)))
    empirical <- empirical / sum(empirical)
    # total variation distance, well within the spread of a correct kernel
    expect_lt(0.5 * sum(abs(empirical - exact)), 0.06)
  }
})

test_that("A discrete prior resolves 'uniform' to the discrete distribution", {
  # the flat lattice law owns every integer value in the support, so a bound off
  # the lattice has no reading
  mk_param$model |>
    add_parameter(rho = list(prior = "uniform", min = 0.5, discr = TRUE)) |>
    evaluate_log_lik() |>
    expect_error("min.*uniform.*discrete")
})

test_that("Parameter defined as a function tracks its dependencies", {
  sigilled <- mk_final_size$model |>
    add_parameter(alpha = 1e-4, beta = hhsize_ref@P^-gamma@P,
                  gamma = list(prior = "uniform"), hhsize_ref = 4)
  # a bare reference encodes the same definition: no sigil is needed to reach
  # a parameter no other holder claims
  bare <- mk_final_size$model |>
    add_parameter(alpha = 1e-4, beta = hhsize_ref^-gamma,
                  gamma = list(prior = "uniform"), hhsize_ref = 4)
  expect_equal(get_parameters(bare), get_parameters(sigilled))
  expect_equal(get_parameters(bare)$func_args[[2]],
               c(x = "hhsize_ref|P", x2 = "minus|x|[gamma|P]"))
  # a sigil reads the same in lower case
  mk_final_size$model |>
    add_parameter(alpha = 1e-4, beta = hhsize_ref@p^-gamma@p,
                  gamma = list(prior = "uniform"), hhsize_ref = 4) |>
    get_parameters() |>
    expect_equal(get_parameters(sigilled))

  # both sampled and function-derived parameters appear in the trace
  draws <- sigilled |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings() |>
    getElement("parameter_draws") |>
    posterior::as_draws_df() |>
    as.data.frame()
  expect_true(all(c("gamma", "beta") %in% names(draws)))
  # every draw of beta equals hhsize_ref^(-gamma) exactly
  expect_equal(draws$beta, 4^(-draws$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hhmodel <- add_latent(mk_param$model, "hidden", model = "normal")

  # parameters do not accept individual-specifc definitons
  hhmodel |>
    add_parameter("rho", prior = "normal", mean = hidden + 1) |>
    expect_error("[Ii]ndividual.*specific")

  # a bare name declared as a latent resolves there, and is refused as such
  hhmodel |>
    add_parameter(rho = hidden^-2) |>
    expect_error("rho.*cannot.*depend.*latent")

  # explicit reference to a data column is also rejected
  mk_param$model |>
    add_data_wide(dplyr::mutate(dummy_household, datum = 1)) |>
    add_parameter(rho = datum@D + 1) |>
    expect_error("rho.*cannot.*depend.*data")
})

test_that("Parameter template import loads the shipped definitions", {
  get_parameters(mk_final_size$template)$name |>
    expect_setequal(c("alpha", "beta", "gamma", "hhsize_ref"))

  # a model shipping no template definition file says so
  mk_final_size$model |>
    add_latent(import = "template") |>
    expect_error("file.*not found")

  # the test model carries no template at all: warn and add nothing
  expect_warning(
    hhmodel <- add_parameter(mk_param$model, import = "template"),
    regexp = "test model"
  )
  expect_identical(nrow(get_parameters(hhmodel)), 0L)
})

test_that("A parameter meets the model definition requirements", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_infectivity_profile(import = "template") |>
    add_latent(import = "template")

  # the household size reference is required strictly positive
  hhmodel |>
    add_parameter(hhsize_ref = 0) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a discrete prior owns the zero bin
  hhmodel |>
    add_parameter(hhsize_ref = list(prior = "poisson", lambda = 4)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a non-integer reference is accepted, constant or sampled
  hhmodel |>
    add_parameter(hhsize_ref = 3.3) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_no_error()
  hhmodel |>
    add_parameter(hhsize_ref = list(prior = "uniform", min = 2, max = 6)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_no_error()
})

test_that("Strict positivity propagates through a parameter function", {
  # the community risk is required strictly positive and carries no
  # discretization requirement, so it reads a function definition directly
  mk <- function(...) {
    mk_final_size$template |>
      add_parameter(rho = list(prior = "gamma", shape = 2, scale = 1), ...) |>
      suppressWarnings() |>
      evaluate_log_lik()
  }

  # positivity follows the operands, one node at a time
  expect_no_error(mk(alpha = exp(-rho@P)))
  expect_no_error(mk(alpha = sqrt(rho@P) / 10))
  expect_no_error(mk(alpha = rho@P^2 * 1e-2))
  expect_no_error(mk(rho_2 = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = min(rho@P, rho_2@P) * 1e-2))
  expect_no_error(mk(rho_2 = list(prior = "normal"),
                     alpha = max(rho@P, rho_2@P) * 1e-2))
  expect_no_error(mk(rho_2 = list(prior = "gamma", shape = 2, scale = 1),
                     alpha = (rho@P + rho_2@P) * 1e-2))

  # a sign-indefinite node is not provably positive, whatever it evaluates to
  mk(alpha = log(rho@P)) |>
    expect_error("not strictly positive")
  mk(alpha = -rho@P) |>
    expect_error("not strictly positive")
  # min over a sign-indefinite operand loses the proof, max keeps it
  mk(rho_2 = list(prior = "normal"), alpha = min(rho@P, rho_2@P) * 1e-2) |>
    expect_error("not strictly positive")
})

test_that("Strict positivity propagates through a shifted prior", {
  # the support a prior places is loc + [min, max], and it is that composed
  # interval the positivity request reads

  # the shift carries an interval spanning both signs past zero
  mk_final_size$template |>
    add_parameter(alpha = list(prior = "uniform", min = -1e-4, max = 1e-4,
                               loc = 1e-3)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_no_error()

  # the same interval, unshifted, keeps no proof
  mk_final_size$template |>
    add_parameter(alpha = list(prior = "uniform", min = -1e-4, max = 1e-4)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a prior hyperparameter is a constant: it cannot read another parameter
  mk_final_size$template |>
    add_parameter(rho = list(prior = "uniform", min = 1e-3, max = 2e-3),
                  alpha = list(prior = "uniform", min = -1e-4, max = 1e-4,
                               loc = rho@P)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("could not be interpreted as a numeric constant")
  # and the refusal quotes the whole expression it could not fold
  mk_final_size$template |>
    add_parameter(rho = list(prior = "uniform", min = 1e-3, max = 2e-3),
                  alpha = list(prior = "uniform", min = -1e-4, max = 1e-4,
                               loc = rho@P * 2)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("[Ee]xpression .prod.*could not be interpreted")
})

test_that("Every prior and function is copied into each chain", {
  # chains sample copies of the definitions, one per prior family
  draws <- mk_param$model |>
    add_parameter(
      rho = list(prior = "beta", shape1 = 2, shape2 = 2),
      rho_2 = list(prior = "rayleigh", scale = 1),
      rho_3 = list(prior = "weibull", shape = 2, scale = 1),
      rho_4 = list(prior = "exponential", rate = 1,
                   proposal = "gamma", proposal_shape = 2),
      rho_5 = list(prior = "poisson", lambda = 3),
      rho_6 = list(prior = "binomial", size = 5, prob = 0.3),
      rho_7 = list(prior = "negative_binomial", size = 5, prob = 0.3),
      rho_8 = list(prior = "geometric", prob = 0.3),
      rho_9 = sin(rho@P), rho_10 = tan(rho@P), rho_11 = exp(rho@P)
    ) |>
    run_inference(steps = 5L, chains = 2L, adapt = 0L, thin = 1L) |>
    suppressWarnings() |>
    get_parameter_draws()

  expect_equal(posterior::nchains(draws), 2L)
  expect_setequal(posterior::variables(draws),
                  c("rho", paste0("rho_", 2:11)))
  # the functions track their operand within every chain
  draws <- as.data.frame(posterior::as_draws_df(draws))
  expect_equal(draws$rho_9, sin(draws$rho))
  expect_equal(draws$rho_11, exp(draws$rho))
})

test_that("An unpinned parameter is drawn from its prior in simulation", {
  # the proposal is an independence sampler parked in a corner of the prior's
  # domain: a generative draw reading it would pin every replicate there
  hhmodel <- data.frame(id = seq_len(20), hhid = 1, status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform", min = 0.1, max = 0.9,
                             proposal = "normal",
                             proposal_mean = 0.88, proposal_sd = 0.005)) |>
    add_latent(hidden = list(func = "rho@P"))

  # a deterministic latent reading the parameter exports the value drawn
  drawn <- run_simulation(hhmodel, nrep = 40) |>
    vapply(function(x) x$hidden[1], numeric(1))

  expect_true(all(drawn >= 0.1 & drawn <= 0.9))
  expect_gt(stats::sd(drawn), 0.15) # uniform spread, not the proposal's 0.005
  expect_lt(mean(drawn), 0.7)       # centred on the prior, not on 0.88
})

test_that("An unpinned parameter is drawn from its prior at evaluation", {
  hhmodel <- mk_final_size$model |>
    add_parameter(beta = 0.1, gamma = 1, hhsize_ref = 4)
  log_lik <- hhmodel |>
    add_parameter(alpha = list(prior = "uniform", min = 0.1, max = 0.9,
                               proposal = "normal",
                               proposal_mean = 0.88, proposal_sd = 0.005)) |>
    (\(x) replicate(20, evaluate_log_lik(x)))()

  # what the whole prior domain is worth, to scale the spread against
  pinned <- function(a) evaluate_log_lik(add_parameter(hhmodel, alpha = !!a))
  span <- abs(pinned(0.9) - pinned(0.1))

  expect_true(all(is.finite(log_lik)))
  expect_gt(stats::sd(log_lik), 0.05 * span)
})

test_that("A fixed seed reproduces parameter draws on a single thread", {
  # the normal prior and the default "rw" proposal both draw standard_normal:
  # sampler state surviving a re-seed desynchronises one run from the next
  draws_at <- function(seed) {
    set.seed(seed)
    mk_param$model |>
      add_parameter("rho", prior = "normal", mean = 2, sd = 1) |>
      run_inference(steps = 50L, chains = 1L, adapt = 0L, thin = 1L) |>
      suppressWarnings() |>
      get_parameter_draws()
  }
  draws <- draws_at(1105L)

  expect_identical(draws, draws_at(1105L))
  # a seed that moves the draws, so the equality above is not vacuous
  expect_false(identical(draws, draws_at(2705L)))
})
