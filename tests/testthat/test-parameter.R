library(testthat)

####  Helpers  ####

# model carrying no parameter of its own
mk_param <- local({
  list(model = new_housepi_model(model = "test"))
})

# transmission model reading its parameters, bare and with its template
mk_final_size <- local({
  model <- flu_1918_frederick |> new_housepi_model("final-size")
  list(model = model, template = model |> add_parameter(import = "template"))
})

# case-ascertained model with every shipped template
mk_case_ascertained <- local({
  model <- sarscov2_2021_icofs |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_infectivity_profile(import = "template") |>
    add_latent(import = "template")
  list(model = model)
})

####  Definitions  ####

test_that("Default parameter definitions are correct", {
  mk_param$model |>
    add_parameter() |>
    expect_error("[Mm]issing.*empty.*definitions")

  # default sample framework
  param <- mk_param$model |>
    add_parameter("rho", prior = "uniform", max = 5) |>
    get_parameters()

  param$name |> expect_identical("rho")
  param$const |> expect_all_na()
  param$discr |> expect_identical(FALSE)
  param$proposal |> expect_all_na()
  param$prior |> expect_identical("uniform")
  param$min |> expect_equal(0)
  param$max |> expect_equal(5)

  # bulk and single modes definitions should be identical
  mk_param$model |>
    add_parameter(rho = list(prior = "uniform", max = 5)) |>
    get_parameters() |>
    expect_identical(param)

  # non-default prior arguments are folded into the prior call
  param <- mk_param$model |>
    add_parameter("rho", prior = "gamma", shape = 2, scale = 1) |>
    get_parameters()
  param$prior |> expect_identical("gamma")
  param$args[[1]] |> expect_equal(c(shape = 2, scale = 1))

  # a user-defined proposal_scale is folded into the proposal call as an
  # argument of the DSL, pinning the scale (no adaptation)
  param <- mk_param$model |>
    add_parameter("rho", prior = "uniform", max = 5,
                  proposal = "rw_log", proposal_scale = 0.5) |>
    get_parameters()
  param$proposal |> expect_identical("rw_log")
  param$proposal_args[[1]] |> expect_equal(c(scale = 0.5))
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
  by_const <- mk_param$model |> add_parameter(rho = list(const = 21))
  get_parameters(by_func) |> expect_identical(get_parameters(by_const))

  # a deterministic definition accepts no proposal
  mk_param$model |>
    add_parameter(rho = list(func = "log", x = 2, proposal = "rw")) |>
    expect_error("[Dd]eterministic parameter.*rho.*cannot be sampled")

  # the spare arguments reach the function, whichever spelling declares it
  hh_model <- mk_param$model |> add_parameter("rho_2", const = 2)
  by_args <- hh_model |> add_parameter("rho", func = "log", x = rho_2)
  by_expr <- hh_model |> add_parameter(rho = log(rho_2))
  get_parameters(by_args) |> expect_identical(get_parameters(by_expr))
  get_parameters(by_args)$func_args[[1]] |> expect_identical(c(x = "rho_2|P"))
})

test_that("A parameter takes a single defining argument", {
  mk_param$model |>
    add_parameter("rho", func = "log", x = 2, prior = "uniform") |>
    expect_error("[Oo]nly one of .func. and .prior.")
  mk_param$model |>
    add_parameter("rho", const = 2, func = "log", x = 1) |>
    expect_error("[Oo]nly one of .const. and .func.")
})

test_that("Spare arguments reach the function they belong to", {
  # neither a prior nor a func takes them
  mk_param$model |>
    add_parameter("rho", const = 3, mean = 1, sd = 2) |>
    expect_error("2 spare arguments for parameter .rho.: .mean. and .sd.")
  mk_param$model |>
    add_parameter("rho", const = 3, mean = 1) |>
    expect_error("[Pp]rovide the .prior. or the .func. they belong to")
  # a func folded into a constant takes none either
  mk_param$model |>
    add_parameter("rho", func = log(3), mean = 1) |>
    expect_error("1 spare argument for parameter .rho.: .mean.")

  # a parameter function is checked against the registry at declaration
  mk_param$model |>
    add_parameter("rho", func = "log", base = 2) |>
    expect_error("Unused argument .base.")
  # a distribution is not a parameter function
  mk_param$model |>
    add_parameter("rho", func = "gamma", shape = 2, scale = 1) |>
    expect_error("Unknown function .gamma.")

  # an argument given twice
  mk_param$model |>
    add_parameter(rho = list(prior = "normal", mean = 0, mean = 5)) |>
    expect_error("[Rr]epeated argument .mean. for parameter .rho.")
})

test_that("Proposal families route correctly", {
  # a distribution name is an independence sampler; its arguments fold in
  param <- mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "normal", proposal_mean = 0, proposal_sd = 1) |>
    get_parameters()
  param$proposal |> expect_identical("normal")
  param$proposal_args[[1]] |> expect_equal(c(mean = 0, sd = 1))

  # sampling from the variable's own distribution uses the "forward" keyword
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 5,
                  proposal = "forward") |>
    get_parameters() |>
    getElement("proposal") |>
    expect_identical("forward")

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
    expect_error(".const. of parameter .rho. must be the value")

  # bulk and single add should be coherent
  hh_model <- mk_param$model |> add_parameter(rho = 3)
  hh_model_2 <- mk_param$model |> add_parameter("rho", const = 3)
  get_parameters(hh_model) |> expect_identical(get_parameters(hh_model_2))

  param <- get_parameters(hh_model)
  param$name |> expect_identical("rho")
  param$const |> expect_equal(3)
  param$proposal |> expect_equal(NA_real_)

  # a constant takes no prior, proposal or bounds
  mk_param$model |>
    add_parameter("rho", const = 3, prior = "uniform") |>
    expect_error("[Oo]nly one of .const. and .prior.")
  mk_param$model |>
    add_parameter("rho", const = 3, proposal = "rw") |>
    expect_error("[Cc]onstant parameter .rho. cannot be sampled")
  mk_param$model |>
    add_parameter("rho", const = 3, min = 0, max = 5) |>
    expect_error("[Cc]onstant parameter .rho. has no distribution to bound")
  # a numeric func alike
  mk_param$model |>
    add_parameter("rho", func = 3, proposal = "rw") |>
    expect_error("[Cc]onstant parameter .rho. cannot be sampled")
  mk_param$model |>
    add_parameter("rho", func = 3, min = 0) |>
    expect_error("[Cc]onstant parameter .rho. has no distribution to bound")
  mk_param$model |>
    add_parameter("rho", func = 3, prior = "uniform") |>
    expect_error("[Oo]nly one of .func. and .prior.")
})

test_that("Parameter options typed as values refuse other types", {
  mk_param$model |>
    add_parameter("rho", prior = sum) |>
    expect_error("[Aa]rgument .prior. must be character, numeric, or logical")
  mk_param$model |>
    add_parameter("rho", prior = "uniform", discr = NULL) |>
    expect_error("[Aa]rgument .discr. must be character, numeric, or logical")
  mk_param$model |>
    add_parameter(rho = list(prior = sum)) |>
    expect_error("[Aa]rgument .prior. must be character, numeric, or logical")
})

test_that("A prior argument outside its domain is rejected at definition", {
  mk_param$model |>
    add_parameter(rho = list(prior = "normal", sd = -1)) |>
    expect_error("sd.*normal.*must be greater than 0")
  mk_param$model |>
    add_parameter(rho = list(prior = "gamma", shape = -2)) |>
    expect_error("shape.*gamma.*must be greater than 0")
})

####  End-to-end  ####

test_that("Unsupported parameter definitions throw errors", {
  # eager rebuild: errors surface at definition, not at run_inference
  data.frame(id = 1, hh_id = 1, status = 1, datum = 1) |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform", min = datum) |>
    expect_error("[Nn]ot supported")
  # a bare name nothing declares carries no source yet: it waits for the wall
  hh_model <- add_parameter(mk_param$model, "rho", prior = "uniform",
                            min = datum)
  get_parameters(hh_model)$min |> expect_identical("datum|?")
  hh_model |>
    evaluate_log_lik() |>
    expect_error("[Cc]ould not resolve the source of .datum.")

  # a nested function missing a non-default argument
  mk_param$model |>
    add_parameter(rho = list(prior = "uniform"), rho_2 = exp(pow(rho))) |>
    evaluate_log_lik() |>
    expect_error("[Nn]on default argument .x2.*cannot be missing")
  mk_param$model |>
    add_parameter(rho = list(prior = "uniform"),
                  rho_2 = exp(pow(x2 = rho))) |>
    evaluate_log_lik() |>
    expect_error("[Aa]rgument .x. of function .pow. cannot be missing")

  # strings outside the naming convention
  mk_param$model |>
    add_parameter("rho", prior = "normal|x") |>
    expect_error("[Pp]rior .normal.x. contained invalid")
  mk_param$model |>
    add_parameter("rho", prior = "uniform", proposal = "rw|x") |>
    expect_error("[Pp]roposal .rw.x. contained invalid")

  # a constant argument outside the function domain
  mk_param$model |>
    add_parameter(rho = list(func = "log", x = -1)) |>
    evaluate_log_lik() |>
    expect_error("[Vv]alue -1 is outside of domain")
})

test_that("A degenerate domain raises an error for parameters", {
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = 3, max = 3) |>
    expect_error("[Dd]egenerate domain.*.min. \\(3\\) equals .max. \\(3\\)")
})

test_that("Invalid proposal name is rejected", {
  mk_param$model |>
    add_parameter("rho", prior = "uniform", proposal = "unknown") |>
    run_inference(steps = 1L, chains = 1L) |>
    expect_error("[Uu]nknown proposal .unknown. for parameter.*add_parameter")
  # a random-walk name is suggested alongside the distributions, the context
  # closes the first line
  mk_param$model |>
    add_parameter("rho", prior = "uniform", proposal = "rww") |>
    expect_error(paste0("[Uu]nknown proposal .rww. for parameter .rho.\\.\\n",
                        ".*[Dd]id you mean .rw."))
})

test_that("Lognormal proposal rejects unsupported domains", {
  # prior with unbounded lower bound
  mk_param$model |>
    add_parameter("rho", prior = "normal", proposal = "rw_log") |>
    expect_error("rw_log.*non-negative minimum.*rho")

  # prior with negative lower bound
  mk_param$model |>
    add_parameter("rho", prior = "uniform", min = -1, proposal = "rw_log") |>
    expect_error("rw_log.*non-negative minimum")

  # discrete lognormal proposal is forbidden
  mk_param$model |>
    add_parameter("rho", prior = "gamma", shape = 2, scale = 1,
                  discr = TRUE, proposal = "rw_log") |>
    expect_error("rw_log.*discrete prior.*rho")

  # location parameter
  draws <- mk_param$model |>
    add_parameter("rho", prior = "gamma", shape = 2,
                  max = 5, loc = -1, proposal = "rw_log") |>
    run_inference(steps = 20L, chains = 1L, adapt = 20L, thin = 1L) |>
    suppressWarnings() |>
    get_parameter_draws()
  (draws[["rho"]] > -1) |> expect_all_true()
  (draws[["rho"]] < 4) |> expect_all_true()
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
                   proposal = "rw", min = !!val$min, max = !!val$max),
      # prior quartiles coincide: the walks start from a unit scale
      rho_4 = list(prior = "poisson", lambda = 0.1, proposal = "rw_unif",
                   min = !!val$min, max = !!val$max),
      rho_5 = list(prior = "poisson", lambda = 0.1, proposal = "rw",
                   min = !!val$min, max = !!val$max)
    ) |>
    run_inference(steps = 50L, thin = 1L, adapt = 50L, chains = 1L) |>
    suppressWarnings() |>
    get_parameter_draws() |>
    draws_to_frame()
  (draws == floor(draws)) |> as.vector() |> expect_all_true()
  (draws >= val$min) |> as.vector() |> expect_all_true()
  (draws <= val$max) |> as.vector() |> expect_all_true()

  # a lattice kernel tuned narrower than one integer value stays pinned to its
  # starting point, every move proposing the value it already holds
  length(unique(draws[, "rho_2"])) |> expect_gt(1L)
  length(unique(draws[, "rho_3"])) |> expect_gt(1L)
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
    (0.5 * sum(abs(empirical - exact))) |> expect_lt(0.06)
  }
})

test_that("Default random walks sample their target", {
  # no proposal declared: each parameter moves in the coordinate of its domain,
  # logit on (0, 1), log above a minimum, log below a maximum
  # no model reads them: the posterior is the prior, the draws match its mean
  # a lost Jacobian draws beta(1, 4), gamma of shape 1, the normal piled at 0
  draws <- mk_param$model |>
    add_parameter(rho = list(prior = "beta", shape1 = 2, shape2 = 5),
                  rho_2 = list(prior = "gamma", shape = 2, scale = 1),
                  rho_3 = list(prior = "normal", mean = 0, sd = 1, max = 0),
                  rho_4 = list(prior = "loguniform", min = 1e-6, max = 1)) |>
    run_inference(steps = 1e4L, thin = 1L, adapt = 500L, chains = 1L) |>
    suppressWarnings() |>
    get_parameter_draws()

  mean(draws[["rho"]]) |> expect_equal(2 / 7, tolerance = 0.06)
  mean(draws[["rho_2"]]) |> expect_equal(2, tolerance = 0.06)
  mean(draws[["rho_3"]]) |> expect_equal(-sqrt(2 / pi), tolerance = 0.06)
  # log of a loguniform on six decades is flat on (log 1e-6, 0)
  log_rho_4 <- log(draws[["rho_4"]])
  mean(log_rho_4) |> expect_equal(log(1e-6) / 2, tolerance = 0.06)
  quartiles <- quantile(log_rho_4, c(0.25, 0.5, 0.75), names = FALSE)
  (abs(quartiles - log(1e-6) * c(0.75, 0.5, 0.25)) < 0.7) |> expect_all_true()
})

test_that("A random walk adapts its scale unless the scale is a number", {
  # an NA scale is no scale: the walk starts from its prior and adapts
  hh_model <- mk_param$model |>
    add_parameter(rho = list(prior = "uniform", max = 5, proposal = "rw",
                             proposal_scale = NA),
                  rho_2 = list(prior = "uniform", max = 5, proposal = "rw",
                               proposal_scale = 0.5))
  hh_model |>
    run_inference(steps = 10L, thin = 1L, adapt = 200L, chains = 1L) |>
    suppressWarnings()

  moves <- inspect_state(hh_model)$.sampler$moves
  scale <- stats::setNames(moves$scale, moves$name)
  scale[["rho_2"]] |> expect_identical(0.5)
  # 1 is the placeholder a walk holds before its first reset
  (scale[["rho"]] == 1) |> expect_false()
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
  hh_model <- mk_final_size$model |>
    add_parameter(alpha = 1e-4, beta = hh_size_ref^-gamma,
                  gamma = list(prior = "uniform"), hh_size_ref = 4)
  get_parameters(hh_model)$func_args[[2]] |>
    expect_identical(c(x = "hh_size_ref|P", x2 = "-|x|[gamma|P]"))

  # both sampled and function-derived parameters appear in the trace
  draws <- hh_model |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings() |>
    get_parameter_draws() |>
    posterior::as_draws_df() |>
    as.data.frame()
  names(draws) |> expect_contains(c("gamma", "beta"))
  # every draw of beta equals hh_size_ref^(-gamma) exactly
  draws$beta |> expect_equal(4^(-draws$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hh_model <- add_latent(mk_param$model, "hidden", model = "normal")

  # parameters do not accept individual-specifc definitons
  hh_model |>
    add_parameter("rho", prior = "normal", mean = hidden + 1) |>
    expect_error("[Ii]ndividual.*specific")

  # a bare name declared as a latent resolves there, and is refused as such
  hh_model |>
    add_parameter(rho = hidden^-2) |>
    expect_error("rho.*cannot.*depend.*latent")

  # a data column is rejected as well
  mk_param$model |>
    add_dataset(dplyr::mutate(dummy_household, datum = 1)) |>
    add_parameter(rho = datum + 1) |>
    expect_error("rho.*cannot.*depend.*data")
})

test_that("Parameter template import loads the shipped definitions", {
  get_parameters(mk_final_size$template)$name |>
    expect_setequal(c("alpha", "beta", "gamma", "hh_size_ref"))

  # a model shipping no template definition file says so
  mk_final_size$model |>
    add_latent(import = "template") |>
    expect_error("file.*not found")

  # the test model carries no template at all: warn and add nothing
  expect_warning(
    hh_model <- add_parameter(mk_param$model, import = "template"),
    "test model"
  )
  get_parameters(hh_model) |> expect_nrow(0L)
})

test_that("A parameter meets the model definition requirements", {
  hh_model <- mk_case_ascertained$model

  # the household size reference is required strictly positive
  hh_model |>
    add_parameter(hh_size_ref = 0) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a discrete prior owns the zero bin
  hh_model |>
    add_parameter(hh_size_ref = list(prior = "poisson", lambda = 4)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("not strictly positive")

  # a non-integer reference is accepted, constant or sampled
  hh_model |>
    add_parameter(hh_size_ref = 3.3) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_all_finite()
  hh_model |>
    add_parameter(hh_size_ref = list(prior = "uniform", min = 2, max = 6)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_all_finite()
})

test_that("A transmission model moves gamma in the block of beta", {
  # block of each parameter after a short run, NA if it moves alone
  read_blocks <- function(hh_model) {
    hh_model |>
      run_inference(steps = 10L, adapt = 10L, chains = 1L, threads = 1L,
                    progress = FALSE, warn = FALSE)
    moves <- inspect_state(hh_model)$.sampler$moves
    stats::setNames(moves$block, moves$name)[is.na(moves$member)]
  }
  hh_model <- mk_final_size$template$copy()

  # beta and gamma sampled with the default proposal
  blocks <- read_blocks(hh_model)
  blocks[c("beta", "gamma")] |> expect_identical(c(beta = 1L, gamma = 1L))
  blocks[["alpha"]] |> expect_all_na()
  blocks <- read_blocks(mk_case_ascertained$model$copy())
  blocks[c("beta", "gamma")] |> expect_identical(c(beta = 1L, gamma = 1L))

  # a rate coupled to beta shares the block
  mask <- flu_2014_matsumoto$hh_id %in% unique(flu_2014_matsumoto$hh_id)[1:300]
  strata <- tibble::tribble(~from_member, ~to_member, ~param_name,
                            "father",     "mother",   ".reference",
                            NA,           NA,         "ctr_rest")
  blocks <- flu_2014_matsumoto[mask, ] |>
    new_housepi_model("final-size") |>
    add_contact_rate(strata = strata) |>
    add_parameter(import = "template") |>
    read_blocks()
  blocks[c("ctr_rest", "beta", "gamma")] |>
    expect_identical(c(ctr_rest = 1L, beta = 1L, gamma = 1L))

  # pinned beta: gamma moves alone
  hh_model |>
    add_parameter(beta = 0.1) |>
    suppressWarnings() |>
    read_blocks() |>
    expect_all_na()
})

test_that("Strict positivity propagates through a parameter function", {
  # the community rate is required strictly positive and carries no
  # discretization requirement, so it reads a function definition directly
  mk <- function(...) {
    mk_final_size$template |>
      add_parameter(rho = list(prior = "gamma", shape = 2, scale = 1), ...) |>
      suppressWarnings() |>
      evaluate_log_lik()
  }

  # positivity follows the operands, one node at a time
  mk(alpha = exp(-rho)) |> expect_all_finite()
  mk(alpha = sqrt(rho) / 10) |> expect_all_finite()
  mk(alpha = rho^2 * 1e-2) |> expect_all_finite()
  mk(rho_2 = list(prior = "gamma", shape = 2, scale = 1),
     alpha = min(rho, rho_2) * 1e-2) |>
    expect_all_finite()
  mk(rho_2 = list(prior = "normal"), alpha = max(rho, rho_2) * 1e-2) |>
    expect_all_finite()
  mk(rho_2 = list(prior = "gamma", shape = 2, scale = 1),
     alpha = (rho + rho_2) * 1e-2) |>
    expect_all_finite()
  # a bounded operand caps the unbounded side of the other
  mk(alpha = min(rho, 1e-2)) |> expect_all_finite()
  mk(alpha = min(1e-2, rho)) |> expect_all_finite()
  mk(rho_2 = list(prior = "normal"), alpha = max(rho_2, 1e-2) * rho) |>
    expect_all_finite()

  # a sign-indefinite node is not provably positive, whatever it evaluates to
  mk(alpha = log(rho)) |>
    expect_error("not strictly positive")
  mk(alpha = -rho) |>
    expect_error("not strictly positive")
  # min over a sign-indefinite operand loses the proof, max keeps it
  mk(rho_2 = list(prior = "normal"), alpha = min(rho, rho_2) * 1e-2) |>
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
    expect_all_finite()

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
                               loc = rho)) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("could not be interpreted as a numeric constant")
  # and the refusal quotes the whole expression it could not fold
  mk_final_size$template |>
    add_parameter(
      rho = list(prior = "uniform", min = 1e-3, max = 2e-3),
      alpha = list(prior = "uniform", min = -1e-4, max = 1e-4, loc = rho * 2)
    ) |>
    suppressWarnings() |>
    evaluate_log_lik() |>
    expect_error("[Ee]xpression .*rho.*could not be interpreted")
})

test_that("Every prior family and nested function samples across chains", {
  # one parameter per prior family, and functions nested as operands
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
      rho_9 = sin(rho), rho_10 = tan(rho), rho_11 = exp(rho),
      rho_12 = log(1 + exp(sub(sin(rho), tan(rho))))
    ) |>
    run_inference(steps = 5L, chains = 2L, adapt = 5L, thin = 1L) |>
    suppressWarnings() |>
    get_parameter_draws()

  posterior::nchains(draws) |> expect_identical(2L)
  posterior::variables(draws) |> expect_setequal(c("rho", paste0("rho_", 2:12)))
  # the functions track their operand within every chain
  draws <- as.data.frame(posterior::as_draws_df(draws))
  draws$rho_9 |> expect_equal(sin(draws$rho))
  draws$rho_11 |> expect_equal(exp(draws$rho))
  draws$rho_12 |> expect_equal(log(1 + exp(sin(draws$rho) - tan(draws$rho))))
})

test_that("An unpinned parameter is drawn from its prior in simulation", {
  # the proposal is an independence sampler parked in a corner of the prior's
  # domain: a generative draw reading it would pin every replicate there
  hh_model <- data.frame(id = seq_len(20), hh_id = 1, status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform", min = 0.1, max = 0.9,
                             proposal = "normal",
                             proposal_mean = 0.88, proposal_sd = 0.005)) |>
    add_latent(hidden = rho)

  # a deterministic latent reading the parameter exports its value
  drawn <- run_simulation(hh_model, nrep = 40) |>
    vapply(\(x) x$hidden[1], numeric(1))

  unique(drawn) |> expect_length(40L)
  (drawn >= 0.1 & drawn <= 0.9) |> expect_all_true()
  stats::sd(drawn) |> expect_gt(0.15) # uniform spread, not the proposal's 0.005
  mean(drawn) |> expect_lt(0.7)            # centred on the prior, not on 0.88
})

test_that("A parameter without a draws column is drawn from its prior", {
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform"),
                  rho_2 = list(prior = "uniform")) |>
    add_latent(hidden = rho, hidden_2 = rho_2)

  simulated <- run_simulation(hh_model, draws = data.frame(rho = rep(0.3, 5)))
  vapply(simulated, \(x) x$hidden[1], numeric(1)) |> expect_all_equal(0.3)
  drawn <- vapply(simulated, \(x) x$hidden_2[1], numeric(1))
  unique(drawn) |> expect_length(5L)
  # a single draw simulates a single frame
  simulated <- run_simulation(hh_model, draws = data.frame(rho = 0.3))
  simulated |> expect_s3_class("data.frame")
  simulated$hidden[1] |> expect_equal(0.3)
  simulated <- run_simulation(hh_model, draws = c(rho = 0.3))
  simulated |> expect_s3_class("data.frame")
  simulated$hidden[1] |> expect_equal(0.3)
  # posterior draws stored as a list simulate one frame per draw
  simulated <- run_simulation(
    hh_model, draws = posterior::as_draws_list(data.frame(rho = c(0.3, 0.4))),
    thin_ess = FALSE
  )
  vapply(simulated, \(x) x$hidden[1], numeric(1)) |> expect_equal(c(0.3, 0.4))
})

test_that("Per-replicate parameter values reach the latents reading them", {
  rho_table <- data.frame(rho = c(0, 10, 20))
  # a deterministic latent and a model argument, both reading rho
  # undeclared: the model is built on the table's first row
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = rho,
               hidden_2 = list(model = "normal", mean = rho, sd = 1e-3))
  # declared constant: the table overrides its value
  hh_model_2 <- hh_model |> add_parameter(rho = 0.3)

  for (simulated in list(run_simulation(hh_model, draws = rho_table),
                         run_simulation(hh_model_2, draws = rho_table))) {
    vapply(simulated, \(x) x$hidden[2], numeric(1)) |>
      expect_equal(rho_table$rho)
    vapply(simulated, \(x) x$hidden_2[2], numeric(1)) |>
      expect_equal(rho_table$rho, tolerance = 1e-2)
  }
})

test_that("A draws table is judged on every value it gives a parameter", {
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "exponential")) |>
    add_latent("hidden", model = "normal", sd = rho)

  # the prior keeps `sd` positive, the table does not
  hh_model |>
    run_simulation(draws = data.frame(rho = -1)) |>
    expect_error("`sd`.*must be greater than 0.*rho")
  # a table value overrides a constant
  hh_model |>
    add_parameter(rho = 1) |>
    suppressWarnings() |>
    run_simulation(draws = data.frame(rho = c(1, -1))) |>
    expect_error("`sd`.*must be greater than 0.*rho")
  # values outside the prior, inside what the latent accepts
  hh_model |>
    add_parameter(rho = list(prior = "uniform", min = 0.1, max = 0.2)) |>
    suppressWarnings() |>
    run_simulation(draws = data.frame(rho = c(0.5, 2))) |>
    vapply(\(x) x$hidden[1], numeric(1)) |>
    expect_all_finite()

  # model requests are judged on the table too
  mk_final_size$template |>
    run_simulation(draws = data.frame(hh_size_ref = c(4, 0))) |>
    expect_error("[Nn]ot strictly positive")
})

test_that("An unpinned parameter is drawn from its prior at evaluation", {
  hh_model <- mk_final_size$model |>
    add_parameter(beta = 0.1, gamma = 1, hh_size_ref = 4)
  log_lik <- hh_model |>
    add_parameter(alpha = list(prior = "uniform", min = 0.1, max = 0.9,
                               proposal = "normal",
                               proposal_mean = 0.88, proposal_sd = 0.005)) |>
    (\(x) replicate(20, evaluate_log_lik(x)))()

  # what the whole prior domain is worth, to scale the spread against
  pinned <- \(a) evaluate_log_lik(add_parameter(hh_model, alpha = !!a))
  span <- abs(pinned(0.9) - pinned(0.1))

  log_lik |> expect_all_finite()
  stats::sd(log_lik) |> expect_gt(0.05 * span)
})

test_that("A fixed seed reproduces parameter draws on a single thread", {
  # the normal prior and the default proposal both draw standard_normal:
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

  draws |> expect_identical(draws_at(1105L))
  # a seed that moves the draws, so the equality above is not vacuous
  identical(draws, draws_at(2705L)) |> expect_false()
})
