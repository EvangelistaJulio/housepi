library(testthat)

aux_add <- function(self, name = NA_character_, value = NULL,
                    status = NA_integer_, const = NA, discr = FALSE,
                    func = NULL, func_args = NULL,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_) {
  self$add_quo(
    name = name, value = rlang::enquo(value),
    status = status, const = const, discr = discr,
    func = rlang::enquo(func),
    func_args = rlang::enquos(func_args),
    prior = prior, rlang::enquos(...),
    min = rlang::enquo(min), max = rlang::enquo(max),
    proposal = proposal
  )
  invisible(self)
}

test_that("Augmented definition error are recognized", {
  HousepiAugmented$new() |>
    aux_add("infec@") |>
    expect_error("invalid characters")
  HousepiAugmented$new() |>
    aux_add("infec") |>
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
  HousepiAugmented$new() |>
    aux_add("test", proposal = "prior", value = NA) |>
    expect_error("prior.*test")
  HousepiAugmented$new() |>
    aux_add("test", value = NA) |>
    expect_error("Missing value")

  augm <- HousepiAugmented$new() |>
    aux_add("test", value = 3, const = TRUE, proposal = "prior")
  expect_true(is.na(augm$defs$proposal))
})

test_that("Augmented variable definitions are correct", {
  # basic constant definition
  augm <- HousepiAugmented$new() |>
    aux_add(name = "dummy", value = 3)
  expect_equal(augm$defs$name, "dummy")
  # status is stored raw (character) to carry labels; NA when generic
  expect_equal(augm$defs$status, NA_character_)
  expect_equal(augm$defs$const, TRUE)
  expect_equal(augm$defs$discr, FALSE)
  expect_equal(augm$defs$value, "3|N")
  expect_equal(augm$defs$proposal, NA_character_)

  # variable definition
  augm <- HousepiAugmented$new() |>
    aux_add(
      name = "dummy", value = expo@D, prior = "uniform",
      min = expo_min@A, max = expo_max@P
    )
  expect_equal(augm$defs$name, "dummy")
  expect_equal(augm$defs$status, NA_character_)
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
  augm <- HousepiAugmented$new() |>
    aux_add(name = "dummy", prior = "norm")

  expect_equal(augm$defs$const, FALSE)
  expect_equal(augm$defs$proposal, "rw")
})

test_that("Augmented independence proposal folds its arguments", {
  augm <- HousepiAugmented$new() |>
    aux_add(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_mean = 0, proposal_sd = 1
    )
  expect_equal(augm$defs$proposal, "normal|mean:sd|{0|N}:{1|N}")

  # random-walk scale is not accepted by an independence sampler
  augm <- HousepiAugmented$new() |>
    aux_add(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_scale = 0.5
    ) |>
    expect_error("random walk")
})

test_that("Independence proposal reject if unsupported prior domain", {
  hhmodel <- data.frame(house_id = "philanthropy", status = 1, x = 3) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_augmented("z",
      prior = "uniform", min = theta@P, max = 5,
      proposal = "gamma", proposal_shape = 2
    ) |>
    evaluate_log_lik() |>
    expect_error("cannot cover")
})

test_that("Independence proposal inherits the prior location shift", {
  # a proposal from another family must sample in the prior's absolute frame:
  # with loc = min@D the support is [min, min + 1], so every draw lands there
  hhmodel <- data.frame(
    id = c("tommy", "gina"), house_id = "house",
    status = 1, min = c(4, 8)
  ) |>
    new_housepi_model("test") |>
    add_augmented("test_date",
      prior = "uniform", loc = min@D,
      proposal = "gamma", proposal_shape = 2
    )
  expect_warning(
    results <- run_inference(hhmodel,
      steps = 20, adapt = 0, chains = 1,
      threads = 1, progress = FALSE,
      id_draws = c("tommy", "gina")
    ),
    regexp = "below 400"
  )
  draws <- as.data.frame(get_augmented_draws(results))
  expect_true(all(draws$tommy_test_date >= 4 & draws$tommy_test_date <= 5))
  expect_true(all(draws$gina_test_date >= 8 & draws$gina_test_date <= 9))
})

test_that("A moving prior parameter keeps augmented draws inside their live support", {
  # regression: a dynamic prior parameter (loc/min) shifts the support without changing
  # the consumer's value, so the acceptance ratio must feel that shift (guarded); else the
  # window drifts off the value and the sampler emits zero-prior draws (~15% pre-fix).
  df <- data.frame(id = seq(2), house_id = "shadow_moses", status = 1)

  # augmented -> augmented: test_date ~ Uniform(loc = otacon@A), support [otacon, otacon + 1],
  # so every draw must satisfy 0 <= test_date - otacon <= 1
  augm <- df |>
    new_housepi_model("test") |>
    add_augmented(
      "otacon"    = list(prior = "uniform"),
      "test_date" = list(prior = "uniform", loc = otacon@A)
    ) |>
    run_inference(
      steps = 1000, adapt = 200, thin = 1, chains = 1,
      threads = 1, progress = FALSE, id_draws = "1"
    ) |>
    suppressMessages()
  d <- as.data.frame(get_augmented_draws(augm))
  gap <- d[["1_test_date"]] - d[["1_otacon"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # parameter -> augmented: z ~ Uniform(min = meryl@P, max = 5), support [meryl, 5],
  # so every draw must satisfy z >= meryl
  par <- df |>
    new_housepi_model("test") |>
    add_parameter("meryl", prior = "uniform") |>
    add_augmented("z", prior = "uniform", min = meryl@P, max = 5) |>
    run_inference(
      steps = 1000, adapt = 200, thin = 1, chains = 1,
      threads = 1, progress = FALSE, id_draws = "1"
    ) |>
    suppressMessages()
  z <- as.data.frame(get_augmented_draws(par))[["1_z"]]
  meryl <- as.data.frame(get_parameter_draws(par))[["meryl"]]
  expect_true(all(z >= meryl))
})

test_that("Augmented function definition is correct", {
  augm <- HousepiAugmented$new()
  aux_add(augm) |>
    expect_error("non-empty")

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
    C = c(
      prior = "uniform",
      min = !!external_value, max = log(infec@A + !!external_value)
    ),
    D = list(
      prior = "uniform",
      min = !!external_value, max = log(infec@A + !!external_value)
    )
  )
  expect_equal(subset(augm$defs, name == "A")$value, str_value)
  expect_equal(
    subset(augm$defs, name == "B")$func,
    sprintf("log|x|{sum|x:x2|{infec|A}:{%s}}", str_value)
  )
  subset(augm$defs, name == "C")
})

###############################
## AUGMENTED : ACTIVE FIELDS ##
###############################

test_that("Active fields are read-only", {
  HousepiAugmented$new() |>
    helper_test_read_only(
      "family", "augm_required", "data_required",
      "group_rates_required", "param_required",
      "observation_names", "house_observ_names",
      "status_names", "status_dictionary"
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
    some_func = list(
      func = "gamma",
      func_args = c(shape = ~shape, scale = ~scale)
    )
  )
  helper_test_aux_load(augm, external, append = FALSE) |>
    expect_error("rlang injection")

  helper_test_aux_load(augm, !!!external)
  expect_equal(augm$defs$func, "gamma|shape:scale|{shape|A}:{scale|A}")
})

test_that("Observation definitions are encoded", {
  augm <- HousepiAugmented$new() |>
    helper_test_aux_load(
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
  expect_identical(augm$observation_names, "onset_date")
})

test_that("Observation node rejects conflicting options", {
  augm <- HousepiAugmented$new()
  # prior and obs are mutually exclusive
  helper_test_aux_load(augm,
    x = list(prior = "normal", obs = "normal"),
    .kind = "obs"
  ) |>
    expect_error("only one between")
  # observation nodes are data-bound, never sampled
  helper_test_aux_load(augm,
    x = list(obs = "normal", proposal = "prior"),
    .kind = "obs"
  ) |>
    expect_warning("[Ii]gnored.*proposal")
})

test_that("Deterministic observation node is encoded as a function obs", {
  augm <- HousepiAugmented$new() |>
    helper_test_aux_load(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # deterministic observation (no 'obs')
    helper_test_aux_load(
      onset_date = floor(infec_time + incub_period),
      .kind = "obs",
      append = TRUE
    )
  d <- augm$defs[augm$defs[["name"]] == "onset_date", ]
  # tagged as observation, carrying a deterministic function and no distribution
  expect_identical(d[["kind"]], "obs")
  expect_true(is.na(d[["prior"]]))
  expect_true(is.na(d[["value"]]))
  expect_match(d[["func"]], "floor")
  # still reported as a generated data column
  expect_identical(augm$observation_names, "onset_date")
})

test_that("Observation and augmented loads each own only their definitions", {
  augm <- HousepiAugmented$new() |>
    helper_test_aux_load(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # non-append observation definition keeps augmented definitions
    helper_test_aux_load(
      onset_date = list(obs = "normal", mean = infec_time@A, sd = 1),
      .kind = "obs", append = FALSE
    )
  expect_setequal(
    augm$defs[["name"]],
    c("incub_period", "infec_time", "onset_date")
  )

  # symmetrically, augmented definitions keep observations
  augm |>
    helper_test_aux_load(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      append = FALSE
    )
  expect_true("onset_date" %in% augm$defs[["name"]])
  expect_identical(augm$observation_names, "onset_date")
})

test_that("Household observations do not interfer with augmented defs", {
  augm <- HousepiAugmented$new() |>
    helper_test_aux_load(
      incub_period = list(prior = "gamma", shape = 2, scale = 1)
    ) |>
    helper_test_aux_load(
      onset_date = list(obs = "normal", mean = incub_period@A, sd = 1),
      .kind = "obs"
    )

  # non-append household observation keeps augmented definitions
  augm |>
    helper_test_aux_load(
      first_test = min(test_date@D),
      .kind = "house", append = FALSE
    )
  augm$defs[["name"]] |>
    expect_setequal(c("incub_period", "onset_date", "first_test"))

  augm$defs[["kind"]][augm$defs[["name"]] == "first_test"] |>
    expect_identical("house")
  augm$defs[["func"]][augm$defs[["name"]] == "first_test"] |>
    expect_identical("min|x|{test_date|D}")

  # non-append load drops only the house rows
  augm |>
    helper_test_aux_load(
      last_test = max(test_date@D),
      .kind = "house", append = FALSE
    )
  expect_false("first_test" %in% augm$defs[["name"]])
  augm$defs[["name"]] |>
    expect_setequal(c("incub_period", "onset_date", "last_test"))
})

test_that("Observation-aware initialization", {
  window <- 0.1
  # augmented z ~ U(0, 1)
  # observation x ~ U(z - w, z + w)
  # if x = 0.5, valid for z in [0.4, 0.6]
  # without observation-aware restart, it would fail with 80% proba per-chain
  # probably of failing on K chains: 1 - 0.2^K
  data.frame(house_id = "shadow_moses", status = 1, x = 0.5) |>
    new_housepi_model("test") |>
    add_augmented("z", prior = "uniform") |>
    add_observation("x",
      obs = "uniform", loc = z,
      min = -!!window, max = !!window
    ) |>
    run_inference(steps = 10L, adapt = 0L, chains = 40L, progress = FALSE) |>
    suppressMessages() |>
    suppressWarnings() |>
    expect_no_error()

  # datum x = 5 is unreachable
  # UPDATE: it's the same as above: new-dataset decorator?
  data.frame(house_id = "outer_heaven", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_augmented("z", prior = "uniform") |>
    add_observation("x",
      obs = "uniform", loc = z,
      min = -!!window, max = !!window
    ) |>
    run_inference() |>
    suppressMessages() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_observation("x", obs = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    suppressMessages() |>
    expect_error("free parameter dependencies")

  # generic error
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", value = 0.3) |>
    add_augmented("z", prior = "uniform") |>
    add_observation("x", obs = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    expect_error("outside of respective domain")
})

test_that("Compressed households forced unique by augmented sampling size the cache", {
  # identical households compress to one eq-class, augmented sampling forces them
  # unique: the per-eq-class cache must grow to match (else out-of-bounds segfault)
  data.frame(status = rep(1, 5)) |>
    new_housepi_model("test") |>
    add_augmented("z", prior = "exponential", rate = 0.5) |>
    run_inference(steps = 10L, adapt = 0L, chains = 2L, progress = FALSE) |>
    suppressMessages() |>
    suppressWarnings() |>
    expect_no_error()
})

test_that("Reading an observation node's column: priors and loops are rejected", {
  # a deterministic readout of an observation is a legitimate downstream value;
  # only a prior reading an observation, or a dependency loop, is illegal.

  # allowed: a deterministic augmented that just reads the observed value
  data.frame(house_id = "shadow_moses", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", obs = "uniform", min = 0, max = 1) |>
    add_augmented("ocelot", value = onset@D) |>
    expect_no_error()

  # rule 1: a prior must not read an observation (illegal in both modes; caught
  # at eager build, so the construction itself errors)
  data.frame(house_id = "outer_heaven", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", obs = "uniform", min = 0, max = 1) |>
    add_augmented("latent", prior = "normal", mean = onset@D) |>
    expect_error("prior.*reads.*observation")

  # rule 2: a deterministic read that closes a loop (the observation's support
  # depends on the reader) is circular
  data.frame(house_id = "metal_gear", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_augmented("ocelot", value = onset@D) |>
    add_observation("onset", obs = "uniform", min = ocelot@A, max = 10) |>
    expect_error("[Cc]yclic")
})

test_that("A simulation override (value + prior) rejects forward dependencies", {
  # a value paired with a prior pins the variable in simulation only, when it is first
  # drawn -- so it may read data or parameters, never an observation generated later
  # nor another latent (which may not be drawn yet).

  # reads an observation
  data.frame(house_id = "shadow_moses", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", obs = "uniform", min = 0, max = 1) |>
    add_augmented("liquid", prior = "normal", mean = 0, sd = 1,
                  proposal = "rw", value = onset@D) |>
    expect_error("[Ss]imulation override.*observation")

  # reads another augmented (no observation present: guard still fires)
  data.frame(house_id = "outer_heaven", status = 1) |>
    new_housepi_model("test") |>
    add_augmented(
      liquid = list(prior = "normal", mean = 0, sd = 1, proposal = "rw",
                    value = solid@A),
      solid = list(prior = "gamma", shape = 2, scale = 1, proposal = "rw")
    ) |>
    expect_error("[Ss]imulation override.*augmented")
})

test_that("Group/contact rate keyed on an observation is rejected in simulation only", {
  # the rate feeds transmission, which the observation is generated from: fine in
  # inference (observed data), circular when generating during simulation.
  hhmodel <- data.frame(house_id = c(1, 1), status = c(1, 0), severity = c(2, 3)) |>
    new_housepi_model("final_size") |>
    add_group_rate("susceptibility", by = "severity") |>
    add_observation("severity", obs = "uniform", min = 0, max = 5) |>
    add_parameter_template()

  hhmodel |>
    run_simulation(nrep = 1L) |>
    suppressMessages() |>
    expect_error("group rate.*observation")
})

######################################
## AUGMENTED : END-TO-END MACHINERY ##
######################################

test_that("Augmented aliasing works as expected", {
  # check that model arguments are correctly integrated
  simulated <- data.frame(id = seq(1, 1000), house_id = 1, status = 1) |>
    new_housepi_model("test") |>
    add_augmented(
      hour = list(prior = "uniform", min = day, max = day + 1),
      day = list(prior = "gamma", shape = 2, scale = 1)
    ) |>
    run_simulation()
  expect_true(all(
    simulated$day <= simulated$hour & simulated$hour <= simulated$day + 1
  ))
})

test_that("A value paired with a prior overrides the draw in simulation only", {
  # simulation pins the sampled variable to its data-bound value (the wide prior is a
  # nuisance law that would scatter it); inference is unaffected -- the value is ignored
  # there, so the draws span the prior instead of collapsing onto the value.
  hhmodel <- data.frame(id = seq_len(400), house_id = 1, status = 1, anchor = 7) |>
    new_housepi_model("test") |>
    add_augmented(
      snake = list(prior = "uniform", min = 0, max = 100,
                   proposal = "rw", value = anchor@D)
    )

  # simulation: every drawn value collapses onto the override
  expect_true(all(run_simulation(hhmodel)$snake == 7))

  # inference: the value does not pin the sampler (draws span the prior, not the override)
  draws <- suppressWarnings( # short smoke run trips the convergence-size diagnostic
    hhmodel |>
      run_inference(steps = 40, chains = 1, id_draws = "1", progress = FALSE) |>
      get_augmented_draws()
  )
  expect_gt(length(unique(draws[["1_snake"]])), 1)
  expect_false(all(draws[["1_snake"]] == 7))
})

test_that("Load augmented list works fine for external lists", {
  hhmodel <- data.frame(
    id = "id1", house_id = "h1",
    status = 1, date = as.Date("2021-01-01")
  ) |>
    new_housepi_model("test")

  # an external list of definitions can be imported
  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  tmp <- add_augmented(hhmodel, import = external_list)$augmented

  hhmodel <- add_augmented(hhmodel, z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})

test_that("add_augmented_template loads", {
  sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_augmented_template() |>
    get_augmented() |>
    nrow() |>
    expect_gt(0L)
})

test_that("Augmented traces export only relevant variables", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(alpha = 1, beta = list(prior = "uniform")) |>
    add_augmented(
      # constants
      blocked = 42,
      # private variables
      .hidden = list(prior = "uniform"),
      .secret = list(prior = "normal"),
      # public non constant function
      layer = .secret - .hidden,
      # equal to parameters
      alpha = alpha@P,
      alpha_beta = alpha@P + beta@P,
      # equal to other augmented objects
      dummy = blocked,
      alias = layer,
      # implicit function
      impl = list(func = "gamma", func_args = c(shape = 3)),
      # constant function
      invar = list(func = "gamma", func_args = c(x = dummy, shape = 3)),
      # function
      expl = list(func = "gamma", func_args = c(x = alias, shape = 3))
    )
  infec_df <- subset(hhmodel$data_wide, status > 0)
  result <- suppressWarnings(run_inference(hhmodel,
    steps = 200, burn = 0, thin = 1, adapt = 0, chains = 1,
    id_draws = infec_df$id
  ))
  expect_equal(
    posterior::variables(result$augmented_draws),
    paste(infec_df$id, c("layer", "alias", "alpha_beta", "expl"), sep = "_")
  )
})

test_that("Augmented degenerate range warns and collapses", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter("d_dog", prior = "uniform") |>
    add_augmented("d_horse", prior = "uniform", min = 5, max = 5)
  expect_warning(
    tryCatch(
      run_inference(hhmodel, steps = 10, chains = 1, adapt = 0, thin = 1),
      error = function(e) NULL
    ),
    regexp = "[Dd]egenerate range|[Cc]ollapsed"
  )
})

test_that("Nested definitions works as expected", {
  hhmodel <- new_housepi_model(model = "test") |>
    # define Box-Muller transformation
    add_augmented(
      norm = list(prior = "normal"),
      u1 = list(prior = "uniform"),
      u2 = list(prior = "uniform"),
      z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    )
  expect_snapshot(print(hhmodel))

  results <- run_inference(hhmodel,
    steps = 1e5, chains = 1, adapt = 0,
    thin = 1, progress = 0
  )
  # no parameters nor priors, only augmented
  expect_equal(posterior::nvariables(results$parameter_draws), 0L)
  expect_null(results$.__enclos_env__$private$param_prior_draws_)
  expect_equal(
    posterior::nvariables(results$augmented_draws),
    dim(hhmodel$augmented)[1]
  )

  # test S3 version of summarise_draws
  out_vars <- posterior::summarise_draws(results)[["variable"]]
  expect_true(all(stringr::str_starts(out_vars, "tommy")))

  # empirical tests for goodness of Gaussian distribution
  sample_norm <- results$augmented_draws$tommy_norm
  sample_boxmul <- results$augmented_draws$tommy_z0

  expect_equal(mean(sample_norm), 0, tolerance = 0.05)
  expect_equal(sd(sample_norm), 1, tolerance = 0.05)
  probs <- c(0.025, 0.5, 0.975)
  q <- setNames(stats::quantile(sample_norm, probs), NULL)
  expect_equal(q, qnorm(probs), tolerance = 0.05)

  expect_equal(mean(sample_boxmul), 0, tolerance = 0.05)
  expect_equal(sd(sample_boxmul), 1, tolerance = 0.05)
  q <- setNames(stats::quantile(sample_boxmul, probs), NULL)
  expect_equal(q, qnorm(probs), tolerance = 0.05)
})

test_that("Nested augmented definitions raise errors", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("test")

  hhmodel |>
    add_augmented(
      .layer = list(prior = "uniform", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      .dummy = list(prior = "uniform"),
      z = .layer
    ) |>
    evaluate_log_lik() |>
    expect_no_error()

  # simple nest
  hhmodel |>
    add_augmented(.dummy = z, z = .dummy) |>
    evaluate_log_lik() |>
    expect_error("cyclic definition")

  # triangular nest
  hhmodel |>
    add_augmented(z = .layer, .layer = .dummy, .dummy = z) |>
    evaluate_log_lik() |>
    expect_error("cyclic definition")

  # first generation of definitions
  hhmodel |>
    add_augmented(
      .layer = list(prior = "normal", status = 1),
      .layer = list(value = .dummy, status = 2),
      #> NO .layer provided for status 3 and 4
      .dummy = list(prior = "lognormal", status = 2),
      z = .layer
    ) |>
    evaluate_log_lik() |>
    expect_error("Could not find a definition.*status")

  # second generation of definitions
  hhmodel |>
    add_augmented(
      .layer = list(prior = "normal", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      #> NO .dummy provided for status 3, 4
      .dummy = list(prior = "lognormal", status = c(2, 3)),
      z = .layer,
      append = FALSE
    ) |>
    evaluate_log_lik() |>
    expect_error("Could not find a definition.*status")
})

test_that("Missing augmented function arguments raise errors", {
  hhmodel <- data.frame(house_id = seq(3), id = seq(3), status = 1) |>
    new_housepi_model("test")

  # non default distribution arguments must be supplied
  hhmodel |>
    add_augmented(name = "z", func = "gamma") |>
    expect_error("Missing non default")
})

test_that("Status-specific augmented variables draw from per-status ranges", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("test") |>
    add_augmented(
      z = list(status = 1, prior = "uniform", min = 0, max = 1),
      z = list(status = 2, prior = "uniform", min = 1, max = 2),
      z = list(status = 3, prior = "uniform", min = 2, max = 3),
      z = list(status = 4, prior = "uniform", min = 3, max = 4)
    )

  # simulation routes draws into per-status supports
  withr::with_seed(1, sim <- run_simulation(hhmodel, nrep = 50))
  sim <- do.call(rbind, sim)
  for (s in seq(4)) {
    vals <- sim$z[sim$status == s]
    expect_true(all(vals >= s - 1 & vals <= s))
  }

  # short-run diagnostic fires (steps = 5 < 400); adapt budget becomes burn-in
  expect_warning(
    results <- hhmodel |>
      run_inference(
        steps = 5, thin = 1, adapt = 2,
        chain = 1, threads = 1,
        id_draws = hhmodel$data_wide$id
      ),
    regexp = "below 400"
  )
  draws <- as.data.frame(get_augmented_draws(results))
  for (s in seq(4)) {
    col <- draws[[as.character(s)]]
    expect_true(all(col >= s - 1 & col <= s))
  }
})

test_that("export_augmented_data works after evaluate_log_lik", {
  hhmodel <- data.frame(house_id = 1, id = c("a", "b"), status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_augmented(name = "z", prior = "uniform")
  evaluate_log_lik(hhmodel)
  result <- hhmodel$export_augmented_data()
  expect_true(is.data.frame(result))
  expect_true("z" %in% colnames(result))
})

test_that("Incubation augmentation recast runs on the test model", {
  # a minimal analogue of the data-augmentation vignette: the incubation gamma
  # enters as an augmented prior, the reported onset day as an observation node.
  set.seed(1)
  n <- 30
  exposure_max <- sample(0:3, n, replace = TRUE)
  infec <- runif(n, 0, exposure_max + 0.99)
  onset <- infec + rgamma(n, shape = 2, scale = 1.1)
  df <- data.frame(
    status = 1, exposure_max = exposure_max,
    onset_min = floor(onset), onset_max = floor(onset)
  )
  df$onset_mid <- (df$onset_min + df$onset_max) / 2 + 0.5

  results <- expect_no_error(
    df |>
      new_housepi_model("test") |>
      add_parameter(
        mean = list(prior = "uniform", min = 0.1, max = 5, proposal = "rw"),
        var = list(prior = "uniform", min = 0.1, max = 5, proposal = "rw")
      ) |>
      add_augmented(
        .infec_time = list(
          prior = "uniform", max = exposure_max@D + 0.99,
          proposal = "rw"
        ),
        incub_period = list(
          prior = "gamma", mean = mean@P, var = var@P,
          min = 0.5, max = 14, proposal = "rw"
        ),
        .onset_time = .infec_time + incub_period
      ) |>
      add_observation("onset_mid",
        obs = "normal",
        mean = .onset_time, sd = 0.3
      ) |>
      run_inference(
        steps = 20, adapt = 0, chains = 1, thin = 1,
        progress = FALSE
      ) |>
      suppressMessages() |>
      suppressWarnings()
  )

  # the gamma incubation prior contributes to the (finite) augmented log-prior
  lprior <- get_complete_draws(results)$lprior
  expect_true(all(is.finite(lprior)))
})
