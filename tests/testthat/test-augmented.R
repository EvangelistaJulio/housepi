library(testthat)

aux_add <- function(self, name = NA_character_,  value = NULL,
                    status = NA_integer_, const = NA, discr = FALSE,
                    func = NULL, func_args = NULL,
                    prior = NA_character_, ...,
                    min = NULL, max = NULL,
                    proposal = NA_character_) {
  self$add_quo(name = name, value = rlang::enquo(value),
               status = status, const = const, discr = discr,
               func = rlang::enquo(func),
               func_args = rlang::enquos(func_args),
               prior = prior, rlang::enquos(...),
               min = rlang::enquo(min), max = rlang::enquo(max),
               proposal = proposal)
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
  aux_add(augm, name = "test", proposal = "prior", value = NA) |>
    expect_error("prior.*test")
  aux_add(augm, name = "test", value = NA) |>
    expect_error("Missing value")

  augm <- HousepiAugmented$new()
  aux_add(augm, name = "test", value = 3, const = TRUE, proposal = "prior")
  expect_true(is.na(augm$defs$proposal))
})

test_that("Augmented variable definitions are correct", {
  # basic constant definition
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", value = 3)
  expect_equal(augm$defs$name, "dummy")
  # status is stored raw (character) to carry labels; NA when generic
  expect_equal(augm$defs$status, NA_character_)
  expect_equal(augm$defs$const, TRUE)
  expect_equal(augm$defs$discr, FALSE)
  expect_equal(augm$defs$value, "3|N")
  expect_equal(augm$defs$proposal, NA_character_)

  # variable definition
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", value = expo@D, prior = "uniform",
          min = expo_min@A, max = expo_max@P)
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
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", prior = "norm")

  expect_equal(augm$defs$const, FALSE)
  expect_equal(augm$defs$proposal, "rw")
})

test_that("Augmented independence proposal folds its arguments", {
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", prior = "uniform", min = 0, max = 5,
          proposal = "normal", proposal_mean = 0, proposal_sd = 1)
  expect_equal(augm$defs$proposal, "normal|mean:sd|{0|N}:{1|N}")

  # random-walk scale is not accepted by an independence sampler
  augm <- HousepiAugmented$new()
  aux_add(augm, name = "dummy", prior = "uniform", min = 0, max = 5,
          proposal = "normal", proposal_scale = 0.5) |>
    expect_error("random walk")
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
    some_func = list(func = "gamma",
                     func_args = c(shape = ~ shape, scale = ~ scale))
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
  helper_test_aux_load(augm, x = list(prior = "normal", obs = "normal"),
                       .kind = "obs") |>
    expect_error("only one between")
  # observation nodes are data-bound, never sampled
  helper_test_aux_load(augm, x = list(obs = "normal", proposal = "prior"),
                       .kind = "obs") |>
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
  expect_setequal(augm$defs[["name"]],
                  c("incub_period", "infec_time", "onset_date"))

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
    add_observation("x", obs = "uniform", loc = z,
                    min = - !!window, max = !!window) |>
    run_inference(steps = 10L, adapt = 0L, chains = 40L, progress = FALSE) |>
    suppressMessages() |>
    suppressWarnings() |>
    expect_no_error()

  # datum x = 5 is unreachable
  # UPDATE: it's the same as above: new-dataset decorator?
  data.frame(house_id = "outer_heaven", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_augmented("z", prior = "uniform") |>
    add_observation("x", obs = "uniform", loc = z,
                    min = - !!window, max = !!window) |>
    run_inference() |>
    suppressMessages() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_observation("x", obs = "uniform", min = - theta@P, max = theta@P) |>
    run_inference() |>
    suppressMessages() |>
    expect_error("free parameter dependencies")

  # generic error
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", value = 0.3) |>
    add_augmented("z", prior = "uniform") |>
    add_observation("x", obs = "uniform", min = - theta@P, max = theta@P) |>
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

#####################################
## AUGMENTED : END-TO-END MACHINERY ##
#####################################

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
      run_inference(steps = 5, thin = 1, adapt = 2,
                    chain = 1, threads = 1,
                    id_draws = hhmodel$data_wide$id),
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
  df <- data.frame(status = 1, exposure_max = exposure_max,
                   onset_min = floor(onset), onset_max = floor(onset))
  df$onset_mid <- (df$onset_min + df$onset_max) / 2 + 0.5

  results <- expect_no_error(
    df |>
      new_housepi_model("test") |>
      add_parameter(
        mean = list(prior = "uniform", min = 0.1, max = 5, proposal = "rw"),
        var = list(prior = "uniform", min = 0.1, max = 5, proposal = "rw")
      ) |>
      add_augmented(
        .infec_time = list(prior = "uniform", max = exposure_max@D + 0.99,
                           proposal = "rw"),
        incub_period = list(prior = "gamma", mean = mean@P, var = var@P,
                            min = 0.5, max = 14, proposal = "rw"),
        .onset_time = .infec_time + incub_period
      ) |>
      add_observation("onset_mid", obs = "normal",
                      mean = .onset_time, sd = 0.3) |>
      run_inference(steps = 20, adapt = 0, chains = 1, thin = 1,
                    progress = FALSE) |>
      suppressMessages() |>
      suppressWarnings()
  )

  # the gamma incubation prior contributes to the (finite) augmented log-prior
  lprior <- get_complete_draws(results)$lprior
  expect_true(all(is.finite(lprior)))
})
