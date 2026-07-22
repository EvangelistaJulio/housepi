library(testthat)

# model carrying the data columns the definition blocks refer to
mk_augm <- local({
  hhmodel <- data.frame(house_id = "h", status = 1, expo = 1,
                        onset_date = 3, test_date = 2) |>
    new_housepi_model("test")
  list(model = hhmodel)
})

test_that("Augmented definition error are recognized", {
  mk_augm$model |>
    add_augmented("infec@") |>
    expect_error("[Ii]nvalid.*characters")
  mk_augm$model |>
    add_augmented("infec") |>
    expect_error("[Ii]ncomplete definition")
  mk_augm$model |>
    add_augmented(name = "") |>
    expect_error("non.?empty")

  # UPDATE: the hint `'!!!'` here seems off
  mk_augm$model |>
    add_augmented(import = list(list(prior = "normal")))
})

test_that("Constant value definition are coherent", {
  augm <- mk_augm$model |> add_augmented(dummy = 3)
  tmp <- mk_augm$model |> add_augmented(name = "dummy", value = 3)
  expect_equal(get_augmented(augm), get_augmented(tmp))

  augm <- mk_augm$model |> add_augmented(dummy = 3@D)
  tmp <- mk_augm$model |> add_augmented(name = "dummy", value = 3@D)
  expect_equal(get_augmented(augm), get_augmented(tmp))
})

test_that("Undefined definitions are controoled", {
  mk_augm$model |>
    add_augmented("test", proposal = "prior", value = NA) |>
    expect_error("[Nn]o .prior. to sample from")
  mk_augm$model |>
    add_augmented("test", value = NA) |>
    expect_error("[Mm]issing value")

  defs <- mk_augm$model |>
    add_augmented("test", value = 3, const = TRUE, proposal = "prior") |>
    get_augmented()
  expect_true(is.na(defs$proposal))
})

test_that("Dependency errors name the family and the intended source", {
  # an unsuffixed name is read as augmented: the error says so and points at
  # the data column / parameter of the same name the user meant
  mk_augm$model |>
    add_observation(start_date = min(test_date)) |>
    expect_error("augmented.*definition.*test_date.")
  mk_augm$model |>
    add_observation(start_date = min(test_date)) |>
    expect_error("Did you mean .test_date@D.")
  mk_augm$model |>
    add_parameter(ocelot = list(prior = "gamma", shape = 2, scale = 1)) |>
    add_observation(onset_date = list(obs = "normal", mean = ocelot, sd = 1)) |>
    expect_error("Did you mean .ocelot@P.")

  # nothing of that name anywhere
  mk_augm$model |>
    add_observation(
      onset_date = list(obs = "normal", mean = liquid@A, sd = 1)
    ) |>
    expect_error("add_augmented")

  mk_augm$model |>
    add_augmented(
      incub = list(prior = "gamma", shape = 2, scale = 1, status = 2),
      otacon = list(prior = "normal", loc = incub, sd = 1)
    ) |>
    expect_error(".*incub.*status 1")
})

test_that("A status left without definitions names the status", {
  # every per-member definition is pinned to a status no individual carries
  mk_augm$model |>
    add_observation(start_date = list(func = min(test_date@D), status = 2)) |>
    expect_error("No.*augmented.*definition.*status 1")

  # no status block to fall back on: the generic definition is the missing one
  data.frame(house_id = "rex", status = 0, test_date = 2) |>
    new_housepi_model("test") |>
    add_observation(
      start_date = list(func = min(test_date@D), status = c(2, 3))
    ) |>
    expect_error("[Nn]o.*generic.*augmented.*definition")
})

test_that("Augmented variable definitions are correct", {
  # basic constant definition
  augm <- mk_augm$model |>
    add_augmented(name = "dummy", value = 3) |>
    get_augmented()
  expect_equal(augm$name, "dummy")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$const, TRUE)
  expect_equal(augm$discr, FALSE)
  expect_equal(augm$value, 3)
  expect_equal(augm$proposal, NA_real_)

  # variable definition: each source resolves against its own family, so the
  # model only accepts it once 'expo_max' is a parameter and 'expo_min' an
  # augmented variable ('expo' being a data column)
  augm <- mk_augm$model |>
    add_parameter("expo_max", prior = "uniform", min = 0, max = 9) |>
    add_augmented("expo_min", prior = "uniform", min = 0, max = 1) |>
    add_augmented(
      name = "dummy", value = expo@D, prior = "uniform",
      min = expo_min@A, max = expo_max@P
    ) |>
    get_augmented() |>
    subset(name == "dummy")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$value, "expo|D")
  expect_equal(augm$prior, "uniform")
  expect_equal(augm$min, "expo_min|A")
  expect_equal(augm$max, "expo_max|P")
})

test_that("Augmented sampling default is correct", {
  augm <- mk_augm$model |>
    add_augmented(name = "dummy", prior = "normal") |>
    get_augmented()

  expect_equal(augm$const, FALSE)
  expect_equal(augm$proposal, "rw")

  # a misspelled family is caught at definition
  mk_augm$model |>
    add_augmented(name = "dummy", prior = "norm") |>
    expect_error("[Uu]nknown distribution")
})

test_that("Augmented independence proposal folds its arguments", {
  augm <- mk_augm$model |>
    add_augmented(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_mean = 0, proposal_sd = 1
    ) |>
    get_augmented()
  expect_equal(augm$proposal, "normal")
  expect_equal(augm$proposal_args[[1]], c(mean = 0, sd = 1))

  # random-walk scale is not accepted by an independence sampler
  mk_augm$model |>
    add_augmented(
      name = "dummy", prior = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_scale = 0.5
    ) |>
    expect_error("random.?walk")
})

test_that("Independence proposal reject if unsupported prior domain", {
  hhmodel <- data.frame(house_id = "philanthropy", status = 1, x = 3) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_augmented("z", prior = "uniform", min = theta@P, max = 5,
                  proposal = "gamma", proposal_shape = 2) |>
    evaluate_log_lik() |>
    expect_error("[Cc]annot cover")
})

test_that("Independence proposal inherits the prior location shift", {
  # a proposal from another family must sample in the prior's absolute frame:
  # with loc = min@D the support is [min, min + 1], so every draw lands there.
  # id_draws = list(house_id =) traces the whole household: both infected members
  # of "house" are stored, exercising the household -> members expansion.
  hhmodel <- data.frame(
    id = c("tommy", "gina"), house_id = "house",
    status = 1, min = c(4, 8)
  ) |>
    new_housepi_model("test") |>
    add_augmented("test_date", prior = "uniform", loc = min@D,
                  proposal = "gamma", proposal_shape = 2)

  results <- hhmodel |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L,
                  threads = 1L, progress = FALSE,
                  id_draws = list(house_id = "house")) |>
    suppressMessages() |>
    suppressWarnings()

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
    add_augmented(otacon = list(prior = "uniform"),
                  test_date = list(prior = "uniform", loc = otacon@A)) |>
    run_inference(steps = 500, adapt = 100, thin = 1, chains = 1,
                  threads = 1, progress = FALSE, id_draws = "1") |>
    suppressMessages() |>
    suppressWarnings()
  d <- as.data.frame(get_augmented_draws(augm))
  gap <- d[["1_test_date"]] - d[["1_otacon"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # parameter -> augmented: z ~ Uniform(min = meryl@P, max = 5), support [meryl, 5],
  # so every draw must satisfy z >= meryl
  par <- df |>
    new_housepi_model("test") |>
    add_parameter("meryl", prior = "uniform") |>
    add_augmented("z", prior = "uniform", min = meryl@P, max = 5) |>
    run_inference(steps = 500L, adapt = 100L, thin = 1L, chains = 1L,
                  threads = 1L, progress = FALSE, id_draws = "1") |>
    suppressMessages() |>
    suppressWarnings()
  z <- as.data.frame(get_augmented_draws(par))[["1_z"]]
  meryl <- as.data.frame(get_parameter_draws(par))[["meryl"]]
  expect_true(all(z >= meryl))
})

test_that("Augmented function definition is correct", {
  base <- mk_augm$model |>
    add_augmented(
      onset = list(prior = "normal", mean = 0, sd = 1),
      incub = list(prior = "gamma", shape = 2, scale = 1),
      infec = list(prior = "normal", mean = 0, sd = 1)
    )

  # difference definition
  augm <- base |>
    add_augmented("dummy", func = onset - incub) |>
    get_augmented() |>
    subset(name == "dummy")
  expect_equal(augm$func, "sub")
  expect_equal(augm$func_args[[1]], c(x = "onset|A", x2 = "incub|A"))
  expect_equal(augm$value, NA_real_)

  # overwriting warnings
  base |>
    add_augmented("dummy", func = onset - incub) |>
    add_augmented("dummy", func = "log", func_args = c(x = infec@A)) |>
    expect_warning("[Oo]verwriting")
  base |>
    add_augmented("dummy", func = onset - incub) |>
    add_augmented("rate", value = 3) |>
    add_augmented("dummy", func = exp(-rate)) |>
    expect_warning("[Oo]verwriting")
})

test_that("Augmented load enable special definitions", {
  base <- mk_augm$model |>
    add_augmented(infec = list(prior = "normal", mean = 0, sd = 1))

  # alternative function call definitions
  augm <- base |>
    add_augmented("dummy", func = "log", func_args = c(x = infec@A))
  tmp <- base |> add_augmented(dummy = log(infec@A))
  expect_equal(get_augmented(tmp), get_augmented(augm))

  expect_warning(
    overwritten <- add_augmented(tmp, name = "dummy", func = log(infec@A)),
    "[Oo]verwriting"
  )
  expect_equal(get_augmented(overwritten), get_augmented(augm))

  # multi-status definition: one row per status, identical otherwise
  multi <- base |>
    add_augmented(dummy = list(func = log(infec@A), status = c(1, 3))) |>
    get_augmented() |>
    subset(name == "dummy")
  expect_equal(multi$status, c("1", "3"))
  multi$status <- NULL
  expect_equal(multi[1, , drop = TRUE], multi[2, , drop = TRUE])
})

test_that("Augmented support rlang double bang", {
  base <- mk_augm$model |>
    add_augmented(infec = list(prior = "normal", mean = 0, sd = 1))
  external_value <- 33
  augm <- base |>
    add_augmented(
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
    ) |>
    get_augmented()

  expect_equal(subset(augm, name == "A")$value, external_value)
  expect_equal(subset(augm, name == "B")$func, "log")
  expect_equal(
    subset(augm, name == "B")$func_args[[1]],
    c(x = sprintf("sum|x:x2|{infec|A}:{%s|N}", external_value))
  )
  # the two spellings of the same bounded definition agree
  cc <- subset(augm, name == "C")
  dd <- subset(augm, name == "D")
  cc$name <- dd$name <- NULL
  expect_equal(cc, dd, ignore_attr = TRUE)
})

test_that("Symbolic definitions are equivalent", {
  base <- mk_augm$model |>
    add_augmented(indiv_shape = list(prior = "uniform", min = 0, max = 2))
  val <- 1
  aux <- base |>
    add_augmented(
      var1 = list(prior = "normal", mean = indiv_shape + 1),
      var2 = list(prior = "normal", mean = indiv_shape + !!val),
      var3 = list(prior = "normal", mean = ~ indiv_shape + 1),
      var4 = list(prior = "normal", mean = ~ indiv_shape + !!val)
    ) |>
    get_augmented() |>
    subset(name != "indiv_shape")
  aux$name <- NULL
  # check that all rows are identical
  expect_equal(nrow(unique(aux)), 1)

  # an external list is loaded through 'import', not as a bare definition
  external <- list(
    some_func = list(prior = "gamma", shape = ~shape, scale = ~scale)
  )
  shapes <- mk_augm$model |>
    add_augmented(
      shape = list(prior = "uniform", min = 1, max = 2),
      scale = list(prior = "uniform", min = 1, max = 2)
    )
  shapes |>
    add_augmented(external, append = FALSE) |>
    expect_error("non.?empty")

  augm <- shapes |>
    add_augmented(import = external) |>
    get_augmented() |>
    subset(name == "some_func")
  expect_equal(augm$prior, "gamma")
  expect_equal(augm$args[[1]], c(shape = "shape|A", scale = "scale|A"))
})

test_that("Observation definitions are encoded", {
  obs <- mk_augm$model |>
    add_augmented(onset_time = list(prior = "normal", mean = 0, sd = 1)) |>
    add_observation(
      onset_date = list(obs = "normal", mean = onset_time@A, sd = 1)
    ) |>
    get_observations()

  # marked as observation, data-bound and constant (never sampled)
  expect_identical(obs[["level"]], "member")
  expect_true(obs[["const"]])
  expect_identical(obs[["value"]], "onset_date|D")
  expect_true(is.na(obs[["proposal"]]))
  # distribution carries the observation arguments
  expect_identical(obs[["prior"]], "normal")
  expect_equal(obs[["args"]][[1]], c(mean = "onset_time|A", sd = "1"))
})

test_that("Observation node rejects conflicting options", {
  base <- mk_augm$model
  # prior and obs are mutually exclusive
  base |>
    add_observation(x = list(prior = "normal", obs = "normal")) |>
    expect_error("[Oo]nly one between")
  # observation nodes are data-bound, never sampled
  base |>
    add_observation(x = list(obs = "normal", proposal = "prior")) |>
    expect_warning("[Ii]gnored.*proposal")
})

test_that("Deterministic observation node is encoded as a function obs", {
  obs <- mk_augm$model |>
    add_augmented(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # deterministic observation (no 'obs')
    add_observation(onset_date = floor(infec_time + incub_period)) |>
    get_observations()

  # tagged as observation, carrying a deterministic function and no distribution
  expect_identical(obs[["level"]], "member")
  expect_true(is.na(obs[["prior"]]))
  expect_true(is.na(obs[["value"]]))
  expect_match(obs[["func"]], "floor")
})

test_that("Observation and augmented loads each own only their definitions", {
  hhmodel <- mk_augm$model |>
    add_augmented(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # non-append observation definition keeps augmented definitions
    add_observation(
      onset_date = list(obs = "normal", mean = incub_period@A, sd = 1),
      append = FALSE
    )
  expect_setequal(
    get_augmented(hhmodel)[["name"]],
    c("incub_period", "infec_time")
  )
  expect_identical(get_observations(hhmodel)[["name"]], "onset_date")

  # symmetrically, a non-append augmented load keeps observations
  hhmodel <- hhmodel |>
    add_augmented(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      append = FALSE
    )
  expect_identical(get_augmented(hhmodel)[["name"]], "incub_period")
  expect_identical(get_observations(hhmodel)[["name"]], "onset_date")
})

test_that("Household observations reject a distribution", {
  hhmodel <- mk_augm$model
  err <- "cannot be sampled"

  # every input branch is guarded, the level is never silently downgraded
  add_observation(hhmodel,
    name = "first_test", obs = "normal", mean = 0, sd = 1, level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    first_test = list(obs = "normal", mean = 0, sd = 1), level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    import = list(first_test = list(obs = "normal", mean = 0, sd = 1)),
    level = "house"
  ) |>
    expect_error(err)
})

test_that("Level is set per definition, the argument being the bulk default", {
  # a per-definition 'level' picks the kind of its own entry, exactly as
  # 'status' does; regression: it used to be swallowed by add_quo's dots and
  # folded into the definition arguments, leaving every entry on the default
  mixed <- mk_augm$model |>
    add_observation(
      start_date = list(func = min(test_date@D), level = "house"),
      end_date = list(func = 21, const = TRUE, level = "member")
    ) |>
    get_observations()
  expect_identical(
    mixed[["level"]][mixed[["name"]] == "start_date"], "house"
  )
  expect_identical(
    mixed[["level"]][mixed[["name"]] == "end_date"], "member"
  )
  # 'level' never reaches the definition arguments
  expect_false("args" %in% names(mixed))

  # the argument is the default: entries not setting their own inherit it
  bulk <- mk_augm$model |>
    add_observation(
      start_date = list(func = min(test_date@D)),
      end_date = list(func = 21, const = TRUE, level = "member"),
      level = "house"
    ) |>
    get_observations()
  expect_identical(bulk[["level"]][bulk[["name"]] == "start_date"], "house")
  expect_identical(bulk[["level"]][bulk[["name"]] == "end_date"], "member")

  # imported definitions honour it too
  imported <- mk_augm$model |>
    add_observation(import = list(
      start_date = list(func = quote(min(test_date@D)), level = "house"),
      end_date = list(func = 21, const = TRUE)
    )) |>
    get_observations()
  expect_identical(
    imported[["level"]][imported[["name"]] == "start_date"], "house"
  )
  expect_identical(
    imported[["level"]][imported[["name"]] == "end_date"], "member"
  )

  # an unknown level is rejected per definition
  mk_augm$model |>
    add_observation(start_date = list(func = 1, level = "liquid")) |>
    expect_error("should be one of")

  # 'level' is meaningless outside observations
  mk_augm$model |>
    add_augmented(ocelot = list(value = 1, level = "house")) |>
    expect_error("reserved to observations")
})

test_that("Household observations do not interfer with augmented defs", {
  hhmodel <- mk_augm$model |>
    add_augmented(incub_period = list(prior = "gamma", shape = 2, scale = 1)) |>
    add_observation(
      onset_date = list(obs = "normal", mean = incub_period@A, sd = 1)
    )

  # non-append household observation keeps augmented definitions
  hhmodel <- add_observation(
    hhmodel, first_test = min(test_date@D), level = "house", append = FALSE
  )
  expect_identical(get_augmented(hhmodel)[["name"]], "incub_period")
  obs <- get_observations(hhmodel)

  # observations are a single family: the member row went with the reset
  expect_setequal(obs[["name"]], "first_test")
  expect_identical(obs[["level"]], "house")
  expect_identical(obs[["func"]], "min")
  expect_identical(obs[["func_args"]][[1]], c(x = "test_date|D"))

  # appending keeps both levels side by side
  hhmodel <- add_observation(
    hhmodel, last_test = max(test_date@D), level = "house"
  )
  obs <- get_observations(hhmodel)
  expect_setequal(obs[["name"]], c("first_test", "last_test"))
  expect_identical(get_augmented(hhmodel)[["name"]], "incub_period")
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
    expect_error("[Ff]ree parameter dependencies")

  # generic error
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", value = 0.3) |>
    add_augmented("z", prior = "uniform") |>
    add_observation("x", obs = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    expect_error("[Oo]utside of respective domain")
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
    expect_error("[Pp]rior.*reads.*observation")

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
                  value = onset@D) |>
    expect_error("[Ss]imulation override.*observation")

  # reads another augmented (no observation present: guard still fires)
  data.frame(house_id = "outer_heaven", status = 1) |>
    new_housepi_model("test") |>
    add_augmented(
      liquid = list(prior = "normal", mean = 0, sd = 1, value = solid@A),
      solid = list(prior = "gamma", shape = 2, scale = 1)
    ) |>
    expect_error("[Ss]imulation override.*augmented")
})

test_that("Group/contact risk keyed on an observation is rejected in simulation only", {
  # the rate feeds transmission, which the observation is generated from: fine in
  # inference (observed data), circular when generating during simulation.
  hhmodel <- data.frame(house_id = c(1, 1), status = c(1, 0), severity = c(2, 3)) |>
    new_housepi_model("final_size") |>
    add_group_risk("susceptibility", by = "severity") |>
    add_observation("severity", obs = "uniform", min = 0, max = 5) |>
    add_parameter_template()

  hhmodel |>
    run_simulation(nrep = 1L) |>
    suppressMessages() |>
    expect_error("[Gg]roup risk.*observation")
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
    add_augmented(snake = list(prior = "uniform", max = 100, value = anchor@D))

  # simulation: every drawn value collapses onto the override
  expect_true(all(run_simulation(hhmodel)$snake == 7))

  # inference: the value does not pin the sampler (draws span the prior, not the override)
  draws <- suppressWarnings( # short smoke run trips the convergence-size diagnostic
    hhmodel |>
      run_inference(steps = 40L, chains = 1L, id_draws = "1", progress = FALSE) |>
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
    expect_error("[Cc]yclic augmented definition")

  # triangular nest
  hhmodel |>
    add_augmented(z = .layer, .layer = .dummy, .dummy = z) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic augmented definition")

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
    expect_error("[Nn]o.*definition.*status")

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
    expect_error("[Nn]o.*definition.*status")
})

test_that("Missing augmented function arguments raise errors", {
  hhmodel <- data.frame(house_id = seq(3), id = seq(3), status = 1) |>
    new_housepi_model("test")

  # non default distribution arguments must be supplied
  hhmodel |>
    add_augmented(name = "z", func = "gamma") |>
    expect_error("[Mm]issing non.?default")
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
  df$onset_mid <- 0.5 * (df$onset_min + df$onset_max + 1.)

  results <- expect_no_error(
    df |>
      new_housepi_model("test") |>
      add_parameter(
        mean = list(prior = "uniform", min = 0.1, max = 5),
        var = list(prior = "uniform", min = 0.1, max = 5)
      ) |>
      add_augmented(
        .infec_time = list(
          prior = "uniform", max = exposure_max@D + 0.99
        ),
        incub_period = list(
          prior = "gamma", mean = mean@P, var = var@P,
          min = 0.5, max = 14
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
