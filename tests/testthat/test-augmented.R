library(testthat)

################################################################################
###################################  HELPERS  ##################################
################################################################################

# model carrying the data columns the definition blocks refer to
mk_augm <- local({
  hhmodel <- data.frame(house_id = "h", status = 1, expo = 1,
                        onset_date = 3, test_date = 2) |>
    new_housepi_model("test")
  list(model = hhmodel)
})

is_whole <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  length(x) > 0 && all(x == floor(x))
}

################################################################################
####################################  TESTS  ###################################
################################################################################

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
    add_augmented(import = list(list(model = "normal"))) |>
    expect_error("non-empty")
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
    add_augmented("test", proposal = "forward", value = NA) |>
    expect_error("[Nn]o .model. to sample from")
  mk_augm$model |>
    add_augmented("test", value = NA) |>
    expect_error("[Mm]issing value")

  defs <- mk_augm$model |>
    add_augmented("test", value = 3, const = TRUE, proposal = "forward") |>
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
    add_observation(onset_date = list(model = "normal", mean = ocelot)) |>
    expect_error("Did you mean .ocelot@P.")

  # nothing of that name anywhere
  mk_augm$model |>
    add_observation(
      onset_date = list(model = "normal", mean = liquid@A, sd = 1)
    ) |>
    expect_error("add_augmented")

  mk_augm$model |>
    add_augmented(
      incub = list(model = "gamma", shape = 2, scale = 1, status = 2),
      otacon = list(model = "normal", loc = incub, sd = 1)
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

test_that("The status process is capped by the engine status storage", {
  weights <- function(k) {
    stats::setNames(as.list(rep(1 / k, k)), paste0("unit_", seq_len(k)))
  }
  do.call(add_status, c(list(mk_augm$model), weights(254))) |>
    expect_silent()
  do.call(add_status, c(list(mk_augm$model), weights(255))) |>
    expect_error("255 statuses.*maximum of 254")
})

test_that("Augmented variable definitions are correct", {
  # basic constant definition
  augm <- mk_augm$model |>
    add_augmented(name = "dummy", value = 3) |>
    get_augmented()
  expect_equal(augm$name, "dummy")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$const, TRUE)
  expect_true(is.na(augm$discr))
  expect_equal(augm$value, 3)
  expect_equal(augm$proposal, NA_real_)

  # variable definition: each source resolves against its own family, so the
  # model only accepts it once 'expo_max' is a parameter and 'expo_min' an
  # augmented variable ('expo' being a data column)
  augm <- mk_augm$model |>
    add_parameter("expo_max", prior = "uniform", min = 0, max = 9) |>
    add_augmented("expo_min", model = "uniform", min = 0, max = 1) |>
    add_augmented(
      name = "dummy", value = expo@D, model = "uniform",
      min = expo_min@A, max = expo_max@P
    ) |>
    get_augmented() |>
    subset(name == "dummy")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$value, "expo|D")
  expect_equal(augm$model, "uniform")
  expect_equal(augm$min, "expo_min|A")
  expect_equal(augm$max, "expo_max|P")
})

test_that("Augmented sampling default is correct", {
  augm <- mk_augm$model |>
    add_augmented(name = "dummy", model = "normal") |>
    get_augmented()

  expect_equal(augm$const, FALSE)
  expect_equal(augm$proposal, "rw")

  # a misspelled family is caught at definition
  mk_augm$model |>
    add_augmented(name = "dummy", model = "norm") |>
    expect_error("[Uu]nknown distribution")
})

test_that("Augmented independence proposal folds its arguments", {
  augm <- mk_augm$model |>
    add_augmented(
      name = "dummy", model = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_mean = 0, proposal_sd = 1
    ) |>
    get_augmented()
  expect_equal(augm$proposal, "normal")
  expect_equal(augm$proposal_args[[1]], c(mean = 0, sd = 1))

  # random-walk scale is not accepted by an independence sampler
  mk_augm$model |>
    add_augmented(
      name = "dummy", model = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_scale = 0.5
    ) |>
    expect_error("random.?walk")
})

test_that("Independence proposal reject if unsupported prior domain", {
  hhmodel <- data.frame(house_id = "philanthropy", status = 1, x = 3) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_augmented("z", model = "uniform", min = theta@P, max = 5,
                  proposal = "gamma", proposal_shape = 2) |>
    evaluate_log_lik() |>
    expect_error("[Cc]annot cover")
})

test_that("Independence proposal inherits the model location shift", {
  hhmodel <- data.frame(id = c("tommy", "gina"), house_id = "house",
                        status = 1, min = c(4, 8)) |>
    new_housepi_model("test") |>
    add_augmented("test_date", model = "uniform", loc = min@D,
                  proposal = "gamma", proposal_shape = 2)

  results <- hhmodel |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L) |>
    suppressMessages() |>
    suppressWarnings()

  draws <- as.data.frame(get_augmented_draws(results))
  expect_true(all(draws$tommy_test_date >= 4 & draws$tommy_test_date <= 5))
  expect_true(all(draws$gina_test_date >= 8 & draws$gina_test_date <= 9))
})

test_that("augmented draws stay inside their live support", {
  df <- data.frame(id = seq(2), house_id = "shadow_moses", status = 1)

  # other augmented variable
  results <- df |>
    new_housepi_model("test") |>
    add_augmented(otacon = list(model = "uniform"),
                  test_date = list(model = "uniform", loc = otacon@A)) |>
    run_inference(steps = 400L, adapt = 0L, thin = 1L, chains = 1L) |>
    suppressMessages() |>
    suppressWarnings()
  d <- as.data.frame(get_augmented_draws(results))
  gap <- d[["1_test_date"]] - d[["1_otacon"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # transitive deterministic augmented variable
  results <- df |>
    new_housepi_model("test") |>
    add_augmented(jaeger = list(model = "uniform"),
                  gray_fox = jaeger,
                  test_date = list(model = "uniform", loc = gray_fox@A)) |>
    run_inference(steps = 500L, adapt = 0L, thin = 1L, chains = 1L) |>
    suppressMessages() |>
    suppressWarnings()
  draws <- as.data.frame(get_augmented_draws(results))
  gap <- draws[["1_test_date"]] - draws[["1_jaeger"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # parameter
  results <- df |>
    new_housepi_model("test") |>
    add_parameter("meryl", prior = "uniform") |>
    add_augmented("eva", model = "uniform", min = meryl@P, max = 5) |>
    run_inference(steps = 500L, adapt = 0L, thin = 1L, chains = 1L) |>
    suppressMessages() |>
    suppressWarnings()

  eva <- as.data.frame(get_augmented_draws(results))[["1_eva"]]
  meryl <- as.data.frame(get_parameter_draws(results))[["meryl"]]
  expect_true(all(eva >= meryl))
})

test_that("Augmented function definition is correct", {
  base <- mk_augm$model |>
    add_augmented(onset = list(model = "normal", mean = 0, sd = 1),
                  incub = list(model = "gamma", shape = 2, scale = 1),
                  infec = list(model = "normal", mean = 0, sd = 1))

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
    add_augmented(infec = list(model = "normal", mean = 0, sd = 1))

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
    add_augmented(infec = list(model = "normal", mean = 0, sd = 1))
  external_value <- 33
  augm <- base |>
    add_augmented(
      A = !!external_value,
      B = log(infec@A + !!external_value),
      C = c(
        model = "uniform",
        min = !!external_value, max = log(infec@A + !!external_value)
      ),
      D = list(
        model = "uniform",
        min = !!external_value, max = log(infec@A + !!external_value)
      )
    ) |>
    get_augmented()

  expect_equal(subset(augm, name == "A")$value, external_value)
  expect_equal(subset(augm, name == "B")$func, "log")
  expect_equal(
    subset(augm, name == "B")$func_args[[1]],
    c(x = sprintf("sum|x:x2|[infec|A]:[%s|N]", external_value))
  )
  # the two spellings of the same bounded definition agree
  cc <- subset(augm, name == "C")
  dd <- subset(augm, name == "D")
  cc$name <- dd$name <- NULL
  expect_equal(cc, dd, ignore_attr = TRUE)
})

test_that("Symbolic definitions are equivalent", {
  base <- mk_augm$model |>
    add_augmented(indiv_shape = list(model = "uniform", min = 0, max = 2))
  val <- 1
  aux <- base |>
    add_augmented(
      var1 = list(model = "normal", mean = indiv_shape + 1),
      var2 = list(model = "normal", mean = indiv_shape + !!val),
      var3 = list(model = "normal", mean = ~ indiv_shape + 1),
      var4 = list(model = "normal", mean = ~ indiv_shape + !!val)
    ) |>
    get_augmented() |>
    subset(name != "indiv_shape")
  aux$name <- NULL
  # check that all rows are identical
  expect_equal(nrow(unique(aux)), 1)

  # an external list is loaded through 'import', not as a bare definition
  external <- list(
    some_func = list(model = "gamma", shape = ~shape, scale = ~scale)
  )
  shapes <- mk_augm$model |>
    add_augmented(
      shape = list(model = "uniform", min = 1, max = 2),
      scale = list(model = "uniform", min = 1, max = 2)
    )
  shapes |>
    add_augmented(external, append = FALSE) |>
    expect_error("non.?empty")

  augm <- shapes |>
    add_augmented(import = external) |>
    get_augmented() |>
    subset(name == "some_func")
  expect_equal(augm$model, "gamma")
  expect_equal(augm$args[[1]], c(shape = "shape|A", scale = "scale|A"))
})

test_that("Observation definitions are encoded", {
  obs <- mk_augm$model |>
    add_augmented(onset_time = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(
      onset_date = list(model = "normal", mean = onset_time@A, sd = 1)
    ) |>
    get_observations()

  # marked as observation, data-bound and constant (never sampled)
  expect_identical(obs[["level"]], "member")
  expect_true(obs[["const"]])
  expect_identical(obs[["value"]], "onset_date|D")
  expect_true(is.na(obs[["proposal"]]))
  # distribution carries the observation arguments
  expect_identical(obs[["model"]], "normal")
  expect_equal(obs[["args"]][[1]], c(mean = "onset_time|A", sd = "1"))
})

test_that("Observation node rejects conflicting options", {
  base <- mk_augm$model
  # observation nodes are data-bound, never sampled
  base |>
    add_observation(x = list(model = "normal", proposal = "forward")) |>
    expect_warning("[Ii]gnored.*proposal")
})

test_that("Deterministic observation node is encoded as a function obs", {
  obs <- mk_augm$model |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # deterministic observation (no 'obs')
    add_observation(onset_date = floor(infec_time + incub_period)) |>
    get_observations()

  # tagged as observation, carrying a deterministic function and no distribution
  expect_identical(obs[["level"]], "member")
  expect_true(is.na(obs[["model"]]))
  expect_true(is.na(obs[["value"]]))
  expect_match(obs[["func"]], "floor")
})

test_that("Observation and augmented loads each own only their definitions", {
  hhmodel <- mk_augm$model |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # non-append observation definition keeps augmented definitions
    add_observation(
      onset_date = list(model = "normal", mean = incub_period@A, sd = 1),
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
      incub_period = list(model = "gamma", shape = 2, scale = 1),
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
    name = "first_test", model = "normal", mean = 0, sd = 1, level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    first_test = list(model = "normal", mean = 0, sd = 1), level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    import = list(first_test = list(model = "normal", mean = 0, sd = 1)),
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
    add_augmented(incub_period = list(model = "gamma", shape = 2, scale = 1)) |>
    add_observation(
      onset_date = list(model = "normal", mean = incub_period@A, sd = 1)
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
    add_augmented("z", model = "uniform") |>
    add_observation("x", model = "uniform", loc = z,
                    min = -!!window, max = !!window) |>
    run_inference(steps = 10L, adapt = 0L, chains = 40L) |>
    suppressMessages() |>
    suppressWarnings() |>
    expect_no_error()

  # datum x = 5 is unreachable
  # UPDATE: it's the same as above: new-dataset decorator?
  data.frame(house_id = "outer_heaven", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_augmented("z", model = "uniform") |>
    add_observation("x", model = "uniform", loc = z,
                    min = -!!window, max = !!window) |>
    run_inference() |>
    suppressMessages() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_observation("x", model = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    suppressMessages() |>
    expect_error("[Ff]ree parameter dependencies")

  # generic error
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", value = 0.3) |>
    add_augmented("z", model = "uniform") |>
    add_observation("x", model = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    expect_error("[Oo]utside of respective domain")
})

test_that("Reading an observation in models and loops are rejected", {

  # an augmented model must not read an observation (illegal in both
  # modes; caught at eager build, so the construction itself errors)
  data.frame(house_id = "outer_heaven", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", model = "uniform", min = 0, max = 1) |>
    add_augmented("latent", model = "normal", mean = onset@D) |>
    expect_error("[Mm]odel.*reads.*observation")

  # a deterministic read that closes a loop (the observation's support
  # depends on the reader) is circular
  data.frame(house_id = "metal_gear", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_augmented("ocelot", value = onset@D) |>
    add_observation("onset", model = "uniform", min = ocelot@A, max = 10) |>
    expect_error("[Cc]yclic")
})

test_that("simulation override (value + model) rejects forward dependencies", {
  # reads an observation
  data.frame(house_id = "shadow_moses", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", model = "uniform", min = 0, max = 1) |>
    # forward dependency: observation before augmented
    add_augmented("liquid", model = "normal", value = onset@D) |>
    expect_error("[Ss]imulation override.*observation")

  # reads another augmented (no observation present: guard still fires)
  data.frame(house_id = "outer_heaven", status = 1) |>
    new_housepi_model("test") |>
    add_augmented(liquid = list(model = "normal", value = solid@A),
                  solid = list(model = "gamma", shape = 2, scale = 1)) |>
    expect_error("[Ss]imulation override.*augmented")
})

test_that("degenerate reduction at member-level is cast to house-level", {
  # member-level reduction by mistake
  expect_warning(
    promoted <- mk_augm$model |>
      add_observation(first_test = min(test_date@D)) |>
      get_observations(),
    regexp = "no-effect.*reduction.*min"
  )
  expect_identical(promoted[["level"]], "house")
  # explicit house-level reduction
  mk_augm$model |>
    add_observation(first_test = min(test_date@D), level = "house") |>
    get_observations() |>
    expect_identical(promoted)

  # redundant reduction raise a error and are not cast
  mk_augm$model |>
    add_observation(first_test = min(test_date@D), level = "house") |>
    add_observation(codec = min(first_test@D)) |>
    expect_error("redundant.*codec")

  # genuine per-member expressions pass normally
  kept <- mk_augm$model |>
    add_observation(liquid = min(test_date@D, 0),
                    snake  = test_date@D,
                    ocelot = log(test_date@D)) |>
    get_observations() |>
    getElement("level") |>
    expect_identical(rep("member", 3))
})

################################################################################
#############################  END-TO-END MACHINERY  ###########################
################################################################################

test_that("Augmented aliasing works as expected", {
  # check that model arguments are correctly integrated
  simulated <- data.frame(id = seq(1, 1000),
                          house_id = rep(seq(10), each = 100),
                          status = 1) |>
    new_housepi_model("test") |>
    add_augmented(
      hour = list(model = "uniform", min = day, max = day + 1),
      day = list(model = "gamma", shape = 2, scale = 1)
    ) |>
    run_simulation()
  expect_true(all(
    simulated$day <= simulated$hour & simulated$hour <= simulated$day + 1
  ))
})

test_that("value paired with model overrides the draw in simulations only", {
  hhmodel <- data.frame(id = seq_len(100), house_id = 1,
                        status = 1, anchor = 7) |>
    new_housepi_model("test") |>
    add_augmented("eli", model = "uniform", max = 100, value = anchor@D)

  # simulation: pins the sampled variable to user value
  expect_true(all(run_simulation(hhmodel)$eli == 7))

  # inference: the value does not pin the sampler and draws span the prior
  draws <- hhmodel |>
    run_inference(steps = 40L, chains = 1L) |>
    get_augmented_draws() |>
    suppressWarnings()
  expect_false(all(draws[["1_eli"]] == 7))
})

test_that("load augmented list works fine for external lists", {
  hhmodel <- data.frame(id = "id1", house_id = "h1",
                        status = 1, date = as.Date("2021-01-01")) |>
    new_housepi_model("test")

  # an external list of definitions can be imported
  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  tmp <- add_augmented(hhmodel, import = external_list)$augmented

  hhmodel <- add_augmented(hhmodel, z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})

test_that("augmented template import loads", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_augmented(pupa = start_date@D + 1)

  # appending keeps the pre-existing definitions beside the template ones
  augm <- hhmodel |>
    add_augmented(import = "template") |>
    get_augmented()
  augm$name |>
    expect_setequal(c("pupa", "infec_time"))

  # dropping resets the store first
  augm <- hhmodel |>
    add_augmented(import = "template", append = FALSE) |>
    get_augmented()
  augm$name |>
    expect_setequal("infec_time")
})

test_that("Augmented traces export only relevant variables", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(alpha = 1, beta = list(prior = "uniform")) |>
    add_augmented(
      # constants
      blocked = 42,
      # private variables
      .hidden = list(model = "uniform"),
      .secret = list(model = "normal"),
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
    add_augmented("d_horse", model = "uniform", min = 5, max = 5)
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
      norm = list(model = "normal"),
      u1 = list(model = "uniform"),
      u2 = list(model = "uniform"),
      z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2)
    )
  expect_snapshot(print(hhmodel))

  results <- hhmodel |>
    run_inference(steps = 1e5, chains = 1, adapt = 0, thin = 1)
  # no parameters nor priors, only augmented
  expect_equal(posterior::nvariables(results$parameter_draws), 0L)
  posterior::nvariables(results$augmented_draws) |>
    expect_equal(dim(hhmodel$augmented)[1])

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
      .layer = list(model = "uniform", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      .dummy = list(model = "uniform"),
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
      .layer = list(model = "normal", status = 1),
      .layer = list(value = .dummy, status = 2),
      #> NO .layer provided for status 3 and 4
      .dummy = list(model = "lognormal", status = 2),
      z = .layer
    ) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")

  # second generation of definitions
  hhmodel |>
    add_augmented(
      .layer = list(model = "normal", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      #> NO .dummy provided for status 3, 4
      .dummy = list(model = "lognormal", status = c(2, 3)),
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
      z = list(status = 1, model = "uniform", min = 0, max = 1),
      z = list(status = 2, model = "uniform", min = 1, max = 2),
      z = list(status = 3, model = "uniform", min = 2, max = 3),
      z = list(status = 4, model = "uniform", min = 3, max = 4)
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
                    chain = 1, id_draws = hhmodel$data_wide$id),
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
    add_augmented(name = "z", model = "uniform")
  evaluate_log_lik(hhmodel)
  result <- hhmodel$export_augmented_data()
  expect_true(is.data.frame(result))
  expect_true("z" %in% colnames(result))
})

test_that("Incubation augmentation recast runs on the test model", {
  # a minimal analogue of the data-augmentation vignette: the incubation gamma
  # enters as an augmented prior, the reported onset day as an observation node.
  n <- 30
  exposure_max <- sample(0:3, n, replace = TRUE)
  infec <- runif(n, 0, exposure_max + 0.99)
  onset <- infec + rgamma(n, shape = 2, scale = 1.1)
  df <- data.frame(
    status = 1, exposure_max = exposure_max,
    onset_min = floor(onset), onset_max = floor(onset)
  )
  df$onset_mid <- 0.5 * (df$onset_min + df$onset_max + 1.)

  results <- df |>
    new_housepi_model("test") |>
    add_parameter(mean = list(prior = "uniform", min = 0.1, max = 5),
                  sd = list(prior = "uniform", min = 0.1, max = 5)) |>
    add_augmented(
      .infec_time = list(model = "uniform", max = exposure_max@D + 0.99),
      incub_period = list(
        model = "gamma", mean = mean@P, sd = sd@P, min = 0.5, max = 14
      ),
      .onset_time = .infec_time + incub_period
    ) |>
    add_observation("onset_mid", model = "normal",
                    mean = .onset_time, sd = 0.3) |>
    run_inference(steps = 20, adapt = 0, chains = 1) |>
    suppressMessages() |>
    suppressWarnings()

  # the gamma incubation prior contributes to the (finite) augmented log-prior
  lprior <- get_complete_draws(results)$lprior
  expect_true(all(is.finite(lprior)))
})

test_that("A surrogate observation overwrites the corresponding data column", {
  hhmodel <- data.frame(house_id = "groznyj", status = c(1, 0, 0),
                        test_date = c(2, 6, NA), start_date = rep(99, 3)) |>
    new_housepi_model("test") |>
    add_observation(
      start_date = list(value = min(test_date@D), level = "house"),
      w = test_date@D + 1
    ) |>
    add_augmented(z = list(value = w@D * 2))

  # the data view reports the derived value before evaluation
  expect_identical(unique(get_data_wide(hhmodel)[["start_date"]]), 2)
  expect_false("w" %in% colnames(get_data_wide(hhmodel)))

  # and after evaluation
  expect_no_error(evaluate_log_lik(hhmodel))
  expect_identical(unique(get_data_wide(hhmodel)[["start_date"]]), 2)
  expect_false("w" %in% colnames(get_data_wide(hhmodel)))

  # simulated observations are for all members, augmented just for infected
  simulated <- run_simulation(hhmodel)
  expect_equal(simulated[["w"]], c(3, 7, NA))
  expect_equal(simulated[["z"]], c(6, NA, NA))

  # FIXME: the requirement reported state "augmented", but it's "observation"
  # a stochastic observation carries no column to read
  #  datum is reported missing, rather than failing on the key it stands for
  data.frame(house_id = "rex", status = 1, otacon = 5) |>
    new_housepi_model("test") |>
    add_augmented(y = list(model = "uniform")) |>
    add_observation(hal = list(model = "normal", mean = y, sd = 1)) |>
    add_observation(otacon = hal@D + 2) |>
    evaluate_log_lik() |>
    expect_error("missing data column")
})

test_that("surrogate columns of both levels build in dependency order", {
  df <- data.frame(house_id = "rassvet", status = c(1, 0, 0),
                   test_date = c(2, 9, 4), window = 0, last_day = 0)
  hhmodel <- df |>
    new_housepi_model("test") |>
    add_observation(window = test_date@D + 3,
                    last_day = list(value = max(window@D), level = "house"))

  # every member carries a value regardless of status
  wide <- get_data_wide(hhmodel)
  expect_identical(wide[["window"]], c(5, 12, 7))
  expect_identical(unique(wide[["last_day"]]), 12)

  # pure surrogate observations with no raw data columns are not exported
  obs_names <- c("window", "last_day")
  df[obs_names] <- NULL
  out_data <- hhmodel |>
    add_data_wide(df) |>
    get_data_wide()
  expect_false(any(colnames(out_data) %in% obs_names))
})

test_that("observation cyclic dependency are captured", {
  data.frame(house_id = "rex", status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_augmented(z = list(model = "lognormal"),
                  y = list(model = "uniform", max = first@D)) |>
    add_observation(t_obs = list(model = "normal", mean = z, sd = 1)) |>
    add_observation("first", value = min(t_obs@D), level = "house") |>
    run_simulation() |>
    expect_error("cyclic dependency")
})

test_that("a surrogate column reading augmented is simulation only", {
  hhmodel <- data.frame(id = seq_len(100), house_id = "deck",
                        status = 1, ReD = 0) |>
    new_housepi_model("test") |>
    add_augmented(z = list(model = "uniform")) |>
    add_observation(ReD = z + 2)

  # the data view has no derived value to report: it keeps the imported datum
  expect_true(all(get_data_wide(hhmodel)[["ReD"]] == 0))

  # observation depends on sampled augmented: not a valid surrogate definition
  evaluate_log_lik(hhmodel) |>
    expect_error("must resolve from data alone")

  # simulation produces it, following the augmented variable as it is drawn
  simulated <- run_simulation(hhmodel)
  expect_equal(simulated[["ReD"]], simulated[["z"]] + 2)
})

test_that("Unspecified discretization is held for resolution", {
  obs <- mk_augm$model |>
    add_observation(test_date = list(model = "normal"),
                    onset_date = list(model = "normal", discr = TRUE),
                    expo = list(model = "normal", discr = FALSE)) |>
    get_observations()

  expect_type(obs$discr, "logical")
  expect_equal(obs$discr, c(NA, TRUE, FALSE))
})

test_that("An observation infers discretization from the column it reads", {
  hhmodel <- data.frame(house_id = rep(c("shadow", "moses"), each = 20),
                        date_col = as.Date("2005-02-28") + rep(0:1, 20),
                        int_col = rep(c(3L, 5L), 20),
                        whole_col = rep(c(3, 5), 20),
                        kept_col = rep(c(3L, 5L), 20),
                        status = 1) |>
    new_housepi_model("test") |>
    add_augmented(z = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      date_col = list(model = "normal", mean = z, sd = 2),
      int_col = list(model = "normal", mean = z, sd = 2),
      # a double column is always continuous
      whole_col = list(model = "normal", mean = z, sd = 2),
      # an explicit declaration is never overridden by the data
      kept_col = list(model = "normal", discr = FALSE, mean = z, sd = 2)
    )

  simulated <- run_simulation(hhmodel)

  # dates and integers are discrete, doubles are not
  expect_true(is_whole(simulated[["date_col"]]))
  expect_true(is_whole(simulated[["int_col"]]))
  expect_false(is_whole(simulated[["whole_col"]]))
  expect_false(is_whole(simulated[["kept_col"]]))

  # a deterministic observation is not influenced
  hhmodel |>
    add_observation(derived = z + 1) |>
    run_simulation() |>
    getElement("derived") |>
    is_whole() |>
    expect_false()
})

test_that("Discretization is resolved every time the data changes", {
  df <- data.frame(house_id = rep(c("naked", "snake"), 20), status = 1)

  hhmodel <- new_housepi_model(model = "test") |>
    add_augmented(z = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(tsuchinoko = list(model = "normal", mean = z, sd = 2))

  # the simulated observation resolves to continuos
  hhmodel |>
    add_data_wide(dplyr::mutate(df, tsuchinoko = rep(c(2.4, 5.1), 20))) |>
    run_simulation() |>
    getElement("tsuchinoko") |>
    is_whole() |>
    expect_false()

  # attaching an integer column changes the resolution to integer
  hhmodel |>
    add_data_wide(dplyr::mutate(df, tsuchinoko = rep(c(2L, 5L), 20))) |>
    run_simulation() |>
    getElement("tsuchinoko") |>
    is_whole() |>
    expect_true()
})
