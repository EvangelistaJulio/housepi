library(testthat)

################################################################################
###################################  HELPERS  ##################################
################################################################################

# model carrying the data columns the definition blocks refer to
mk_latent <- local({
  hhmodel <- data.frame(house_id = "h", status = 1, expo = 1,
                        onset_date = 3, test_date = 2) |>
    new_housepi_model("test")
  list(model = hhmodel)
})

on_integer_lattice <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  length(x) > 0 && all(x == floor(x))
}

################################################################################
####################################  TESTS  ###################################
################################################################################

test_that("Latent definition error are recognized", {
  mk_latent$model |>
    add_latent("infec@") |>
    expect_error("[Ii]nvalid.*characters")
  mk_latent$model |>
    add_latent("infec") |>
    expect_error("[Ii]ncomplete definition")
  mk_latent$model |>
    add_latent(name = "") |>
    expect_error("non.?empty")

  # UPDATE: the hint `'!!!'` here seems off
  mk_latent$model |>
    add_latent(import = list(list(model = "normal"))) |>
    expect_error("non-empty")
})

test_that("Constant value definition are coherent", {
  augm <- mk_latent$model |> add_latent(dummy = 3)
  tmp <- mk_latent$model |> add_latent(name = "dummy", value = 3)
  expect_equal(get_latent(augm), get_latent(tmp))

  augm <- mk_latent$model |> add_latent(dummy = 3@D)
  tmp <- mk_latent$model |> add_latent(name = "dummy", value = 3@D)
  expect_equal(get_latent(augm), get_latent(tmp))
})

test_that("Undefined definitions are controoled", {
  mk_latent$model |>
    add_latent("test", proposal = "forward", value = NA) |>
    expect_error("[Nn]o .model. to sample from")
  mk_latent$model |>
    add_latent("test", value = NA) |>
    expect_error("[Mm]issing value")

  defs <- mk_latent$model |>
    add_latent("test", value = 3, const = TRUE, proposal = "forward") |>
    get_latent()
  expect_true(is.na(defs$proposal))
})

test_that("Passing a prior raise a hint on model", {
  mk_latent$model |>
    add_latent(infec_time = list(prior = "uniform", min = -1, max = 2)) |>
    expect_error("[Ss]pare.*prior.*infec_time")
  mk_latent$model |>
    add_latent(infec_time = list(prior = "uniform", min = -1, max = 2)) |>
    expect_error("[Dd]id you mean.*model")
})

test_that("A func definition is deterministic", {
  # a plain number is a constant value, whichever argument carried it
  by_func <- mk_latent$model |> add_latent(dummy = list(func = 21))
  by_value <- mk_latent$model |> add_latent(dummy = list(value = 21))
  expect_equal(get_latent(by_func), get_latent(by_value))

  # deterministic definitions reject model or proposal
  mk_latent$model |>
    add_latent(dummy = list(func = "log", func_args = c(x = 2),
                            model = "normal", mean = 0, sd = 1)) |>
    expect_error("only one.*func.*model")
  mk_latent$model |>
    add_latent(dummy = list(func = "log", func_args = c(x = 2),
                            proposal = "rw")) |>
    expect_error("[Dd]eterministic.*dummy.*cannot be sampled")
})

test_that("Dependency errors name the family and the intended source", {
  # definitions are left unresolved until a wall
  mk_latent$model |>
    add_observation(start_date = min(test_date)) |>
    expect_silent()

  # an unsuffixed name is read as latent: the error says so and points at
  # the data column / parameter of the same name the user meant
  mk_latent$model |>
    add_observation(start_date = min(test_date)) |>
    evaluate_log_lik() |>
    expect_error("latent.*definition.*test_date.")
  mk_latent$model |>
    add_observation(start_date = min(test_date)) |>
    evaluate_log_lik() |>
    expect_error("[Dd]id you mean.*test_date@D")
  mk_latent$model |>
    add_parameter(ddog = list(prior = "gamma", shape = 2, scale = 1)) |>
    add_observation(onset_date = list(model = "normal", mean = ddog)) |>
    evaluate_log_lik() |>
    expect_error("[Dd]id you mean.*ddog@P")

  # nothing of that name anywhere
  mk_latent$model |>
    add_observation(onset_date = list(model = "normal", mean = liquid@A)) |>
    evaluate_log_lik() |>
    expect_error("add_latent")

  mk_latent$model |>
    add_latent(
      incub = list(model = "gamma", shape = 2, scale = 1, status = 2),
      otacon = list(model = "normal", loc = incub, sd = 1)
    ) |>
    evaluate_log_lik() |>
    expect_error(".*incub.*status 1")
})

test_that("A status left without definitions names the status", {
  # definitions are left unresolved until a wall
  mk_latent$model |>
    add_observation(start_date = list(func = min(test_date@D), status = 2)) |>
    add_observation(start_date = list(func = min(test_date@D), status = 1)) |>
    expect_silent()

  # every per-member definition is pinned to a status no individual carries
  mk_latent$model |>
    add_observation(start_date = list(func = min(test_date@D), status = 2)) |>
    evaluate_log_lik() |>
    expect_error("No.*latent.*definition.*status 1")

  # no status block to fall back on: the generic definition is the missing one
  data.frame(house_id = "rex", status = 0, test_date = 2) |>
    new_housepi_model("test") |>
    add_observation(
      start_date = list(func = min(test_date@D), status = c(2, 3))
    ) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*generic.*latent.*definition")
})

test_that("The status process is capped by the engine status storage", {
  weights <- function(k) {
    stats::setNames(as.list(rep(1 / k, k)), paste0("unit_", seq_len(k)))
  }
  do.call(add_status, c(list(mk_latent$model), weights(254))) |>
    expect_silent()
  do.call(add_status, c(list(mk_latent$model), weights(255))) |>
    expect_error("255 statuses.*maximum of 254")
})

test_that("Latent variable definitions are correct", {
  # basic constant definition
  augm <- mk_latent$model |>
    add_latent(name = "dummy", value = 3) |>
    get_latent()
  expect_equal(augm$name, "dummy")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$const, TRUE)
  expect_equal(augm$discr, FALSE)
  expect_equal(augm$value, 3)
  expect_equal(augm$proposal, NA_real_)

  # variable definition: each source resolves against its own family, so the
  # model only accepts it once 'expo_max' is a parameter and 'expo_min' an
  # latent variable ('expo' being a data column)
  augm <- mk_latent$model |>
    add_parameter("expo_max", prior = "uniform", min = 0, max = 9) |>
    add_latent("expo_min", model = "uniform", min = 0, max = 1) |>
    add_latent(
      name = "dummy", value = expo@D, model = "uniform",
      min = expo_min@A, max = expo_max@P
    ) |>
    get_latent() |>
    subset(name == "dummy")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$value, "expo|D")
  expect_equal(augm$model, "uniform")
  expect_equal(augm$min, "expo_min|A")
  expect_equal(augm$max, "expo_max|P")
})

test_that("Latent sampling default is correct", {
  augm <- mk_latent$model |>
    add_latent(name = "dummy", model = "normal") |>
    get_latent()

  expect_equal(augm$const, FALSE)
  expect_equal(augm$proposal, "rw")

  # a misspelled family is caught at definition
  mk_latent$model |>
    add_latent(name = "dummy", model = "norm") |>
    expect_error("[Uu]nknown distribution")
})

test_that("Latent independence proposal folds its arguments", {
  augm <- mk_latent$model |>
    add_latent(
      name = "dummy", model = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_mean = 0, proposal_sd = 1
    ) |>
    get_latent()
  expect_equal(augm$proposal, "normal")
  expect_equal(augm$proposal_args[[1]], c(mean = 0, sd = 1))

  # random-walk scale is not accepted by an independence sampler
  mk_latent$model |>
    add_latent(
      name = "dummy", model = "uniform", min = 0, max = 5,
      proposal = "normal", proposal_scale = 0.5
    ) |>
    expect_error("random.?walk")
})

test_that("Independence proposal reject if unsupported prior domain", {
  hhmodel <- data.frame(house_id = "philanthropy", status = 1, x = 3) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_latent("z", model = "uniform", min = theta@P, max = 5,
               proposal = "gamma", proposal_shape = 2) |>
    evaluate_log_lik() |>
    expect_error("[Cc]annot cover")
})

test_that("Independence proposal inherits the model location shift", {
  hhmodel <- data.frame(id = c("tommy", "gina"), house_id = "house",
                        status = 1, min = c(4, 8)) |>
    new_housepi_model("test") |>
    add_latent("test_date", model = "uniform", loc = min@D,
               proposal = "gamma", proposal_shape = 2)

  results <- hhmodel |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  draws <- as.data.frame(get_latent_draws(results))
  expect_true(all(draws$tommy_test_date >= 4 & draws$tommy_test_date <= 5))
  expect_true(all(draws$gina_test_date >= 8 & draws$gina_test_date <= 9))
})

test_that("latent draws stay inside their live support", {
  df <- data.frame(id = seq(2), house_id = "shadow_moses", status = 1)

  # other latent variable
  results <- df |>
    new_housepi_model("test") |>
    add_latent(otacon = list(model = "uniform"),
               test_date = list(model = "uniform", loc = otacon@A)) |>
    run_inference(steps = 400L, adapt = 0L, thin = 1L, chains = 1L) |>
    suppressWarnings()
  d <- as.data.frame(get_latent_draws(results))
  gap <- d[["1_test_date"]] - d[["1_otacon"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # transitive deterministic latent variable
  results <- df |>
    new_housepi_model("test") |>
    add_latent(jaeger = list(model = "uniform"),
               gray_fox = jaeger,
               test_date = list(model = "uniform", loc = gray_fox@A)) |>
    run_inference(steps = 400L, adapt = 0L, thin = 1L, chains = 1L) |>
    suppressWarnings()
  draws <- as.data.frame(get_latent_draws(results))
  gap <- draws[["1_test_date"]] - draws[["1_jaeger"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # parameter
  results <- df |>
    new_housepi_model("test") |>
    add_parameter("meryl", prior = "uniform") |>
    add_latent("eva", model = "uniform", min = meryl@P, max = 5) |>
    run_inference(steps = 400L, adapt = 0L, thin = 1L, chains = 1L) |>
    suppressWarnings()

  eva <- as.data.frame(get_latent_draws(results))[["1_eva"]]
  meryl <- as.data.frame(get_parameter_draws(results))[["meryl"]]
  expect_true(all(eva >= meryl))
})

test_that("Latent function definition is correct", {
  base <- mk_latent$model |>
    add_latent(onset = list(model = "normal", mean = 0, sd = 1),
               incub = list(model = "gamma", shape = 2, scale = 1),
               infec = list(model = "normal", mean = 0, sd = 1))

  # difference definition
  augm <- base |>
    add_latent("dummy", func = onset - incub) |>
    get_latent() |>
    subset(name == "dummy")
  expect_equal(augm$func, "sub")
  expect_equal(augm$func_args[[1]], c(x = "onset|A", x2 = "incub|A"))
  expect_equal(augm$value, NA_real_)

  # overwriting warnings
  base |>
    add_latent("dummy", func = onset - incub) |>
    add_latent("dummy", func = "log", func_args = c(x = infec@A)) |>
    expect_warning("[Oo]verwriting")
  base |>
    add_latent("dummy", func = onset - incub) |>
    add_latent("rate", value = 3) |>
    add_latent("dummy", func = exp(-rate)) |>
    expect_warning("[Oo]verwriting")
})

test_that("Latent load enable special definitions", {
  base <- mk_latent$model |>
    add_latent(infec = list(model = "normal", mean = 0, sd = 1))

  # alternative function call definitions
  augm <- base |>
    add_latent("dummy", func = "log", func_args = c(x = infec@A))
  tmp <- base |> add_latent(dummy = log(infec@A))
  expect_equal(get_latent(tmp), get_latent(augm))

  expect_warning(
    overwritten <- add_latent(tmp, name = "dummy", func = log(infec@A)),
    "[Oo]verwriting"
  )
  expect_equal(get_latent(overwritten), get_latent(augm))

  # multi-status definition: one row per status, identical otherwise
  multi <- base |>
    add_latent(dummy = list(func = log(infec@A), status = c(1, 3))) |>
    get_latent() |>
    subset(name == "dummy")
  expect_equal(multi$status, c("1", "3"))
  multi$status <- NULL
  expect_equal(multi[1, , drop = TRUE], multi[2, , drop = TRUE])
})

test_that("Latent support rlang double bang", {
  base <- mk_latent$model |>
    add_latent(infec = list(model = "normal", mean = 0, sd = 1))
  external_value <- 33
  augm <- base |>
    add_latent(
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
    get_latent()

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
  base <- mk_latent$model |>
    add_latent(indiv_shape = list(model = "uniform", min = 0, max = 2))
  val <- 1
  aux <- base |>
    add_latent(
      var1 = list(model = "normal", mean = indiv_shape + 1),
      var2 = list(model = "normal", mean = indiv_shape + !!val),
      var3 = list(model = "normal", mean = ~ indiv_shape + 1),
      var4 = list(model = "normal", mean = ~ indiv_shape + !!val)
    ) |>
    get_latent() |>
    subset(name != "indiv_shape")
  aux$name <- NULL
  # check that all rows are identical
  expect_equal(nrow(unique(aux)), 1)

  # an external list is loaded through 'import', not as a bare definition
  external <- list(
    some_func = list(model = "gamma", shape = ~shape, scale = ~scale)
  )
  shapes <- mk_latent$model |>
    add_latent(
      shape = list(model = "uniform", min = 1, max = 2),
      scale = list(model = "uniform", min = 1, max = 2)
    )
  shapes |>
    add_latent(external, append = FALSE) |>
    expect_error("non.?empty")

  augm <- shapes |>
    add_latent(import = external) |>
    get_latent() |>
    subset(name == "some_func")
  expect_equal(augm$model, "gamma")
  expect_equal(augm$args[[1]], c(shape = "shape|A", scale = "scale|A"))
})

test_that("Observation definitions are encoded", {
  obs <- mk_latent$model |>
    add_latent(onset_time = list(model = "normal", mean = 0, sd = 1)) |>
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
  base <- mk_latent$model
  # observation nodes are data-bound, never sampled
  base |>
    add_observation(x = list(model = "normal", proposal = "forward")) |>
    expect_warning("[Ii]gnored.*proposal")
})

test_that("Deterministic observation node is encoded as a function obs", {
  obs <- mk_latent$model |>
    add_latent(
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

test_that("Observation and latent loads each own only their definitions", {
  hhmodel <- mk_latent$model |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # non-append observation definition keeps latent definitions
    add_observation(
      onset_date = list(model = "normal", mean = incub_period@A, sd = 1),
      append = FALSE
    )
  expect_setequal(
    get_latent(hhmodel)[["name"]],
    c("incub_period", "infec_time")
  )
  expect_identical(get_observations(hhmodel)[["name"]], "onset_date")

  # symmetrically, a non-append latent load keeps observations
  hhmodel <- hhmodel |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      append = FALSE
    )
  expect_identical(get_latent(hhmodel)[["name"]], "incub_period")
  expect_identical(get_observations(hhmodel)[["name"]], "onset_date")
})

test_that("Household observations reject a distribution", {
  hhmodel <- mk_latent$model
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

test_that("Household observations reject a status", {
  hhmodel <- mk_latent$model
  err <- "cannot be status specific"

  # every input branch is guarded, at the bulk default level and per definition
  add_observation(hhmodel,
    name = "raiden", func = min(test_date@D), status = 1, level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    raiden = list(func = min(test_date@D), status = 1), level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    import = list(raiden = list(func = quote(min(test_date@D)), status = 1)),
    level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    raiden = list(func = min(test_date@D), status = 1, level = "house")
  ) |>
    expect_error(err)
  add_observation(hhmodel, raiden = list(func = min(test_date@D)),
                  status = 1, level = "house") |>
    expect_error(err)

  # a member-level declaration keeps its status, and the no-effect reduction
  # that would promote it to household level is not taken when one is set
  member <- hhmodel |>
    add_observation(raiden = list(func = min(test_date@D), status = 1)) |>
    get_observations()
  expect_identical(member[["level"]], "member")
  expect_identical(as.integer(member[["status"]]), 1L)

  # a bulk 'status =' reaches every entry, single and bulk modes alike
  bulk <- hhmodel |>
    add_observation(raiden = list(func = min(test_date@D)), status = 1) |>
    get_observations()
  expect_identical(as.integer(bulk[["status"]]), 1L)
  imported <- hhmodel |>
    add_observation(import = list(raiden = list(func = quote(min(test_date@D)))),
                    status = 1) |>
    get_observations()
  expect_identical(as.integer(imported[["status"]]), 1L)
})

test_that("A name is declared under a single kind", {
  hhmodel <- mk_latent$model
  err <- "already declared as"

  # latent and observations share one store: every pair collides, either way
  hhmodel |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(ocelot = list(func = test_date@D)) |>
    expect_error(err)
  hhmodel |>
    add_observation(ocelot = list(func = test_date@D)) |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1)) |>
    expect_error(err)
  hhmodel |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(ocelot = min(test_date@D), level = "house") |>
    expect_error(err)
  hhmodel |>
    add_observation(ocelot = min(test_date@D), level = "house") |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1)) |>
    expect_error(err)
  hhmodel |>
    add_observation(ocelot = list(func = test_date@D)) |>
    add_observation(ocelot = min(test_date@D), level = "house") |>
    expect_error(err)

  # a status-specific declaration is no loophole: the store is keyed by name
  hhmodel |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1, status = 1)) |>
    add_observation(ocelot = min(test_date@D), level = "house") |>
    expect_error(err)

  # the reserved status node is a declaration like any other
  hhmodel |>
    add_latent(.status = list(value = 1, const = TRUE)) |>
    add_status(a = 1) |>
    expect_error(err)

  # siblings of one kind stay legal: statuses of the same name, and overwrites
  hhmodel |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1, status = 1)) |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1, status = 2)) |>
    expect_silent()
  hhmodel |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 1)) |>
    add_latent(ocelot = list(model = "normal", mean = 0, sd = 2)) |>
    expect_warning("Overwriting")

  # an observation standing in for a like-named data column is the surrogate
  # design, not a clash: only declarations collide with each other
  hhmodel |>
    add_observation(test_date = list(func = onset_date@D)) |>
    expect_silent()
})

test_that("Level is set per definition, the argument being the bulk default", {
  # a per-definition 'level' picks the kind of its own entry, exactly as
  # 'status' does; regression: it used to be swallowed by add_quo's dots and
  # folded into the definition arguments, leaving every entry on the default
  mixed <- mk_latent$model |>
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
  bulk <- mk_latent$model |>
    add_observation(
      start_date = list(func = min(test_date@D)),
      end_date = list(func = 21, const = TRUE, level = "member"),
      level = "house"
    ) |>
    get_observations()
  expect_identical(bulk[["level"]][bulk[["name"]] == "start_date"], "house")
  expect_identical(bulk[["level"]][bulk[["name"]] == "end_date"], "member")

  # imported definitions honour it too
  imported <- mk_latent$model |>
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
  mk_latent$model |>
    add_observation(start_date = list(func = 1, level = "liquid")) |>
    expect_error("should be one of")

  # 'level' is meaningless outside observations
  mk_latent$model |>
    add_latent(ocelot = list(value = 1, level = "house")) |>
    expect_error("reserved to observations")
})

test_that("Household observations do not interfer with latent defs", {
  hhmodel <- mk_latent$model |>
    add_latent(incub_period = list(model = "gamma", shape = 2, scale = 1)) |>
    add_observation(
      onset_date = list(model = "normal", mean = incub_period@A, sd = 1)
    )

  # non-append household observation keeps latent definitions
  hhmodel <- add_observation(
    hhmodel, first_test = min(test_date@D), level = "house", append = FALSE
  )
  expect_identical(get_latent(hhmodel)[["name"]], "incub_period")
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
  expect_identical(get_latent(hhmodel)[["name"]], "incub_period")
})

test_that("Observation-aware initialization", {
  window <- 0.1
  # latent z ~ U(0, 1)
  # observation x ~ U(z - w, z + w)
  # if x = 0.5, valid for z in [0.4, 0.6]
  # without observation-aware restart, it would fail with 80% proba per-chain
  # probably of failing on K chains: 1 - 0.2^K
  data.frame(house_id = "shadow_moses", status = 1, x = 0.5) |>
    new_housepi_model("test") |>
    add_latent("z", model = "uniform") |>
    add_observation("x", model = "uniform", loc = z,
                    min = -!!window, max = !!window) |>
    run_inference(steps = 10L, adapt = 0L, chains = 40L) |>
    suppressWarnings() |>
    expect_no_error()

  # datum x = 5 is unreachable
  # UPDATE: it's the same as above: new-dataset decorator?
  data.frame(house_id = "outer_heaven", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_latent("z", model = "uniform") |>
    add_observation("x", model = "uniform", loc = z,
                    min = -!!window, max = !!window) |>
    run_inference() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", prior = "uniform") |>
    add_observation("x", model = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    expect_error("[Ff]ree parameter dependencies")

  # generic error
  data.frame(house_id = "metal_gear", status = 1, x = 5) |>
    new_housepi_model("test") |>
    add_parameter("theta", value = 0.3) |>
    add_latent("z", model = "uniform") |>
    add_observation("x", model = "uniform", min = -theta@P, max = theta@P) |>
    run_inference() |>
    expect_error("[Oo]utside of respective domain")
})

test_that("Reading an observation in models and loops are rejected", {

  # a latent model must not read an observation (illegal in both
  # modes; caught at eager build, so the construction itself errors)
  data.frame(house_id = "outer_heaven", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", model = "uniform", min = 0, max = 1) |>
    add_latent("latent", model = "normal", mean = onset@D) |>
    expect_error("[Mm]odel.*reads.*observation")

  # a deterministic read that closes a loop (the observation's support
  # depends on the reader) is circular
  data.frame(house_id = "metal_gear", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_latent("ocelot", value = onset@D) |>
    add_observation("onset", model = "uniform", min = ocelot@A, max = 10) |>
    expect_error("[Cc]yclic")
})

test_that("simulation override (value + model) rejects forward dependencies", {
  # reads an observation
  data.frame(house_id = "shadow_moses", status = 1, onset = 0.5) |>
    new_housepi_model("test") |>
    add_observation("onset", model = "uniform", min = 0, max = 1) |>
    # forward dependency: observation before latent
    add_latent("liquid", model = "normal", value = onset@D) |>
    expect_error("[Ss]imulation override.*observation")

  # reads another latent (no observation present: guard still fires)
  data.frame(house_id = "outer_heaven", status = 1) |>
    new_housepi_model("test") |>
    add_latent(liquid = list(model = "normal", value = solid@A),
               solid = list(model = "gamma", shape = 2, scale = 1)) |>
    expect_error("[Ss]imulation override.*latent")
})

test_that("degenerate reduction at member-level is cast to house-level", {
  # member-level reduction by mistake
  expect_warning(
    promoted <- mk_latent$model |>
      add_observation(first_test = min(test_date@D)) |>
      get_observations(),
    regexp = "no-effect.*reduction.*min"
  )
  expect_identical(promoted[["level"]], "house")
  # explicit house-level reduction
  mk_latent$model |>
    add_observation(first_test = min(test_date@D), level = "house") |>
    get_observations() |>
    expect_identical(promoted)

  # redundant reduction raise a error and are not cast
  mk_latent$model |>
    add_observation(first_test = min(test_date@D), level = "house") |>
    add_observation(codec = min(first_test@D)) |>
    expect_error("redundant.*codec")

  # genuine per-member expressions pass normally
  kept <- mk_latent$model |>
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

test_that("An undeclared producer defers the build to the wall", {
  # parameters land in the engine on their own decorator, so a definition
  # reading one that is not declared yet waits instead of failing a key lookup
  raiden <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_latent(z = list(model = "gamma", shape = shape@P, scale = scale@P))
  expect_silent(raiden)

  # the wall names the missing parameters, not the key that could not be found
  sms <- tryCatch(evaluate_log_lik(raiden), error = conditionMessage)
  expect_match(sms, "missing parameter definitions")
  expect_false(grepl("status-specific", sms))
})

test_that("A household observation with no datum yet reads as NA", {
  # 'start_date' reduces over an observation the run has still to draw: there
  # is nothing to reduce, so the data view reports NA rather than reading the
  # empty reduction as a datum failing the discreteness 'case_ascertained' asks
  hhmodel <- data.frame(
    id = c("snake", "meryl", "campbell", "naomi"),
    house_id = c("shadow", "shadow", "moses", "moses"),
    status = NA, end_date = 21
  ) |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_observation(
      onset_date = list(model = "normal", discr = TRUE, mean = 3, sd = 1),
      start_date = list(value = min(onset_date@D), level = "house")
    )

  df <- expect_silent(get_data_wide(hhmodel))
  expect_true(all(is.na(df$start_date)))
})

test_that("A surrogate column proves the requirement its model asks", {
  # a surrogate filling 'start_date' answers the discreteness 'case_ascertained'
  # asks of the column, and answers for what it can return, not for the value it
  # happens to settle on
  # a household-level definition reduces over its members, so its root is the
  # reduction; the member level takes the algebra directly
  mk <- function(defin, level = "house") {
    data.frame(
      id = c("snake", "meryl", "campbell", "naomi"),
      house_id = c("shadow", "shadow", "moses", "moses"),
      status = c(1, 0, 1, 0), end_date = 16, half = 0.5
    ) |>
      new_housepi_model("case_ascertained", quiet = TRUE) |>
      add_observation("start_date", value = !!defin, level = level) |>
      get_data_wide()
  }

  # integer-valued cells fold to the integer lattice through the algebra
  expect_silent(mk(quote(min(end_date@D))))
  expect_silent(mk(quote(sum(end_date@D))))
  expect_silent(mk(quote(end_date@D * 2 - 1), level = "member"))

  # two halves sum to an integer, and are refused all the same: the column
  # answers for what the definition can return, not for the value it settles on
  mk(quote(sum(half@D))) |>
    expect_error("not discrete-valued")
  mk(quote(end_date@D / 4), level = "member") |>
    expect_error("not discrete-valued")
  mk(quote(sqrt(end_date@D)), level = "member") |>
    expect_error("not discrete-valued")

  # rounding restores the proof
  expect_silent(mk(quote(round(end_date@D / 4)), level = "member"))
  expect_silent(mk(quote(ceil(sqrt(end_date@D))), level = "member"))

  # the refusal names the definition and the place it was built for
  sms <- tryCatch(mk(quote(sum(half@D))), error = conditionMessage)
  expect_match(sms, "[Hh]ousehold observation")
  expect_match(sms, "shadow")
  sms <- tryCatch(mk(quote(end_date@D / 4), level = "member"),
                  error = conditionMessage)
  expect_match(sms, "individual")
  expect_match(sms, "snake")
})

test_that("A discrete model resolves 'uniform' to the discrete distribution", {
  outer_heaven <- data.frame(
    id = c("snake", "meryl", "otacon"), house_id = "shadow_moses",
    status = c(1, 1, 0), start_date = 0, end_date = 10
  )
  build <- function(min) {
    outer_heaven |>
      new_housepi_model("case_ascertained", discr = TRUE, quiet = TRUE) |>
      add_parameter(import = "template") |>
      add_latent(infec_time = list(model = "uniform", min = !!min, max = 3)) |>
      add_infectivity_profile("uniform", max = 4) |>
      evaluate_log_lik()
  }

  expect_silent(build(0))

  # the flat lattice law owns every integer value in the window, so a bound that
  # is not one has no reading, and the error names the latent that carries it
  sms <- tryCatch(build(0.5), error = conditionMessage)
  expect_match(sms, "min.*uniform.*discrete")
  expect_match(sms, "infec_time")
})

test_that("A location off the integers is refused when discrete", {
  hhmodel <- data.frame(house_id = 1, id = c("a", "b"), status = c(1, 0)) |>
    new_housepi_model("test")

  # a sampled location would shift the lattice under the chain at every step:
  # the fix is on the declaration it reads, named in place of `loc`
  hhmodel |>
    add_parameter("ocelot", prior = "uniform", min = 0, max = 1) |>
    add_latent("quiet", model = "poisson", lambda = 2, loc = ocelot@P) |>
    evaluate_log_lik() |>
    expect_error("[Cc]hange parameter.*ocelot.*discrete")
  # a fractional constant is refused just the same, but carries no declaration
  # to flag: only its value can move
  hhmodel |>
    add_latent(name = "quiet", model = "poisson", lambda = 2, loc = 0.5) |>
    evaluate_log_lik() |>
    expect_error("[Ss]et.*loc.*integer")
})

test_that("Latent aliasing works as expected", {
  # check that model arguments are correctly integrated
  simulated <- data.frame(id = seq(1, 1000),
                          house_id = rep(seq(10), each = 100),
                          status = 1) |>
    new_housepi_model("test") |>
    add_latent(hour = list(model = "uniform", min = day, max = day + 1),
               day = list(model = "gamma", shape = 2, scale = 1)) |>
    run_simulation()
  expect_true(all(
    simulated$day <= simulated$hour & simulated$hour <= simulated$day + 1
  ))
})

test_that("value paired with model overrides the draw in simulations only", {
  hhmodel <- data.frame(id = seq_len(100), house_id = 1,
                        status = 1, anchor = 7) |>
    new_housepi_model("test") |>
    add_latent("eli", model = "uniform", max = 100, value = anchor@D)

  # simulation: pins the sampled variable to user value
  expect_true(all(run_simulation(hhmodel)$eli == 7))

  # inference: the value does not pin the sampler and draws span the prior
  draws <- hhmodel |>
    run_inference(steps = 40L, chains = 1L) |>
    get_latent_draws() |>
    suppressWarnings()
  expect_false(all(draws[["1_eli"]] == 7))
})

test_that("load latent list works fine for external lists", {
  hhmodel <- data.frame(id = "id1", house_id = "h1",
                        status = 1, date = as.Date("2021-01-01")) |>
    new_housepi_model("test")

  # an external list of definitions can be imported
  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  tmp <- add_latent(hhmodel, import = external_list)$latent

  hhmodel <- add_latent(hhmodel, z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$latent)
})

test_that("latent template import loads", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_latent(pupa = start_date@D + 1)

  # appending keeps the pre-existing definitions beside the template ones
  augm <- hhmodel |>
    add_latent(import = "template") |>
    get_latent()
  augm$name |>
    expect_setequal(c("pupa", "infec_time"))

  # dropping resets the store first
  augm <- hhmodel |>
    add_latent(import = "template", append = FALSE) |>
    get_latent()
  augm$name |>
    expect_setequal("infec_time")
})

test_that("Latent traces export only relevant variables", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(alpha = 1, beta = list(prior = "uniform")) |>
    add_latent(
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
      # equal to other latent objects
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
    steps = 200L, burn = 0L, thin = 1L, adapt = 0L, chains = 1L,
    id_draws = infec_df$id
  ))
  expect_equal(
    posterior::variables(result$latent_draws),
    paste(infec_df$id, c("layer", "alias", "alpha_beta", "expl"), sep = "_")
  )
})

test_that("A latent degenerate domain warns and collapses", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter("d_dog", prior = "uniform") |>
    add_latent("d_horse", model = "uniform", min = 5, max = 5)
  expect_warning(
    tryCatch(
      run_inference(hhmodel, steps = 10L, chains = 1L, adapt = 0L, thin = 1L),
      error = function(e) NULL
    ),
    regexp = "[Dd]egenerate domain|[Cc]ollapsed"
  )
})

test_that("Nested definitions works as expected", {
  hhmodel <- new_housepi_model(model = "test") |>
    # define Box-Muller transformation
    add_latent(norm = list(model = "normal"),
               u1 = list(model = "uniform"),
               u2 = list(model = "uniform"),
               z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2))
  expect_snapshot(print(hhmodel))

  results <- hhmodel |>
    run_inference(steps = 1e5L, chains = 1L, adapt = 0L, thin = 1L)
  # no parameters nor priors, only latent
  expect_equal(posterior::nvariables(results$parameter_draws), 0L)
  posterior::nvariables(results$latent_draws) |>
    expect_equal(dim(hhmodel$latent)[1])

  # test S3 version of summarise_draws
  out_vars <- posterior::summarise_draws(results)[["variable"]]
  expect_true(all(stringr::str_starts(out_vars, "tommy")))

  # empirical tests for goodness of Gaussian distribution
  sample_norm <- results$latent_draws$tommy_norm
  sample_boxmul <- results$latent_draws$tommy_z0

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

test_that("Nested latent definitions raise errors", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("test")

  hhmodel |>
    add_latent(.layer = list(model = "uniform", status = 1),
               .layer = list(value = .dummy, status = c(2, 3, 4)),
               .dummy = list(model = "uniform"),
               z = .layer) |>
    evaluate_log_lik() |>
    expect_no_error()

  # simple nest
  hhmodel |>
    add_latent(.dummy = z, z = .dummy) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic latent definition")

  # triangular nest
  hhmodel |>
    add_latent(z = .layer, .layer = .dummy, .dummy = z) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic latent definition")

  # first generation of definitions
  hhmodel |>
    add_latent(.layer = list(model = "normal", status = 1),
               .layer = list(value = .dummy, status = 2),
               #> NO .layer provided for status 3 and 4
               .dummy = list(model = "lognormal", status = 2),
               z = .layer) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")

  # second generation of definitions
  hhmodel |>
    add_latent(.layer = list(model = "normal", status = 1),
               .layer = list(value = .dummy, status = c(2, 3, 4)),
               #> NO .dummy provided for status 3, 4
               .dummy = list(model = "lognormal", status = c(2, 3)),
               z = .layer,
               append = FALSE) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")
})

test_that("Missing latent function arguments raise errors", {
  hhmodel <- data.frame(house_id = seq(3), id = seq(3), status = 1) |>
    new_housepi_model("test")

  # non default distribution arguments must be supplied
  hhmodel |>
    add_latent(name = "z", func = "gamma") |>
    expect_error("[Mm]issing non.?default")
})

test_that("Status-specific latent variables draw from per-status ranges", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("test") |>
    add_latent(z = list(status = 1, model = "uniform", min = 0, max = 1),
               z = list(status = 2, model = "uniform", min = 1, max = 2),
               z = list(status = 3, model = "uniform", min = 2, max = 3),
               z = list(status = 4, model = "uniform", min = 3, max = 4))

  # simulation routes draws into per-status supports
  withr::with_seed(1, sim <- run_simulation(hhmodel, nrep = 50))
  sim <- do.call(rbind, sim)
  for (s in seq(4)) {
    vals <- sim$z[sim$status == s]
    expect_true(all(vals >= s - 1 & vals <= s))
  }

  # short-run diagnostic fires (steps = 5 < 400)
  expect_warning(
    results <- hhmodel |>
      run_inference(steps = 5L, thin = 1L, adapt = 2L,
                    chain = 1L, id_draws = hhmodel$data_wide$id),
    regexp = "below 400"
  )
  draws <- as.data.frame(get_latent_draws(results))
  for (s in seq(4)) {
    col <- draws[[as.character(s)]]
    expect_true(all(col >= s - 1 & col <= s))
  }
})

test_that("A latent bound to a distribution returns its density, not a draw", {
  # a distribution called as a func gives back a mass: reading its support here
  # would type the column integer and truncate every mass to zero
  hhmodel <- data.frame(house_id = rep(seq(3), each = 2), id = seq(6),
                        status = 1) |>
    new_housepi_model("test") |>
    add_latent(k = list(model = "uniform", min = 0, max = 5, discr = TRUE),
               w = list(func = "poisson", func_args = c(x = k, lambda = 3)))

  sim <- run_simulation(hhmodel, nrep = 1)
  expect_type(sim$k, "integer")
  expect_type(sim$w, "double")
  expect_equal(sim$w, dpois(sim$k, 3))
})

test_that("export_augmented_data works after evaluate_log_lik", {
  hhmodel <- data.frame(house_id = 1, id = c("a", "b"), status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent(name = "z", model = "uniform")
  evaluate_log_lik(hhmodel)
  result <- hhmodel$export_augmented_data()
  expect_true(is.data.frame(result))
  expect_true("z" %in% colnames(result))
})

test_that("Incubation augmentation recast runs on the test model", {
  # a minimal analogue of the data-augmentation vignette: the incubation gamma
  # enters as a latent prior, the reported onset day as an observation node.
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
    add_latent(
      .infec_time = list(model = "uniform", max = exposure_max@D + 0.99),
      incub_period = list(
        model = "gamma", mean = mean@P, sd = sd@P, min = 0.5, max = 14
      ),
      .onset_time = .infec_time + incub_period
    ) |>
    add_observation("onset_mid", model = "normal",
                    mean = .onset_time, sd = 0.3) |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  # the gamma incubation prior contributes to the (finite) latent log-prior
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
    add_latent(z = list(value = w@D * 2))

  # the data view reports the derived value before evaluation
  expect_identical(unique(get_data_wide(hhmodel)[["start_date"]]), 2)
  expect_false("w" %in% colnames(get_data_wide(hhmodel)))

  # and after evaluation
  expect_no_error(evaluate_log_lik(hhmodel))
  expect_identical(unique(get_data_wide(hhmodel)[["start_date"]]), 2)
  expect_false("w" %in% colnames(get_data_wide(hhmodel)))

  # simulated observations are for all members, latent just for infected
  simulated <- run_simulation(hhmodel)
  expect_equal(simulated[["w"]], c(3, 7, NA))
  expect_equal(simulated[["z"]], c(6, NA, NA))

  # FIXME: the requirement reported state "latent", but it's "observation"
  # a stochastic observation carries no column to read
  #  datum is reported missing, rather than failing on the key it stands for
  data.frame(house_id = "rex", status = 1, otacon = 5) |>
    new_housepi_model("test") |>
    add_latent(y = list(model = "uniform")) |>
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
    add_latent(z = list(model = "lognormal"),
               y = list(model = "uniform", max = first@D)) |>
    add_observation(t_obs = list(model = "normal", mean = z, sd = 1)) |>
    add_observation("first", value = min(t_obs@D), level = "house") |>
    run_simulation() |>
    expect_error("cyclic dependency")
})

test_that("a surrogate column reading latent is simulation only", {
  hhmodel <- data.frame(id = seq_len(100), house_id = "deck",
                        status = 1, ReD = 0) |>
    new_housepi_model("test") |>
    add_latent(z = list(model = "uniform")) |>
    add_observation(ReD = z + 2)

  # the data view has no derived value to report: it keeps the imported datum
  expect_true(all(get_data_wide(hhmodel)[["ReD"]] == 0))

  # observation depends on sampled latent: not a valid surrogate definition
  evaluate_log_lik(hhmodel) |>
    expect_error("must resolve from data alone")

  # simulation produces it, following the latent variable as it is drawn
  simulated <- run_simulation(hhmodel)
  expect_equal(simulated[["ReD"]], simulated[["z"]] + 2)
})

test_that("Unspecified discretization reads back resolved to continuous", {
  obs <- mk_latent$model |>
    add_observation(test_date = list(model = "normal"),
                    onset_date = list(model = "normal", discr = TRUE),
                    expo = list(model = "normal", discr = FALSE)) |>
    get_observations()

  expect_type(obs$discr, "logical")
  expect_equal(obs$discr, c(FALSE, TRUE, FALSE))
})

test_that("An observation infers discretization from the column it reads", {
  hhmodel <- data.frame(house_id = rep(c("shadow", "moses"), each = 20),
                        date_col = as.Date("2005-02-28") + rep(0:1, 20),
                        int_col = rep(c(3L, 5L), 20),
                        whole_col = rep(c(3, 5), 20),
                        kept_col = rep(c(3L, 5L), 20),
                        status = 1) |>
    new_housepi_model("test") |>
    add_latent(z = list(model = "uniform", min = 2, max = 8)) |>
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
  expect_true(on_integer_lattice(simulated[["date_col"]]))
  expect_true(on_integer_lattice(simulated[["int_col"]]))
  expect_false(on_integer_lattice(simulated[["whole_col"]]))
  expect_false(on_integer_lattice(simulated[["kept_col"]]))

  # a deterministic observation is not influenced
  hhmodel |>
    add_observation(derived = z + 1) |>
    run_simulation() |>
    getElement("derived") |>
    on_integer_lattice() |>
    expect_false()
})

test_that("Simulated columns carry the discretization of their definition", {
  hhmodel <- data.frame(house_id = rep(c("otacon", "meryl"), each = 3),
                        id = seq(6), status = 1) |>
    new_housepi_model("test") |>
    add_latent(z = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      codec = list(model = "normal", discr = TRUE, mean = z, sd = 1),
      nomad = list(model = "normal", mean = z, sd = 1),
      first_codec = list(value = min(codec@D), level = "house")
    )

  withr::with_seed(5, simulated <- run_simulation(hhmodel))

  # a discrete declaration, and the reduction folding it, export as integers
  expect_type(simulated[["codec"]], "integer")
  expect_type(simulated[["first_codec"]], "integer")
  # a continuous observation and a continuous latent stay double
  expect_type(simulated[["nomad"]], "double")
  expect_type(simulated[["z"]], "double")

  # feeding the simulated data back resolves the same discretization
  simulated |>
    new_housepi_model("test") |>
    add_latent(z = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(codec = list(model = "normal", mean = z, sd = 1)) |>
    get_observations() |>
    getElement("discr") |>
    expect_true()
})

test_that("A column a status declares continuous keeps its fractions", {
  hhmodel <- data.frame(house_id = rep(c("otacon", "meryl"), each = 2),
                        id = seq(4), status = rep(c(1, 2), 2)) |>
    new_housepi_model("test") |>
    add_latent(z = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      codec = list(status = 1, model = "normal", discr = TRUE, mean = z),
      codec = list(status = 2, model = "normal", mean = z)
    )

  withr::with_seed(3, simulated <- run_simulation(hhmodel))

  # one column, two declarations: the continuous one has values to lose
  expect_type(simulated[["codec"]], "double")
  expect_false(on_integer_lattice(simulated[["codec"]][simulated$status == 2]))
})

test_that("Discretization is resolved every time the data changes", {
  df <- data.frame(house_id = rep(c("naked", "snake"), 20), status = 1)

  hhmodel <- new_housepi_model(model = "test") |>
    add_latent(z = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(tsuchinoko = list(model = "normal", mean = z, sd = 2))

  # observations resolve to continuos in simulations and definitions
  hhmodel <- hhmodel |>
    add_data_wide(dplyr::mutate(df, tsuchinoko = rep(c(2.4, 5.1), 20)))
  get_observations(hhmodel)$discr |>
    expect_false()
  on_integer_lattice(run_simulation(hhmodel)$tsuchinoko) |>
    expect_false()

  # observations resolve to discrete in simulations and definitions
  hhmodel <- hhmodel |>
    add_data_wide(dplyr::mutate(df, tsuchinoko = rep(c(2L, 5L), 20)))
  get_observations(hhmodel)$discr |>
    expect_true()
  on_integer_lattice(run_simulation(hhmodel)$tsuchinoko) |>
    expect_true()
})
