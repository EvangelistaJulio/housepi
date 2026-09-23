library(testthat)

# model carrying the data columns the definition blocks refer to
mk_latent <- local({
  hhmodel <- data.frame(hhid = "h1", status = 1, datum = 1,
                        record = 3, datum_2 = 2) |>
    new_housepi_model("test")
  list(model = hhmodel)
})

# saturated outbreak: every household infects all of its members
mk_outbreak <- local({
  hhdata <- sarscov2_2021_icofs |>
    dplyr::filter(hhid %in% utils::head(unique(hhid), 50)) |>
    dplyr::mutate(start_date = 0, end_date = 21)
  hhmodel <- hhdata |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-12, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    )
  list(data = hhdata, model = hhmodel)
})

# case-ascertained household sampling incubation periods and infection times
mk_onset <- local({
  hhdata <- subset(helper_onset_house$data, select = -infec_time)
  hhmodel <- hhdata |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1, max = 14),
      infec_time = list(model = "uniform", min = -14, max = end_date@D)
    ) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)
  list(data = hhdata, model = hhmodel)
})

# two infected statuses under the template definitions
mk_status <- local({
  hhdata <- data.frame(id = c("m1", "m2", "m3", "m4"), hhid = "h1",
                       status = c(1, 2, 0, 0), start_date = 0, end_date = 21)
  hhmodel <- hhdata |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_latent(import = "template")
  list(data = hhdata, model = hhmodel)
})

# every case_ascertained parameter pinned
mk_pinned <- local({
  hhmodel <- new_housepi_model(model = "case-ascertained") |>
    add_parameter(import = helper_onset_house$param) |>
    add_latent(infec_time = infec_time@D) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)
  list(model = hhmodel)
})

on_integer_lattice <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  length(x) > 0 && all(x == floor(x))
}

# draws of a run as a data frame of its latent traces
latent_frame <- function(results) {
  results |>
    get_latent_draws() |>
    posterior::as_draws_matrix() |>
    as.data.frame()
}

##############################
###  LATENT : DEFINITIONS  ###
##############################

test_that("Latent definition errors are recognized", {
  # a decorator carrying no definition at all
  mk_latent$model |>
    add_latent() |>
    expect_error("[Mm]issing or empty latent definitions")
  mk_latent$model |>
    add_latent(import = list()) |>
    expect_error("[Mm]issing or empty latent definitions")
  # not valid definitions
  mk_latent$model |>
    add_latent("hidden@") |>
    expect_error("[Ii]nvalid.*characters")
  mk_latent$model |>
    add_latent("hidden") |>
    expect_error("[Ii]ncomplete definition")
  mk_latent$model |>
    add_latent("") |>
    expect_error("non.?empty")
  # non default distribution arguments must be supplied
  mk_latent$model |>
    add_latent("hidden", func = "gamma") |>
    expect_error("[Mm]issing non.?default")
  # options typed as values are protected
  regexp <- "[Aa]rgument.*must be.*[Gg]ot.*instead"
  mk_latent$model |>
    add_latent("hidden", model = "uniform", discr = NULL) |>
    expect_error(regexp)
  mk_latent$model |>
    add_latent("hidden", model = "uniform", discr = sum) |>
    expect_error(regexp)
  mk_latent$model |>
    add_latent("hidden", model = "normal", status = sum) |>
    expect_error(regexp)
  mk_latent$model |>
    add_latent(hidden = list(model = "uniform", discr = NULL)) |>
    expect_error(regexp)

  mk_latent$model |>
    add_latent(import = list(list(model = "normal"))) |>
    expect_error("non-empty")
})

test_that("Undefined definitions are controlled", {
  mk_latent$model |>
    add_latent("hidden", proposal = "forward", value = NA) |>
    expect_error("[Nn]o .model. to sample from")
  mk_latent$model |>
    add_latent("hidden", value = NA) |>
    expect_error("[Mm]issing value")
  # a plain value declared variable has no model to draw from
  mk_latent$model |>
    add_latent("hidden", value = 3, const = FALSE) |>
    expect_error("[Nn]o .model. to sample from")

  defs <- mk_latent$model |>
    add_latent("hidden", value = 3, const = TRUE, proposal = "forward") |>
    get_latent()
  expect_true(is.na(defs$proposal))
})

test_that("Passing a prior raises a hint on model", {
  sms <- mk_latent$model |>
    add_latent(hidden = list(prior = "uniform", min = -1, max = 2)) |>
    tryCatch(error = conditionMessage)
  expect_match(sms, "[Ss]pare.*prior.*hidden")
  expect_match(sms, "[Dd]id you mean.*model")
})

test_that("A func definition is deterministic", {
  # a plain number is a constant value, whichever argument carried it
  by_func <- mk_latent$model |> add_latent(hidden = list(func = 21))
  by_value <- mk_latent$model |> add_latent(hidden = list(value = 21))
  expect_equal(get_latent(by_func), get_latent(by_value))

  # deterministic definitions reject model or proposal
  mk_latent$model |>
    add_latent("hidden", func = "log", x = 2,
               model = "normal", mean = 0, sd = 1) |>
    expect_error("only one.*func.*model")
  mk_latent$model |>
    add_latent("hidden", func = "log", x = 2, proposal = "rw")   |>
    expect_error("[Dd]eterministic.*hidden.*cannot be sampled")
  # nor a value alongside the function
  mk_latent$model |>
    add_latent("hidden", value = 2, func = "log", x = 1) |>
    expect_error("[Oo]nly one between .value. and .func.")
})

test_that("A name carrying no source sigil is resolved, not guessed", {
  # a bare name reaches the data column declaring it
  obs <- mk_latent$model |>
    add_observation(record_2 = min(datum_2)) |>
    get_observations()
  expect_equal(obs[["func_args"]][[1]], c(x = "datum_2|D"))

  # and a bare name reaches the parameter declaring it
  obs <- mk_latent$model |>
    add_parameter(rho = list(prior = "gamma", shape = 2, scale = 1)) |>
    add_observation(record = list(model = "normal", mean = rho, sd = 1)) |>
    get_observations()
  expect_equal(obs[["args"]][[1]], c(mean = "rho|P", sd = "1"))

  # a name two holders answer to is settled by the sigil, never by the engine
  both <- mk_latent$model |>
    add_parameter(datum = list(prior = "gamma", shape = 2, scale = 1))
  sms <- both |>
    add_observation(record = list(model = "normal", mean = datum, sd = 1)) |>
    evaluate_log_lik() |>
    tryCatch(error = conditionMessage)
  expect_match(sms, "[Cc]ould not resolve the source of .datum.")
  expect_match(sms, "[Tt]ag it .datum@D. or .datum@P.")
  both |>
    add_observation(record = list(model = "normal", mean = datum@P, sd = 1)) |>
    evaluate_log_lik() |>
    expect_no_error()

  # nothing of that name anywhere: the defect is named, with the near miss
  sms <- mk_latent$model |>
    add_observation(record = list(model = "normal", mean = datums, sd = 1)) |>
    evaluate_log_lik() |>
    tryCatch(error = conditionMessage)
  expect_match(sms, "[Cc]ould not resolve the source of .datums.")
  expect_match(sms, "[Dd]id you mean .datum.")

  # an explicitly mis-tagged sigil still reaches the engine hint
  mk_latent$model |>
    add_observation(record = list(model = "normal", mean = hidden@A)) |>
    evaluate_log_lik() |>
    expect_error("add_latent")
  # a latent read naming a column or a parameter hints at the sigil it misses
  mk_latent$model |>
    add_latent("hidden", model = "normal", mean = datum@A) |>
    evaluate_log_lik() |>
    expect_error("[Dd]id you mean .datum@D.")
  mk_latent$model |>
    add_parameter(rho = 1) |>
    add_latent(hidden = list(model = "normal", mean = rho@A)) |>
    evaluate_log_lik() |>
    expect_error("[Dd]id you mean .rho@P.")
})

test_that("A forward consumer waits resolves after declaration", {
  obs <- mk_latent$model |>
    add_observation(record = list(model = "normal", mean = rho, sd = 1)) |>
    get_observations()
  expect_equal(obs[["args"]][[1]], c(mean = "rho|?", sd = "1"))

  # the consumer is encoded before the parameter exists, and resolves after
  obs <- mk_latent$model |>
    add_observation(record = list(model = "normal", mean = rho, sd = 1)) |>
    add_parameter(rho = list(prior = "gamma", shape = 2, scale = 1)) |>
    get_observations()
  expect_equal(obs[["args"]][[1]], c(mean = "rho|P", sd = "1"))
})

test_that("A status left without definitions names the status", {
  # definitions are left unresolved until a wall
  mk_latent$model |>
    add_latent(hidden_2 = list(func = min(datum_2@D), status = 2),
               hidden_2 = list(func = min(datum_2@D), status = 1)) |>
    get_latent() |>
    nrow() |>
    expect_equal(2L)

  # every per-member definition is pinned to a status no individual carries
  mk_latent$model |>
    add_latent("hidden_2", func = min(datum_2@D), status = 2) |>
    evaluate_log_lik() |>
    expect_error("No.*latent.*definition.*status 1")

  # no status block to fall back on: the generic definition is the missing one
  data.frame(hhid = "h1", status = 0, datum_2 = 2) |>
    new_housepi_model("test") |>
    add_latent("hidden_2", func = min(datum_2@D), status = c(2, 3)) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*generic.*latent.*definition")
})

test_that("Latent variable definitions are correct", {
  # basic constant definition
  augm <- mk_latent$model |>
    add_latent("hidden", value = 3) |>
    get_latent()
  expect_equal(augm$name, "hidden")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$const, TRUE)
  expect_equal(augm$discr, FALSE)
  expect_equal(augm$value, 3)
  expect_equal(augm$proposal, NA_real_)
  # a bare value is the same constant, whatever its source
  expect_equal(get_latent(add_latent(mk_latent$model, hidden = 3)), augm)
  mk_latent$model |>
    add_latent(hidden = 3@D) |>
    get_latent() |>
    expect_equal(get_latent(add_latent(mk_latent$model, "hidden", value = 3@D)))

  # a model with no value is sampled by random walk
  augm <- mk_latent$model |>
    add_latent("hidden", model = "normal") |>
    get_latent()
  expect_equal(augm$const, FALSE)
  expect_equal(augm$proposal, "rw")
  # a misspelled family is caught at definition
  mk_latent$model |>
    add_latent("hidden", model = "norm") |>
    expect_error("[Uu]nknown distribution")

  # variable definition: each source resolves against its own family, so the
  # model only accepts it once 'rho' is a parameter and 'layer_1' a latent
  # variable ('datum' being a data column)
  augm <- mk_latent$model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 9) |>
    add_latent("layer_1", model = "uniform", min = 0, max = 1) |>
    add_latent("hidden", value = datum@D, model = "uniform",
               min = layer_1@A, max = rho@P) |>
    get_latent("hidden")
  expect_equal(augm$status, NA_real_)
  expect_equal(augm$value, "datum|D")
  expect_equal(augm$model, "uniform")
  expect_equal(augm$min, "layer_1|A")
  expect_equal(augm$max, "rho|P")
})

test_that("Latent independence proposal folds its arguments", {
  augm <- mk_latent$model |>
    add_latent("hidden", model = "uniform", min = 0, max = 5,
               proposal = "normal", proposal_mean = 0, proposal_sd = 1) |>
    get_latent()
  expect_equal(augm$proposal, "normal")
  expect_equal(augm$proposal_args[[1]], c(mean = 0, sd = 1))

  # random-walk scale is not accepted by an independence sampler
  mk_latent$model |>
    add_latent("hidden", model = "uniform", min = 0, max = 5,
               proposal = "normal", proposal_scale = 0.5) |>
    expect_error("random.?walk")
})

test_that("Independence proposal rejects an unsupported prior domain", {
  mk_latent$model |>
    add_parameter("rho", prior = "uniform") |>
    add_latent("hidden", model = "uniform", min = rho@P, max = 5,
               proposal = "gamma", proposal_shape = 2) |>
    evaluate_log_lik() |>
    expect_error("[Cc]annot cover")
})

test_that("Latent function definition is correct", {
  hhmodel <- mk_latent$model |>
    add_latent(layer_1 = list(model = "normal", mean = 0, sd = 1),
               layer_2 = list(model = "gamma", shape = 2, scale = 1))

  # difference definition
  augm <- hhmodel |>
    add_latent("hidden", func = layer_1 - layer_2) |>
    get_latent("hidden")
  expect_equal(augm$func, "sub")
  expect_equal(augm$func_args[[1]], c(x = "layer_1|A", x2 = "layer_2|A"))
  expect_equal(augm$value, NA_real_)

  # alternative function call definitions
  by_args <- hhmodel |>
    add_latent("hidden", func = "log", x = layer_1@A)
  by_expr <- hhmodel |> add_latent(hidden = log(layer_1@A))
  expect_equal(get_latent(by_expr), get_latent(by_args))

  # an expression carries its own arguments: spare ones have no reading
  hhmodel |>
    add_latent("hidden", func = log(layer_1@A), x = layer_1@A) |>
    expect_error("[Ee]xpression already carries its arguments")

  expect_warning(
    overwritten <- by_expr |> add_latent("hidden", func = log(layer_1@A)),
    "[Oo]verwriting"
  )
  expect_equal(get_latent(overwritten), get_latent(by_args))

  # multi-status definition: one row per status, identical otherwise
  multi <- hhmodel |>
    add_latent("hidden", func = log(layer_1@A), status = c(1, 3)) |>
    get_latent("hidden")
  expect_equal(multi$status, c("1", "3"))
  multi$status <- NULL
  expect_equal(multi[1, , drop = TRUE], multi[2, , drop = TRUE])
})

test_that("Symbolic definitions are equivalent", {
  hhmodel <- mk_latent$model |>
    add_latent("layer_1", model = "uniform", min = 0, max = 2)

  # rlang double bang injects external values
  external_value <- 33
  hhmodel_2 <- hhmodel |>
    add_latent(
      hidden = !!external_value,
      hidden_2 = log(layer_1@A + !!external_value),
      hidden_3 = c(
        model = "uniform",
        min = !!external_value, max = log(layer_1@A + !!external_value)
      ),
      hidden_4 = list(
        model = "uniform",
        min = !!external_value, max = log(layer_1@A + !!external_value)
      )
    )

  get_latent(hhmodel_2, "hidden")$value |>
    expect_equal(external_value)
  get_latent(hhmodel_2, "hidden_2")$func |>
    expect_equal("log")
  get_latent(hhmodel_2, "hidden_2")$func_args[[1]] |>
    expect_equal(c(x = sprintf("sum|x:x2|[layer_1|A]:[%s|N]", external_value)))

  # the two spellings of the same bounded definition agree
  by_vector <- get_latent(hhmodel_2, "hidden_3")
  by_list <- get_latent(hhmodel_2, "hidden_4")
  by_vector$name <- by_list$name <- NULL
  expect_equal(by_vector, by_list, ignore_attr = TRUE)

  # a bare name, a bang and a formula encode the same read
  val <- 1
  augm <- hhmodel |>
    add_latent(
      hidden = list(model = "normal", mean = layer_1 + 1),
      hidden_2 = list(model = "normal", mean = layer_1 + !!val),
      hidden_3 = list(model = "normal", mean = ~ layer_1 + 1),
      hidden_4 = list(model = "normal", mean = ~ layer_1 + !!val)
    ) |>
    get_latent() |>
    subset(name != "layer_1")
  augm$name <- NULL
  # every row is identical
  expect_equal(nrow(unique(augm)), 1)
})

test_that("An external list of latent definitions can be imported", {
  # an external list is loaded through 'import', not as a bare definition
  external <- list(
    hidden = list(model = "gamma", shape = ~hidden_2, scale = ~hidden_3)
  )
  hhmodel <- mk_latent$model |>
    add_latent(hidden_2 = list(model = "uniform", min = 1, max = 2),
               hidden_3 = list(model = "uniform", min = 1, max = 2))
  hhmodel |>
    add_latent(external, append = FALSE) |>
    expect_error("non.?empty")

  augm <- hhmodel |>
    add_latent(import = external) |>
    get_latent("hidden")
  expect_equal(augm$model, "gamma")
  expect_equal(augm$args[[1]], c(shape = "hidden_2|A", scale = "hidden_3|A"))

  # formulas import as the bare expressions they wrap
  mk_latent$model |>
    add_latent(import = list(layer_1 = ~ datum@D,
                             hidden = ~ sqrt(log(layer_1)))) |>
    get_latent() |>
    expect_equal(get_latent(add_latent(mk_latent$model, layer_1 = datum@D,
                                       hidden = sqrt(log(layer_1)))))
})

test_that("Latent template import loads", {
  hhmodel <- mk_outbreak$model |>
    add_latent(hidden = start_date@D + 1, append = FALSE)

  # appending keeps the pre-existing definitions beside the template ones
  hhmodel |>
    add_latent(import = "template") |>
    get_latent() |>
    getElement("name") |>
    expect_setequal(c("hidden", "infec_time"))

  # dropping resets the store first
  hhmodel |>
    add_latent(import = "template", append = FALSE) |>
    get_latent() |>
    getElement("name") |>
    expect_setequal("infec_time")
})

test_that("Nested latent definitions raise errors", {
  hhmodel <- data.frame(hhid = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("test")

  hhmodel |>
    add_latent(layer_2 = list(model = "uniform"),
               layer_1 = list(model = "uniform", status = 1),
               layer_1 = list(value = layer_2, status = c(2, 3, 4)),
               hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_no_error()

  # simple nest
  hhmodel |>
    add_latent(layer_1 = hidden, hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic latent definition")

  # triangular nest
  hhmodel |>
    add_latent(hidden = layer_1, layer_1 = layer_2, layer_2 = hidden) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic latent definition")

  # first generation of definitions: no layer_1 for status 3 and 4
  hhmodel |>
    add_latent(layer_2 = list(model = "lognormal", status = 2),
               layer_1 = list(model = "normal", status = 1),
               layer_1 = list(value = layer_2, status = 2),
               hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")

  # second generation of definitions: no layer_2 for status 4
  hhmodel |>
    add_latent(layer_2 = list(model = "lognormal", status = c(2, 3)),
               layer_1 = list(model = "normal", status = 1),
               layer_1 = list(value = layer_2, status = c(2, 3, 4)),
               hidden = layer_1,
               append = FALSE) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")
})

test_that("A discrete latent resolves 'uniform' to the discrete distribution", {
  hhmodel <- data.frame(hhid = "h1", id = c("m1", "m2"),
                        status = c(1, 1)) |>
    new_housepi_model("test")

  hhmodel |>
    add_latent("hidden", model = "uniform", min = 0, max = 3, discr = TRUE) |>
    evaluate_log_lik() |>
    expect_no_error()

  # the flat lattice law owns every integer value in the window, so a bound that
  # is not one has no reading, and the error names the latent that carries it
  sms <- hhmodel |>
    add_latent("hidden", model = "uniform", min = 0.5, max = 3, discr = TRUE) |>
    evaluate_log_lik() |>
    tryCatch(error = conditionMessage)
  expect_match(sms, "min.*uniform.*discrete")
  expect_match(sms, "hidden")
})

test_that("A location off the integers is refused when discrete", {
  hhmodel <- data.frame(hhid = 1, id = c("m1", "m2"), status = c(1, 0)) |>
    new_housepi_model("test")

  # a sampled location would shift the lattice under the chain at every step:
  # the fix is on the declaration it reads, named in place of `loc`
  hhmodel |>
    add_parameter("rho", prior = "uniform", min = 0, max = 1) |>
    add_latent("hidden", model = "poisson", lambda = 2, loc = rho@P) |>
    evaluate_log_lik() |>
    expect_error("[Cc]hange parameter.*rho.*discrete")
})

test_that("A latent bound to a distribution returns its density, not a draw", {
  # a distribution called as a func gives back a mass: reading its support here
  # would type the column integer and truncate every mass to zero
  simulated <- data.frame(hhid = rep(seq(3), each = 2), id = seq(6),
                          status = 1) |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "uniform", max = 5, discr = TRUE),
               hidden = list(func = "poisson", x = layer_1, lambda = 3)) |>
    run_simulation(nrep = 1)

  expect_type(simulated$layer_1, "integer")
  expect_type(simulated$hidden, "double")
  expect_equal(simulated$hidden, dpois(simulated$layer_1, 3))
})

test_that("Bare column reads inherit the discretization of their column", {
  # sources resolve before `discr` does, so a read carrying no sigil reaches
  # the column and inherits its storage type, exactly as an `@D` read does
  augm <- data.frame(hhid = rep(c("h1", "h2"), each = 2),
                     datum = as.Date("2005-02-28") + rep(0:1, 2),
                     datum_2 = rep(c(3.5, 5.5), 2),
                     status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "normal", value = datum),
               hidden_2 = list(model = "normal", value = datum@D),
               hidden_3 = list(model = "normal", value = datum_2)) |>
    get_latent()

  expect_equal(augm$value, c("datum|D", "datum|D", "datum_2|D"))
  # the date column is discrete, the double one is not
  expect_equal(augm$discr, c(TRUE, TRUE, FALSE))
})

###################################
###  OBSERVATION : DEFINITIONS  ###
###################################

test_that("Observation definitions are encoded", {
  obs <- mk_latent$model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(
      record = list(model = "normal", mean = hidden@A, sd = 1)
    ) |>
    get_observations()

  # marked as observation, data-bound and constant (never sampled)
  expect_identical(obs[["level"]], "member")
  expect_true(obs[["const"]])
  expect_identical(obs[["value"]], "record|D")
  expect_true(is.na(obs[["proposal"]]))
  # distribution carries the observation arguments
  expect_identical(obs[["model"]], "normal")
  expect_equal(obs[["args"]][[1]], c(mean = "hidden|A", sd = "1"))

  # a deterministic observation carries a function and no distribution
  obs <- mk_latent$model |>
    add_latent(layer_1 = list(model = "gamma", shape = 2, scale = 1),
               hidden = -layer_1) |>
    add_observation(record = floor(hidden + layer_1)) |>
    get_observations()

  expect_identical(obs[["level"]], "member")
  expect_true(is.na(obs[["model"]]))
  expect_true(is.na(obs[["value"]]))
  expect_match(obs[["func"]], "floor")
})

test_that("Observation node rejects conflicting options", {
  # observation nodes are data-bound, never sampled: the options are dropped
  ignored <- function(...) {
    expect_warning(hhmodel <- add_observation(mk_latent$model, ...),
                   "[Ii]gnored.*const.*proposal")
    get_observations(hhmodel)
  }
  obs <- ignored("record", model = "normal", const = FALSE,
                 proposal = "forward")
  expect_true(obs[["const"]])
  expect_true(is.na(obs[["proposal"]]))

  # a deterministic observation drops them alike
  obs <- ignored(record_2 = list(func = min(datum_2@D), status = 1,
                                 const = TRUE, proposal = "rw"))
  expect_identical(obs[["func"]], "min")
  expect_false(obs[["const"]])
  expect_true(is.na(obs[["proposal"]]))

  # a decorator carrying no definition at all
  mk_latent$model |>
    add_observation() |>
    expect_error("[Mm]issing or empty observation definitions")
})

test_that("Household observations reject a distribution", {
  err <- "cannot be sampled"

  # the level is never silently downgraded
  mk_latent$model |>
    add_observation(record_2 = list(model = "normal", mean = 0, sd = 1),
                    level = "house") |>
    expect_error(err)
})

test_that("Household observations reject a status", {
  hhmodel <- mk_latent$model
  err <- "cannot be status specific"

  # guarded at the bulk default level and per definition
  add_observation(hhmodel,
    record_2 = list(func = min(datum_2@D), status = 1), level = "house"
  ) |>
    expect_error(err)
  add_observation(hhmodel,
    record_2 = list(func = min(datum_2@D), status = 1, level = "house")
  ) |>
    expect_error(err)
  add_observation(hhmodel, record_2 = list(func = min(datum_2@D)),
                  status = 1, level = "house") |>
    expect_error(err)

  # a member-level declaration keeps its status, and the no-effect reduction
  # that would promote it to household level is not taken when one is set
  member <- hhmodel |>
    add_observation(record_2 = list(func = min(datum_2@D), status = 1)) |>
    get_observations()
  expect_identical(member[["level"]], "member")
  expect_identical(as.integer(member[["status"]]), 1L)

  # a bulk 'status =' reaches every entry, single and bulk modes alike
  bulk <- hhmodel |>
    add_observation(record_2 = list(func = min(datum_2@D)), status = 1) |>
    get_observations()
  expect_identical(as.integer(bulk[["status"]]), 1L)
  imported <- hhmodel |>
    add_observation(
      import = list(record_2 = list(func = quote(min(datum_2@D)))),
      status = 1
    ) |>
    get_observations()
  expect_identical(as.integer(imported[["status"]]), 1L)
})

test_that("A name is declared under a single kind", {
  hhmodel <- mk_latent$model
  err <- "already declared as"

  # latent and observations share one store: every pair collides, either way
  hhmodel |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(hidden = list(func = datum_2@D)) |>
    expect_error(err)
  hhmodel |>
    add_observation(hidden = list(func = datum_2@D)) |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    expect_error(err)
  hhmodel |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(hidden = min(datum_2@D), level = "house") |>
    expect_error(err)
  hhmodel |>
    add_observation(hidden = min(datum_2@D), level = "house") |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    expect_error(err)
  hhmodel |>
    add_observation(hidden = list(func = datum_2@D)) |>
    add_observation(hidden = min(datum_2@D), level = "house") |>
    expect_error(err)

  # a status-specific declaration is no loophole: the store is keyed by name
  hhmodel |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1, status = 1)) |>
    add_observation(hidden = min(datum_2@D), level = "house") |>
    expect_error(err)

  # the reserved status node is named by add_status, never by the user
  hhmodel |>
    add_latent(.status = list(value = 1, const = TRUE)) |>
    expect_error("cannot start with a dot")

  # siblings of one kind stay legal: statuses of the same name, and overwrites
  hhmodel |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1, status = 1)) |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1, status = 2)) |>
    get_latent() |>
    nrow() |>
    expect_equal(2L)
  hhmodel |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 2)) |>
    expect_warning("Overwriting")

  # an observation standing in for a like-named data column is the surrogate
  # design, not a clash: only declarations collide with each other
  hhmodel |>
    add_observation(datum_2 = list(func = record@D)) |>
    get_observations() |>
    getElement("name") |>
    expect_identical("datum_2")
})

test_that("Level is set per definition, the argument being the bulk default", {
  level_of <- function(obs) {
    stats::setNames(obs[["level"]], obs[["name"]])[c("record_2", "record_3")]
  }
  expected <- c(record_2 = "house", record_3 = "member")

  # a per-definition 'level' picks the kind of its own entry, exactly as
  # 'status' does, and never reaches the definition arguments
  mixed <- mk_latent$model |>
    add_observation(
      record_2 = list(func = min(datum_2@D), level = "house"),
      record_3 = list(func = 21, level = "member")
    ) |>
    get_observations()
  expect_identical(level_of(mixed), expected)
  expect_false("args" %in% names(mixed))

  # the argument is the default: entries not setting their own inherit it
  mk_latent$model |>
    add_observation(
      record_2 = list(func = min(datum_2@D)),
      record_3 = list(func = 21, level = "member"),
      level = "house"
    ) |>
    get_observations() |>
    level_of() |>
    expect_identical(expected)

  # imported definitions honour it too
  mk_latent$model |>
    add_observation(import = list(
      record_2 = list(func = quote(min(datum_2@D)), level = "house"),
      record_3 = list(func = 21)
    )) |>
    get_observations() |>
    level_of() |>
    expect_identical(expected)

  # an unknown level is rejected per definition
  mk_latent$model |>
    add_observation(record_2 = list(func = 1, level = "unknown")) |>
    expect_error("should be one of")

  # 'level' is meaningless outside observations
  mk_latent$model |>
    add_latent(hidden = list(value = 1, level = "house")) |>
    expect_error("reserved to observations")
})

test_that("Observation and latent loads each keep the other's definitions", {
  hhmodel <- mk_latent$model |>
    add_latent(hidden = list(model = "gamma", shape = 2, scale = 1)) |>
    add_observation(record = list(model = "normal", mean = hidden@A, sd = 1))

  # non-append household observation keeps latent definitions
  hhmodel <- add_observation(
    hhmodel, record_2 = min(datum_2@D), level = "house", append = FALSE
  )
  expect_identical(get_latent(hhmodel)[["name"]], "hidden")
  obs <- get_observations(hhmodel)

  # observations are a single family: the member row went with the reset
  expect_setequal(obs[["name"]], "record_2")
  expect_identical(obs[["level"]], "house")
  expect_identical(obs[["func"]], "min")
  expect_identical(obs[["func_args"]][[1]], c(x = "datum_2|D"))

  # appending keeps both levels side by side
  hhmodel <- add_observation(
    hhmodel, record_3 = max(datum_2@D), level = "house"
  )
  expect_setequal(get_observations(hhmodel)[["name"]],
                  c("record_2", "record_3"))
  expect_identical(get_latent(hhmodel)[["name"]], "hidden")

  # symmetrically, a non-append latent load keeps observations
  hhmodel <- add_latent(hhmodel, append = FALSE,
                        layer_1 = list(model = "gamma", shape = 2, scale = 1))
  expect_identical(get_latent(hhmodel)[["name"]], "layer_1")
  expect_setequal(get_observations(hhmodel)[["name"]],
                  c("record_2", "record_3"))
})

test_that("A degenerate reduction at member level is cast to house level", {
  # member-level reduction by mistake
  expect_warning(
    promoted <- mk_latent$model |>
      add_observation(record_2 = min(datum_2@D)) |>
      get_observations(),
    regexp = "no-effect.*reduction.*min"
  )
  expect_identical(promoted[["level"]], "house")
  # explicit house-level reduction
  mk_latent$model |>
    add_observation(record_2 = min(datum_2@D), level = "house") |>
    get_observations() |>
    expect_identical(promoted)

  # redundant reductions raise an error and are not cast
  mk_latent$model |>
    add_observation(record_2 = min(datum_2@D), level = "house") |>
    add_observation(record_3 = min(record_2@D)) |>
    expect_error("redundant.*record_3")

  # genuine per-member expressions pass normally
  mk_latent$model |>
    add_observation(record_2 = min(datum_2@D, 0),
                    record_3 = datum_2@D,
                    record_4 = log(datum_2@D)) |>
    get_observations() |>
    getElement("level") |>
    expect_identical(rep("member", 3))
})

test_that("A missing observation leaves the likelihood", {
  # an observation value is a datum: where the cell is empty there is nothing
  # to condition on, so the term drops out instead of entering as NA
  mk <- function(x) {
    data.frame(id = c("m1", "m2"), hhid = "h1", status = 1, record = x) |>
      new_housepi_model("test") |>
      add_latent("hidden", value = 0.5) |>
      add_observation("record", model = "uniform", loc = hidden@A,
                      min = -0.4, max = 0.4)
  }

  # both members observed: two terms of a uniform over a window of 0.8
  expect_equal(evaluate_log_lik(mk(c(0.5, 0.5))), 2 * log(1 / 0.8))

  # one datum absent: one term left, and the run says which column and how many
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA))),
                 "[Mm]issing for 1 of 2 infected individuals")
  expect_equal(log_lik, log(1 / 0.8))

  # the domain still binds where the datum is there
  expect_equal(evaluate_log_lik(mk(c(0.5, 9))), -Inf)

  # no datum anywhere: the model declares a density nothing feeds
  mk(c(NA_real_, NA_real_)) |>
    evaluate_log_lik() |>
    expect_error("missing for every infected individual")

  # a categorical column reaches the engine as integer codes: its empty cell
  # must code to no level, or the density would read the code of a non-value
  mk_categorical <- function(x) {
    data.frame(id = c("m1", "m2"), hhid = "h1", status = 1, record = x) |>
      new_housepi_model("test") |>
      add_observation("record", model = "uniform", min = 0, max = 3)
  }
  expect_warning(log_lik <- evaluate_log_lik(mk_categorical(c(TRUE, NA))),
                 "[Mm]issing for 1 of 2 infected individuals")
  expect_equal(log_lik, log(1 / 3))
  mk_categorical(c(NA, NA)) |>
    evaluate_log_lik() |>
    expect_error("missing for every infected individual")

  # an observation is judged over the members carrying it: the status-1 member
  # never declares a status-2 observation, so its empty cell is not missing
  data.frame(id = c("m1", "m2"), hhid = "h1",
             status = c(1, 2), record = c(NA, 0.5)) |>
    new_housepi_model("test") |>
    add_latent("hidden", value = 0.5) |>
    add_observation(record = list(model = "uniform", loc = hidden@A,
                                  min = -0.4, max = 0.4, status = 2)) |>
    evaluate_log_lik() |>
    expect_equal(log(1 / 0.8))
})

test_that("A latent goes with the observation that reads it", {
  # a latent nothing else reads is an unused parameter once its observation has
  # left the likelihood: marginalizing it integrates its model to one, so its
  # density goes with the datum. A flat model keeps the total free of the draw.
  mk <- function(x, x_2 = NULL, hop = FALSE) {
    hhdata <- data.frame(id = c("m1", "m2"), hhid = "h1",
                         status = 1, record = x)
    if (!is.null(x_2)) hhdata$record_2 <- x_2
    hhmodel <- hhdata |>
      new_housepi_model("test") |>
      add_latent(hidden = list(model = "uniform", min = 0, max = 2))
    hhmodel <- if (hop) {
      hhmodel |>
        add_latent(hidden_2 = hidden@A) |>
        add_observation("record", model = "uniform", loc = hidden_2@A,
                        min = -5, max = 5)
    } else {
      add_observation(hhmodel, "record", model = "uniform", loc = hidden@A,
                      min = -5, max = 5)
    }
    if (is.null(x_2)) hhmodel
    else add_observation(hhmodel, "record_2", model = "uniform",
                         loc = hidden@A, min = -5, max = 5)
  }

  # both observed: two flat latent densities and two observation terms
  expect_equal(evaluate_log_lik(mk(c(0.5, 0.5))),
               2 * log(1 / 2) + 2 * log(1 / 10))

  # m2 carries no datum: its observation and the latent feeding it both drop
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA))),
                 "[Mm]issing for 1 of 2")
  expect_equal(log_lik, log(1 / 2) + log(1 / 10))

  # the same one hop away, through a deterministic readout
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA), hop = TRUE)),
                 "[Mm]issing for 1 of 2")
  expect_equal(log_lik, log(1 / 2) + log(1 / 10))

  # a second observation still reads it, so m2 keeps its draw
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA), x_2 = c(0.5, 0.5))),
                 "[Mm]issing for 1 of 2")
  expect_equal(log_lik, 2 * log(1 / 2) + 3 * log(1 / 10))

  # simulations generate every observation, so they prune nothing
  simulated <- run_simulation(mk(c(0.5, 0.5)))
  expect_false(anyNA(simulated$record))
})

test_that("Reading an observation in models and loops are rejected", {
  # a latent model must not read an observation (illegal in both
  # modes; caught at eager build, so the construction itself errors)
  data.frame(hhid = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_observation("record", model = "uniform", min = 0, max = 1) |>
    add_latent("hidden", model = "normal", mean = record@D) |>
    expect_error("[Mm]odel.*reads.*observation")

  # a deterministic read that closes a loop (the observation's support
  # depends on the reader) is circular
  data.frame(hhid = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_latent("hidden", value = record@D) |>
    add_observation("record", model = "uniform", min = hidden@A, max = 10) |>
    expect_error("[Cc]yclic")
})

test_that("A simulation override (value + model) rejects observation reads", {
  hhmodel <- data.frame(hhid = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_observation("record", model = "uniform", min = 0, max = 1)

  # the override is applied before the observation is generated
  hhmodel |>
    add_latent("hidden", model = "normal", value = record@D) |>
    expect_error("[Ss]imulation override.*observation")

  # a deterministic readout does not make the observation available any earlier,
  # so the same read one hop away is rejected too
  hhmodel |>
    add_latent(layer_1 = record@D,
               hidden = list(model = "normal", value = layer_1@A)) |>
    expect_error("[Ss]imulation override.*record.*through.*layer_1")

  # a latent read is legal, so a downstream one closes an ordinary loop
  data.frame(hhid = "h1", status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "normal", value = hidden_2@A),
               hidden_2 = list(model = "normal", mean = hidden@A, sd = 1)) |>
    expect_error("[Cc]ycl")
})

test_that("A surrogate column reading latent is simulation only", {
  hhmodel <- data.frame(id = seq_len(10), hhid = "h1",
                        status = c(rep(1, 9), 0), record = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record = hidden + 2)

  # the data view has no derived value to report: it keeps the imported datum
  expect_true(all(get_data_wide(hhmodel)[["record"]] == 0))

  # observation depends on sampled latent: not a valid surrogate definition
  evaluate_log_lik(hhmodel) |>
    expect_error("must resolve from data alone")

  # simulation produces it for infected members, following the latent draw
  simulated <- run_simulation(hhmodel)
  infected <- simulated$status > 0
  expect_equal(sum(names(simulated) == "record"), 1)
  expect_equal(simulated$record[infected], simulated$hidden[infected] + 2)
  expect_true(all(is.na(simulated$record[!infected])))
})

test_that("A surrogate observation overwrites the corresponding data column", {
  hhmodel <- data.frame(hhid = "h1", status = c(1, 0, 0),
                        datum_2 = c(2, 6, NA), record_2 = rep(99, 3)) |>
    new_housepi_model("test") |>
    add_observation(
      record_2 = list(value = min(datum_2@D), level = "house"),
      record_3 = datum_2@D + 1
    ) |>
    add_latent(hidden = list(value = record_3@D * 2))

  # the data view reports the derived value before evaluation
  expect_identical(unique(get_data_wide(hhmodel)[["record_2"]]), 2)
  expect_false("record_3" %in% colnames(get_data_wide(hhmodel)))

  # and after evaluation
  expect_no_error(evaluate_log_lik(hhmodel))
  expect_identical(unique(get_data_wide(hhmodel)[["record_2"]]), 2)
  expect_false("record_3" %in% colnames(get_data_wide(hhmodel)))

  # simulated observations are for all members, latent just for infected
  simulated <- run_simulation(hhmodel)
  expect_equal(simulated[["record_3"]], c(3, 7, NA))
  expect_equal(simulated[["hidden"]], c(6, NA, NA))

  # FIXME: the requirement reported state "latent", but it's "observation"
  # a stochastic observation carries no column to read: the datum is reported
  # missing, rather than failing on the key it stands for
  data.frame(hhid = "h1", status = 1, record = 5) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record_2 = list(model = "normal", mean = hidden, sd = 1)) |>
    add_observation(record = record_2@D + 2) |>
    evaluate_log_lik() |>
    expect_error("missing data column")
})

test_that("Surrogate columns of both levels build in dependency order", {
  hhdata <- data.frame(hhid = "h1", status = c(1, 0, 0),
                       datum = c(2, 9, 4), record = 0, record_2 = 0)
  hhmodel <- hhdata |>
    new_housepi_model("test") |>
    add_observation(record = datum@D + 3,
                    record_2 = list(value = max(record@D), level = "house"))

  # every member carries a value regardless of status
  hhdata_2 <- get_data_wide(hhmodel)
  expect_identical(hhdata_2[["record"]], c(5, 12, 7))
  expect_identical(unique(hhdata_2[["record_2"]]), 12)

  # pure surrogate observations with no raw data columns are not exported
  hhmodel |>
    add_data_wide(hhdata[c("hhid", "status", "datum")]) |>
    get_data_wide() |>
    colnames() |>
    intersect(c("record", "record_2")) |>
    expect_length(0)
})

test_that("A household observation with no datum yet reads as NA", {
  # 'start_date' reduces over an observation the run has still to draw: there
  # is nothing to reduce, so the data view reports NA rather than reading the
  # empty reduction as a datum failing the discreteness 'case_ascertained' asks
  hhdata <- data.frame(
    id = c("m1", "m2", "m3", "m4"), hhid = c("h1", "h1", "h2", "h2"),
    status = NA, end_date = 21
  ) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_observation(
      onset_date = list(model = "normal", discr = TRUE, mean = 3, sd = 1),
      start_date = list(value = min(onset_date@D), level = "house")
    ) |>
    get_data_wide()

  expect_true(all(is.na(hhdata$start_date)))
})

test_that("A surrogate column proves the requirement its model asks", {
  # a surrogate filling 'start_date' answers the discreteness 'case_ascertained'
  # asks of the column, and answers for what it can return, not for the value it
  # happens to settle on. A household-level definition reduces over its members,
  # so its root is the reduction; the member level takes the algebra directly
  mk <- function(defin, level = "house") {
    data.frame(
      id = c("m1", "m2", "m3", "m4"), hhid = c("h1", "h1", "h2", "h2"),
      status = c(1, 0, 1, 0), end_date = 16, datum = 0.5
    ) |>
      new_housepi_model("case-ascertained", quiet = TRUE) |>
      add_observation("start_date", value = !!defin, level = level) |>
      get_data_wide()
  }

  # integer-valued cells fold to the integer lattice through the algebra
  expect_no_error(mk(quote(min(end_date@D))))
  expect_no_error(mk(quote(sum(end_date@D))))
  expect_no_error(mk(quote(end_date@D * 2 - 1), level = "member"))

  # two halves sum to an integer, and are refused all the same: the column
  # answers for what the definition can return, not for the value it settles on
  mk(quote(sum(datum@D))) |>
    expect_error("not discrete-valued")
  mk(quote(end_date@D / 4), level = "member") |>
    expect_error("not discrete-valued")
  mk(quote(sqrt(end_date@D)), level = "member") |>
    expect_error("not discrete-valued")

  # rounding restores the proof
  expect_no_error(mk(quote(round(end_date@D / 4)), level = "member"))
  expect_no_error(mk(quote(ceil(sqrt(end_date@D))), level = "member"))

  # the refusal names the definition and the place it was built for
  sms <- tryCatch(mk(quote(sum(datum@D))), error = conditionMessage)
  expect_match(sms, "[Hh]ousehold observation")
  expect_match(sms, "h1")
  sms <- tryCatch(mk(quote(end_date@D / 4), level = "member"),
                  error = conditionMessage)
  expect_match(sms, "individual")
  expect_match(sms, "m1")
})

test_that("Data-bound household observation rejects latent dependence", {
  # a household observation standing in for data resolves from data alone
  mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -start_date)) |>
    add_observation(start_date = min(test_date@D) + alpha@P,
                    level = "house") |>
    evaluate_log_lik() |>
    expect_error("must resolve from data alone")
})

############################
## STATUS : DEFINITIONS   ##
############################

test_that("A status process needs probabilities, within the status storage", {
  mk_latent$model |>
    add_status() |>
    expect_error("[Nn]o status probability")

  # the engine stores at most 254 infected statuses
  weights <- function(k) {
    stats::setNames(as.list(rep(1 / k, k)), paste0("unit_", seq_len(k)))
  }
  do.call(add_status, c(list(mk_latent$model), weights(254))) |>
    expect_no_error()
  do.call(add_status, c(list(mk_latent$model), weights(255))) |>
    expect_error("255 statuses.*maximum of 254")
})

test_that("Status process rejects latent or parameter probabilities", {
  # probabilities must be a pure function of data/numerics
  mk_status$model |>
    add_status(symptomatic = alpha@P, asymptomatic = 3) |>
    evaluate_log_lik() |>
    expect_error("only on data or numeric")
})

test_that("Escaped label colliding with an infected status is rejected", {
  mk_status$data |>
    new_housepi_model("case-ascertained", quiet = TRUE,
                      escaped = "symptomatic") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    expect_error("collides with the escaped-status label")
})

test_that("Mixing integer and character status declarations is rejected", {
  # one call carries both an integer and a character status: caught at the wall
  mk_status$model |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1,
                          status = "symptomatic"),
      onset = list(value = 1, status = 1)
    ) |>
    evaluate_log_lik() |>
    expect_error("[Mm]ix")
})

test_that("Status is a constant in inference", {
  # refined statuses already in the data; infec_time is read from data so the
  # density is deterministic and the only difference is the status process.
  hhmodel <- helper_onset_house$data |>
    transform(status = c(1, 2, 1, 0, 0, 0, 0)) |>
    add_data_wide(hhmodel = mk_pinned$model)

  # the status process drops out of the density (status is fixed data)
  hhmodel |>
    add_status(symptomatic = 0.7, asymptomatic = 0.3) |>
    evaluate_log_lik() |>
    expect_equal(evaluate_log_lik(hhmodel))
})

###########################
## LATENT : END-TO-END   ##
###########################

test_that("An undeclared producer defers the build to the wall", {
  # parameters land in the engine on their own decorator, so a definition
  # reading one that is not declared yet waits instead of failing a key lookup
  hhmodel <- mk_latent$model |>
    add_latent(hidden = list(model = "gamma", shape = rho@P, scale = rho_2@P))

  # the wall names the missing parameters, not the key that could not be found
  sms <- tryCatch(evaluate_log_lik(hhmodel), error = conditionMessage)
  expect_match(sms, "missing parameter definitions")
  expect_false(grepl("status-specific", sms))
})

test_that("Independence proposal inherits the model location shift", {
  results <- data.frame(id = c("m1", "m2"), hhid = "h1",
                        status = 1, datum = c(4, 8), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", loc = datum@D,
               proposal = "gamma", proposal_shape = 2) |>
    helper_flat_sink(loc = hidden@A) |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  draws <- as.data.frame(get_latent_draws(results))
  expect_true(all(draws$m1_hidden >= 4 & draws$m1_hidden <= 5))
  expect_true(all(draws$m2_hidden >= 8 & draws$m2_hidden <= 9))
})

test_that("Truncated latent models keep their normalizer in the log-likelihood trace", {
  # bounds read data, constants and a sampled parameter
  hhdata <- data.frame(id = c("m1", "m2"), hhid = c("h1", "h2"),
                       status = 1, datum = c(3, 6), obs = 0)
  results <- hhdata |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform", min = 1, max = 3)) |>
    add_latent(layer_1 = list(model = "gamma", shape = 2, scale = 1.5, max = datum@D),
               hidden = list(model = "gamma", shape = rho@P, scale = 1, min = 0.5, max = 8),
               hidden_2 = list(model = "gamma", shape = 1.7, scale = 1, max = 5)) |>
    helper_flat_sink(loc = layer_1@A + hidden@A + hidden_2@A) |>
    run_inference(steps = 50L, adapt = 0L, thin = 1L, chains = 1L,
                  id_draws = c("m1", "m2")) |>
    suppressWarnings()

  draws <- latent_frame(results)
  rho <- as.data.frame(posterior::as_draws_matrix(get_parameter_draws(results)))$rho
  log_lik <- as.data.frame(posterior::as_draws_matrix(get_log_lik_draws(results)))
  log_dtrunc <- function(x, shape, scale, min, max) {
    dgamma(x, shape, scale = scale, log = TRUE) -
      log(pgamma(max, shape, scale = scale) - pgamma(min, shape, scale = scale))
  }
  # one member per household, the flat sink adds log(1 / 2e6)
  for (k in seq(2)) {
    col <- function(name) draws[[sprintf("m%d_%s", k, name)]]
    expected <- log_dtrunc(col("layer_1"), 2, 1.5, 0, hhdata$datum[k]) +
      log_dtrunc(col("hidden"), rho, 1, 0.5, 8) +
      log_dtrunc(col("hidden_2"), 1.7, 1, 0, 5) - log(2e6)
    expect_equal(log_lik[[sprintf("h%d", k)]], expected, tolerance = 1e-8)
  }
})

test_that("Latent draws stay inside their live support", {
  hhdata <- data.frame(id = seq(2), hhid = "h1", status = 1, obs = 0)
  run_short <- function(hhmodel) {
    hhmodel |>
      run_inference(steps = 100L, adapt = 0L, thin = 1L, chains = 1L) |>
      suppressWarnings()
  }

  # other latent variable
  draws <- hhdata |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "uniform"),
               hidden = list(model = "uniform", loc = layer_1@A)) |>
    helper_flat_sink(loc = hidden@A) |>
    run_short() |>
    latent_frame()
  gap <- draws[["1_hidden"]] - draws[["1_layer_1"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # transitive deterministic latent variable
  draws <- hhdata |>
    new_housepi_model("test") |>
    add_latent(layer_2 = list(model = "uniform"),
               layer_1 = layer_2,
               hidden = list(model = "uniform", loc = layer_1@A)) |>
    helper_flat_sink(loc = hidden@A) |>
    run_short() |>
    latent_frame()
  gap <- draws[["1_hidden"]] - draws[["1_layer_2"]]
  expect_true(all(gap >= 0 & gap <= 1))

  # parameter
  results <- hhdata |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform") |>
    add_latent("hidden", model = "uniform", min = rho@P, max = 5) |>
    helper_flat_sink(loc = hidden@A) |>
    run_short()
  hidden <- latent_frame(results)[["1_hidden"]]
  rho <- as.data.frame(get_parameter_draws(results))[["rho"]]
  expect_true(all(hidden >= rho))
})

test_that("Draws requested where no latent is sampled are ignored", {
  # a member carrying no infected status draws no latent
  data.frame(id = c("m1", "m2"), hhid = "h1", status = c(1, 0), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    helper_flat_sink(loc = hidden@A) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L,
                  id_draws = c("m1", "m2")) |>
    expect_warning("[Ii]gnored draws of non infected individual .m2.") |>
    suppressWarnings()

  # a model sampling no latent at all has nothing to trace
  mk_latent$model |>
    add_parameter(rho = list(prior = "uniform")) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, id_draws = "m1") |>
    expect_warning("[Ii]gnored.*id_draws.*no latent variables") |>
    suppressWarnings()
})

test_that("Latent definitions are integrated in simulations", {
  # the downstream latent is declared first: declaration order does not matter
  simulated <- data.frame(id = seq(1, 100),
                          hhid = rep(seq(10), each = 10),
                          status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform",
                             min = layer_1, max = layer_1 + 1),
               layer_1 = list(model = "gamma", shape = 2, scale = 1)) |>
    run_simulation()
  expect_true(all(simulated$layer_1 <= simulated$hidden &
                    simulated$hidden <= simulated$layer_1 + 1))
})

test_that("A value paired with model overrides the draw in simulations only", {
  hhmodel <- data.frame(id = seq_len(100), hhid = 1,
                        status = 1, datum = 7, obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", max = 100, value = datum@D) |>
    helper_flat_sink(loc = hidden@A)

  # simulation: pins the sampled variable to user value
  expect_true(all(run_simulation(hhmodel)$hidden == 7))

  # inference: the value does not pin the sampler and draws span the prior
  draws <- hhmodel |>
    run_inference(steps = 40L, chains = 1L) |>
    get_latent_draws() |>
    suppressWarnings()
  expect_false(all(draws[["1_hidden"]] == 7))

  # a pin composed on another latent binds live: the relation holds member by
  # member, and follows the draw rather than collapsing to a single value
  simulated <- data.frame(id = seq_len(100), hhid = 1,
                          status = 1, datum = 7) |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "gamma", shape = 2, scale = 1),
               hidden = list(model = "uniform", max = 100,
                             value = datum@D - layer_1@A)) |>
    run_simulation()
  expect_equal(simulated$hidden, 7 - simulated$layer_1)
  expect_gt(length(unique(simulated$layer_1)), 1)
})

test_that("Latent traces export only relevant variables", {
  # a constant carries no draw, an implicit function cannot be called on a
  # value, and what no sink reads is not sampled at all
  hhmodel <- data.frame(id = c("m1", "m2"), hhid = "h1",
                        status = 1, record = 1) |>
    new_housepi_model(model = "test") |>
    add_parameter(rho = 1, rho_2 = list(prior = "uniform")) |>
    add_latent(
      # constants
      blocked = 42,
      # sampled variables
      hidden = list(model = "uniform"),
      hidden_2 = list(model = "normal"),
      # public non constant function
      layer = hidden_2 - hidden,
      # equal to parameters
      rho = rho@P,
      rho_sum = rho@P + rho_2@P,
      # equal to other latent objects
      dummy = blocked,
      alias = layer,
      # implicit function
      impl = list(func = "gamma", shape = 3),
      # constant function
      invar = list(func = "gamma", x = dummy, shape = 3),
      # function
      expl = list(func = "gamma", x = alias, shape = 3)
    )
  observed <- add_observation(hhmodel, "record", model = "normal",
                              mean = expl@A + rho_sum@A, sd = 1)

  traced <- function(hhmodel) {
    hhmodel |>
      run_inference(steps = 100L, burn = 0L, thin = 1L, adapt = 0L,
                    chains = 1L, id_draws = c("m1", "m2")) |>
      suppressWarnings() |>
      get_latent_draws() |>
      suppressWarnings() |>
      posterior::variables()
  }

  # the observation reads expl and rho_sum, so their producers are drawn too
  expect_setequal(
    traced(observed),
    paste(rep(c("m1", "m2"), each = 6),
          c("hidden", "hidden_2", "layer", "rho_sum", "alias", "expl"),
          sep = "_")
  )

  # with no observation nothing reads any of them: the model has no sink
  expect_length(traced(hhmodel), 0)
})

test_that("A latent degenerate domain warns and collapses", {
  hhmodel <- mk_latent$model |>
    add_parameter("rho", prior = "uniform") |>
    add_latent("hidden", model = "uniform", min = 5, max = 5)
  # steps meet the diagnostic threshold, so the collapse warning stands alone
  expect_warning(
    tryCatch(
      run_inference(hhmodel, steps = 400L, chains = 1L, adapt = 0L, thin = 1L),
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

  # a flat observation on each chain leaves both posteriors equal to the model
  results <- hhmodel |>
    add_data_wide(transform(hhmodel$data_wide, obs = 0, obs_norm = 0)) |>
    helper_flat_sink(loc = z0@A) |>
    helper_flat_sink(loc = norm@A, column = "obs_norm") |>
    run_inference(steps = 1e5L, chains = 1L, adapt = 0L, thin = 1L)
  # no parameters nor priors, only latent
  expect_equal(posterior::nvariables(results$parameter_draws), 0L)
  posterior::nvariables(results$latent_draws) |>
    expect_equal(dim(hhmodel$latent)[1])

  # the S3 method of summarise_draws dispatches on the results
  out_vars <- posterior::summarise_draws(results)[["variable"]]
  expect_true(all(startsWith(out_vars, "tommy")))

  # empirical tests for goodness of Gaussian distribution
  probs <- c(0.025, 0.5, 0.975)
  for (draws in list(results$latent_draws$tommy_norm,
                     results$latent_draws$tommy_z0)) {
    expect_equal(mean(draws), 0, tolerance = 0.05)
    expect_equal(sd(draws), 1, tolerance = 0.05)
    expect_equal(unname(stats::quantile(draws, probs)), qnorm(probs),
                 tolerance = 0.05)
  }
})

test_that("Status-specific latent variables draw from per-status ranges", {
  hhmodel <- data.frame(hhid = 1, id = seq(4), status = seq(4), obs = 0) |>
    new_housepi_model("test") |>
    add_latent(
      hidden = list(status = 1, model = "uniform", min = 0, max = 1),
      hidden = list(status = 2, model = "uniform", min = 1, max = 2),
      hidden = list(status = 3, model = "uniform", min = 2, max = 3),
      hidden = list(status = 4, model = "uniform", min = 3, max = 4)
    ) |>
    helper_flat_sink(loc = hidden@A)

  # simulation routes draws into per-status supports
  simulated <- do.call(rbind, run_simulation(hhmodel, nrep = 50))
  for (s in seq(4)) {
    draws <- simulated$hidden[simulated$status == s]
    expect_true(all(draws >= s - 1 & draws <= s))
  }

  # inference traces every member within its status support
  results <- hhmodel |>
    run_inference(steps = 5L, thin = 1L, adapt = 2L,
                  chains = 1L, id_draws = hhmodel$data_wide$id) |>
    suppressWarnings()
  draws <- as.data.frame(get_latent_draws(results))
  for (s in seq(4)) {
    trace <- draws[[paste0(s, "_hidden")]]
    expect_length(trace, 5L)
    expect_true(all(trace >= s - 1 & trace <= s))
  }
})

test_that("Evaluating a model allocates its latent storage", {
  hhmodel <- data.frame(hhid = 1, id = c("m1", "m2"), status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform")
  evaluate_log_lik(hhmodel)
  augmented <- inspect_augmented_data(hhmodel)
  expect_true(is.data.frame(augmented))
  expect_true("hidden" %in% colnames(augmented))
})

################################
## OBSERVATION : END-TO-END   ##
################################

test_that("Observation-aware initialization", {
  window <- 0.1
  # latent hidden ~ U(0, 1)
  # observation record ~ U(hidden - w, hidden + w)
  # if record = 0.5, valid for hidden in [0.4, 0.6]
  # without observation-aware restart, it would fail with 80% proba per-chain
  # probability of failing on K chains: 1 - 0.2^K
  hhmodel <- data.frame(hhid = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    add_observation("record", model = "uniform", loc = hidden,
                    min = -!!window, max = !!window)
  hhmodel |>
    run_inference(steps = 10L, adapt = 0L, chains = 40L) |>
    suppressWarnings() |>
    expect_no_error()

  # datum record = 5 is unreachable
  hhdata <- data.frame(hhid = "h1", status = 1, record = 5)
  hhmodel |>
    add_data_wide(hhdata) |>
    run_inference() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  hhdata |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform") |>
    add_observation("record", model = "uniform", min = -rho@P, max = rho@P) |>
    run_inference() |>
    expect_error("[Ff]ree parameter dependencies")

  # generic error
  hhdata |>
    new_housepi_model("test") |>
    add_parameter("rho", value = 0.3) |>
    add_latent("hidden", model = "uniform") |>
    add_observation("record", model = "uniform", min = -rho@P, max = rho@P) |>
    run_inference() |>
    expect_error("[Oo]utside of respective domain")
})

test_that("Members missing different observations still initialize", {
  # m1 carries a record and no record_2, m2 the reverse, m3 neither: the
  # ordinary censoring pattern must reach the sampler, not the
  # domain-feasibility wall of the observation-aware restart
  data.frame(
    id = c("m1", "m2", "m3"), hhid = "h1", status = 1,
    record = c(0.5, NA, NA), record_2 = c(NA, 0.5, NA)
  ) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    add_observation(
      record = list(model = "uniform", loc = hidden, min = -0.1, max = 0.1),
      record_2 = list(model = "uniform", loc = hidden, min = -0.1, max = 0.1)
    ) |>
    run_inference(steps = 10L, adapt = 0L, chains = 2L) |>
    suppressWarnings() |>
    expect_no_error()
})

test_that("An observation infers discretization from the column it reads", {
  hhmodel <- data.frame(hhid = rep(c("h1", "h2"), each = 20),
                        date_col = as.Date("2005-02-28") + rep(0:1, 20),
                        int_col = rep(c(3L, 5L), 20),
                        whole_col = rep(c(3, 5), 20),
                        kept_col = rep(c(3L, 5L), 20),
                        status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      date_col = list(model = "normal", mean = hidden, sd = 2),
      int_col = list(model = "normal", mean = hidden, sd = 2),
      # a double column is always continuous
      whole_col = list(model = "normal", mean = hidden, sd = 2),
      # an explicit declaration is never overridden by the data
      kept_col = list(model = "normal", discr = FALSE, mean = hidden, sd = 2)
    )

  simulated <- run_simulation(hhmodel)

  # dates and integers are discrete, doubles are not
  expect_true(on_integer_lattice(simulated[["date_col"]]))
  expect_true(on_integer_lattice(simulated[["int_col"]]))
  expect_false(on_integer_lattice(simulated[["whole_col"]]))
  expect_false(on_integer_lattice(simulated[["kept_col"]]))

  # a deterministic observation is not influenced
  hhmodel |>
    add_observation(record = hidden + 1) |>
    run_simulation() |>
    getElement("record") |>
    on_integer_lattice() |>
    expect_false()
})

test_that("Simulated columns carry the discretization of their definition", {
  hhmodel <- data.frame(hhid = rep(c("h1", "h2"), each = 3),
                        id = seq(6), status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      record = list(model = "normal", discr = TRUE, mean = hidden, sd = 1),
      record_2 = list(model = "normal", mean = hidden, sd = 1),
      record_3 = list(value = min(record@D), level = "house")
    )

  simulated <- run_simulation(hhmodel)

  # a discrete declaration, and the reduction folding it, export as integers
  expect_type(simulated[["record"]], "integer")
  expect_type(simulated[["record_3"]], "integer")
  # a continuous observation and a continuous latent stay double
  expect_type(simulated[["record_2"]], "double")
  expect_type(simulated[["hidden"]], "double")

  # feeding the simulated data back resolves the same discretization
  simulated |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1)) |>
    get_observations() |>
    getElement("discr") |>
    expect_true()
})

test_that("A column a status declares continuous keeps its fractions", {
  simulated <- data.frame(hhid = rep(c("h1", "h2"), each = 2),
                          id = seq(4), status = rep(c(1, 2), 2)) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      record = list(status = 1, model = "normal", discr = TRUE, mean = hidden),
      record = list(status = 2, model = "normal", mean = hidden)
    ) |>
    run_simulation()

  # one column, two declarations: the continuous one has values to lose
  expect_type(simulated[["record"]], "double")
  expect_false(on_integer_lattice(simulated[["record"]][simulated$status == 2]))
})

test_that("Discretization is resolved every time the data changes", {
  hhdata <- data.frame(hhid = rep(c("h1", "h2"), 20), status = 1)
  hhmodel <- new_housepi_model(model = "test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 2))

  # observations resolve to continuous in simulations and definitions
  hhmodel <- hhmodel |>
    add_data_wide(dplyr::mutate(hhdata, record = rep(c(2.4, 5.1), 20)))
  expect_false(get_observations(hhmodel)$discr)
  expect_false(on_integer_lattice(run_simulation(hhmodel)$record))

  # observations resolve to discrete in simulations and definitions
  hhmodel <- hhmodel |>
    add_data_wide(dplyr::mutate(hhdata, record = rep(c(2L, 5L), 20)))
  expect_true(get_observations(hhmodel)$discr)
  expect_true(on_integer_lattice(run_simulation(hhmodel)$record))
})

test_that("Observation cyclic dependencies are captured", {
  data.frame(hhid = "h1", status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "lognormal"),
               hidden_2 = list(model = "uniform", max = record_2@D)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1)) |>
    add_observation("record_2", value = min(record@D), level = "house") |>
    run_simulation() |>
    expect_error("cyclic dependency")
})

test_that("A deterministic observation constrains its ancestors", {
  # an observation whose location is a deterministic node
  # (onset_time = infec_time + incub_period) must propagate the constraint to
  # the sampled ancestors. Otherwise moving infec_time / incub_period skips the
  # observation ratio and the latent onset drifts away from the datum.
  onset <- mk_onset$model |>
    add_latent(onset_time = infec_time + incub_period) |>
    # tightly informative observation reached through the deterministic node
    add_observation(
      "onset_date", model = "normal", mean = onset_time, sd = 0.2
    ) |>
    run_inference(steps = 400L, thin = 1L, adapt = 100L, chains = 1L) |>
    suppressWarnings() |>
    latent_frame()

  # the sampled onset (infec_time + incub_period) tracks the observed onset_date
  recon <- unname(colMeans(onset[, grep("onset_time", names(onset))]))
  observed <- mk_onset$data$onset_date[mk_onset$data$status > 0]
  expect_lt(max(abs(recon - observed)), 0.5)
})

test_that("Discretized observation node initializes finite", {
  # a discretized normal observation whose latent onset starts (from the wide
  # infec_time prior) far below the observed onset_date must not read a
  # spurious -Inf at chain init (upper-tail CDF-difference cancellation)
  results <- mk_onset$model |>
    add_observation("onset_date", model = "normal", discr = TRUE, sd = 1,
                    mean = infec_time + incub_period) |>
    run_inference(steps = 50L, thin = 1L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  expect_true(all(is.finite(get_log_lik_draws(results)$h1)))
})

test_that("A stochastic observation is simulated for infected members only", {
  simulated <- mk_outbreak$model |>
    add_latent(onset_time = infec_time + incub_period) |>
    # one definition set drives simulation too: no hand-written obs_process
    add_observation("onset_date", model = "normal", mean = onset_time) |>
    run_simulation()

  # onset_date is generated for infected members and only those
  infected <- simulated[simulated$status > 0, ]
  expect_false(anyNA(infected$onset_date))
  expect_true(all(is.na(simulated$onset_date[simulated$status == 0])))

  # generated onset_date is a draw from normal(infec_time + incub_period, 1)
  resid <- infected$onset_date - (infected$infec_time + infected$incub_period)
  expect_lt(abs(mean(resid)), 0.2)
  expect_lt(abs(sd(resid) - 1), 0.2)
})

test_that("Household observation derives from simulated observations", {
  # both the observation (test_date) and the surrogate column (start_date) are
  # produced internally: neither is provided in the input data
  reps <- mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -test_date, -start_date)) |>
    add_observation("test_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1) |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    run_simulation(nrep = 2)

  # each repetition is self-consistent: start_date is the per-household na.rm
  # min of that repetition's generated test_date, one value per household
  # replicated across rows, never NA (every house has a case)
  for (simulated in reps) {
    expect_false(anyNA(simulated$start_date))
    per_house <- tapply(simulated$start_date, simulated$hhid, unique)
    expect_true(all(lengths(per_house) == 1))
    infected <- simulated[simulated$status > 0, ]
    expected <- tapply(infected$test_date, infected$hhid, min)
    expect_equal(as.numeric(unlist(per_house)[names(expected)]),
                 as.numeric(expected))
  }
})

test_that("Household observations aggregate a data column across members", {
  # one household carries no test_date at all: its na.rm reduction is empty
  empty_house <- mk_outbreak$data$hhid[[1]]
  hhdata <- dplyr::mutate(
    mk_outbreak$data,
    test_date = dplyr::if_else(hhid == empty_house, as.Date(NA), test_date)
  )

  simulated <- mk_outbreak$model |>
    add_data_wide(hhdata) |>
    add_observation(
      first_test = min(test_date@D),
      last_test = max(test_date@D),
      level = "house"
    ) |>
    run_simulation()

  # one value per household, replicated across every member row
  per_house <- function(col) {
    values <- tapply(simulated[[col]], simulated$hhid, unique)
    expect_true(all(lengths(values) == 1))
    unlist(values)[sort(names(values))]
  }
  # values equal the na.rm min/max of test_date over each household's members,
  # NA when every member is missing
  expected <- function(f) {
    reduce <- function(x) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
    values <- tapply(as.numeric(hhdata$test_date), hhdata$hhid, reduce)
    values[sort(names(values))]
  }
  expect_equal(per_house("first_test"), expected(min))
  expect_equal(per_house("last_test"), expected(max))
  expect_true(is.na(per_house("first_test")[[empty_house]]))
})

test_that("Simulation leaves overwritten data columns pristine", {
  # test_date is provided data and a generated observation: simulation
  # overwrites it with generated values, then must restore the input
  hhmodel <- mk_outbreak$model |>
    add_observation("test_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1)

  before <- inspect_augmented_data(hhmodel)
  invisible(run_simulation(hhmodel, nrep = 2))
  after <- inspect_augmented_data(hhmodel)

  before <- before[order(before$id), ]
  after <- after[order(after$id), ]
  # the overwritten column is restored, and no simulation-added column lingers
  expect_equal(as.numeric(after$test_date), as.numeric(before$test_date))
  expect_setequal(names(after), names(before))
})

test_that("A latent reads a static data-bound household observation", {
  # start_date stands in for provided (static) test_date, read by infec_time
  simulated <- mk_outbreak$model |>
    add_data_wide(mk_outbreak$data |>
                    dplyr::mutate(test_date = 5, end_date = 26) |>
                    dplyr::select(-start_date)) |>
    add_latent(
      infec_time = list(model = "uniform", min = start_date@D - 5,
                        max = start_date@D + 5),
      append = FALSE
    ) |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    run_simulation()

  # start_date is materialized statically (= 5) before the index is drawn
  expect_true(all(simulated$start_date == 5))
})

test_that("Simulation rejects a latent reading a live household observation", {
  # infec_time depends on start_date, a surrogate of the generated test_date:
  # an unsatisfiable cycle, rejected with a clear message
  mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -test_date, -start_date)) |>
    add_latent(
      infec_time = list(model = "uniform", min = start_date@D - 5,
                        max = start_date@D + 5),
      append = FALSE
    ) |>
    add_observation("test_date", model = "normal", mean = infec_time, sd = 1) |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    run_simulation() |>
    expect_error("[Cc]yclic")
})

test_that("Data-bound household observation supplies a required data column", {
  # the dates of the dataset itself, where start_date is the first test date
  hhdata <- sarscov2_2021_icofs |>
    dplyr::filter(hhid %in% utils::head(unique(hhid), 40))
  hhmodel <- new_housepi_model(model = "case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")

  # 'start_date' provided as an ordinary data column
  provided <- add_data_wide(hhmodel, hhdata)

  # 'start_date' is passed as a surrogate observation
  surrogate <- hhmodel |>
    add_data_wide(dplyr::select(hhdata, -start_date)) |>
    add_observation(start_date = min(test_date@D), level = "house")

  # the two agree on the pointwise sum: the unpinned template draws the same
  # values under the same seed
  set.seed(42)
  log_lik <- evaluate_log_lik(provided)
  set.seed(42)
  evaluate_log_lik(surrogate) |>
    expect_equal(log_lik)
})

###########################
## STATUS : END-TO-END   ##
###########################

test_that("Status process generates a labelled status factor in simulation", {
  simulated <- mk_outbreak$model |>
    # categorical drawn first; labels (declaration order) name the statuses
    add_status(symptomatic = 0.2, asymptomatic = 0.8) |>
    run_simulation()

  # status is a factor: code 0 -> escaped, 1 -> symptomatic, 2 -> asymptomatic
  expect_s3_class(simulated$status, "factor")
  expect_identical(levels(simulated$status),
                   c("escaped", "symptomatic", "asymptomatic"))
  infected <- as.character(simulated$status[simulated$status != "escaped"])
  # outcome split follows the probabilities (0.2 symptomatic)
  expect_lt(abs(mean(infected == "symptomatic") - 0.2), 0.08)

  # feed the factor-status output back for inference: factorization maps the
  # labels to the same codes (escaped -> 0), so the density is finite
  hhmodel <- simulated |>
    dplyr::select(id, hhid, status, infec_time, start_date, end_date) |>
    add_data_wide(hhmodel = mk_pinned$model)
  expect_true(is.finite(evaluate_log_lik(hhmodel)))

  # the data-derived dictionary matches the simulated infected labels
  expect_identical(levels(get_data_wide(hhmodel)$status),
                   c("escaped", "symptomatic", "asymptomatic"))
})

test_that("Status process drives status-specific latent blocks", {
  simulated <- mk_outbreak$model |>
    # status-specific block: symptomatic incubation is longer
    add_latent(status = "symptomatic",
               incub_period = list(model = "gamma", shape = 8, scale = 1)) |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  by_status <- tapply(simulated$incub_period,
                      as.character(simulated$status), mean, na.rm = TRUE)
  # the resolved status selects the block: shape 8 symptomatic, 2 asymptomatic
  expect_gt(by_status[["symptomatic"]], 5)
  expect_lt(by_status[["asymptomatic"]], 4)
})

test_that("Escaped label is configurable", {
  simulated <- mk_outbreak$data |>
    new_housepi_model("case-ascertained", escaped = "healthy") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-12, beta = 0.001, gamma = 1, hhsize_ref = 4) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  # weak transmission leaves escapees, labelled per the construction argument
  expect_true("healthy" %in% levels(simulated$status))
  expect_true(any(as.character(simulated$status) == "healthy"))
})

test_that("Status-specific blocks alone cover every drawn status", {
  # simulation draws the status, so the data carries none: the blocks a drawn
  # status can land on are the declared ones, never the observed ones
  mk <- function(..., weights = c(0.5, 0.5)) {
    data.frame(
      id = c("m1", "m2", "m3", "m4"), hhid = c("h1", "h1", "h2", "h2"),
      status = NA, start_date = 0, end_date = 21
    ) |>
      new_housepi_model("case-ascertained", quiet = TRUE) |>
      add_parameter(import = "template") |>
      add_infectivity_profile(import = "template") |>
      add_status(symptomatic = weights[1], asymptomatic = weights[2]) |>
      add_latent(...)
  }

  # every declared status carries its own block, so no generic one is needed
  covered <- mk(
    infec_time = list(value = 0, status = 1),
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14, status = 1),
    onset = list(func = infec_time + incub, status = 1),
    infec_time = list(value = 0, status = 2),
    incub = list(model = "gamma", mean = 3, sd = 1, max = 14, status = 2),
    onset = list(func = infec_time + incub, status = 2)
  )
  expect_s3_class(suppressWarnings(run_simulation(covered)), "data.frame")

  # a status left without a block is named itself, not reported as a missing
  # generic definition the model never had to declare. Every infected member
  # draws that status, so the report does not depend on how the draw falls
  uncovered <- mk(
    infec_time = list(value = 0, status = 1),
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14, status = 1),
    onset = list(func = infec_time + incub, status = 1),
    weights = c(0, 1)
  )
  sms <- tryCatch(suppressWarnings(run_simulation(uncovered)),
                  error = conditionMessage)
  expect_match(sms, "status 2")
  expect_false(grepl("generic latent definition", sms))

  # a generic block on its own still serves every drawn status
  generic <- mk(
    infec_time = 0,
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14),
    onset = list(func = infec_time + incub)
  )
  expect_s3_class(suppressWarnings(run_simulation(generic)), "data.frame")
})
