library(testthat)

# model carrying the data columns the definition blocks refer to
mk_latent <- local({
  hh_model <- data.frame(hh_id = "h1", status = 1, datum = 1,
                         record = 3, datum_2 = 2) |>
    new_housepi_model("test")
  list(model = hh_model)
})

# saturated outbreak: every household infects all of its members
mk_outbreak <- local({
  hhdata <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 50)) |>
    dplyr::mutate(end_date = 21)
  hh_model <- hhdata |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-12, beta = 100, gamma = 1, hh_size_ref = 4) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    )
  list(data = hhdata, model = hh_model)
})

# case-ascertained household sampling incubation periods and infection times
mk_onset <- local({
  hhdata <- dplyr::select(helper_onset_house$data, -infec_time)
  hh_model <- hhdata |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1, max = 14),
      infec_time = list(min = -14, max = end_date)
    ) |>
    add_infectivity_profile(func = "gamma", shape = shape, scale = scale)
  list(data = hhdata, model = hh_model)
})

# two infected statuses under the template definitions
mk_status <- local({
  hhdata <- data.frame(id = c("m1", "m2", "m3", "m4"), hh_id = "h1",
                       status = c(1, 2, 0, 0), end_date = 21)
  hh_model <- hhdata |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")
  list(data = hhdata, model = hh_model)
})

# whether on an integer lattice
on_integer_lattice <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  length(x) > 0 && all(x == floor(x))
}

# draws of a run as a data frame of its latent traces
latent_frame <- \(x) draws_to_frame(get_latent_draws(x))

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

  # strings outside the naming convention
  mk_latent$model |>
    add_latent(hidden = list(model = "normal", loc = "a|b", sd = 1)) |>
    expect_error("[Vv]alue .a.b. contained invalid")
  mk_latent$model |>
    add_latent(hidden = list(model = "normal|x", loc = 0, sd = 1)) |>
    expect_error("[Mm]odel .normal.x. contained invalid")
  mk_latent$model |>
    add_latent(hidden = list(model = "normal", loc = 0, sd = 1,
                             proposal = "rw|x")) |>
    expect_error("[Pp]roposal .rw.x. contained invalid")
  mk_latent$model |>
    add_latent(hidden = list(model = "normal", `l:c` = 0, sd = 1)) |>
    expect_error("[Aa]rgument name .l.c. contained invalid")
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
  defs$proposal |> expect_all_na()
})

test_that("A latent without a model is sampled within its bounds", {
  # bounds alone sample with the default random walk
  latent <- mk_status$model |>
    add_latent(infec_time = list(min = -1, max = end_date), append = FALSE) |>
    get_latent()
  latent$proposal |> expect_identical("rw")
  latent$const |> expect_false()
  latent$model |> expect_all_na()

  # a value next to a proposal pins the simulation only
  latent <- mk_status$model |>
    add_latent(infec_time = list(value = 0, min = -1, max = end_date,
                                 proposal = "rw"),
               append = FALSE) |>
    get_latent()
  latent$value |> expect_identical(0)
  latent$const |> expect_false()

  # a constant value stays constant, whatever else the latent declares
  latent <- mk_status$model |>
    add_latent(infec_time = list(value = 3, const = TRUE, min = 0, max = 5,
                                 proposal = "rw"),
               append = FALSE) |>
    get_latent()
  latent$const |> expect_true()
  latent$proposal |> expect_all_na()
})

test_that("A latent without a model rejects incomplete definitions", {
  # a value with bounds reads either as a constant or as a simulation value
  mk_latent$model |>
    add_latent(hidden = list(value = 0, min = 0, max = 5)) |>
    expect_error("value and bounds, but no model or proposal")
  # a forward proposal draws from the model
  mk_latent$model |>
    add_latent(hidden = list(min = 0, max = 5, proposal = "forward")) |>
    expect_error("[Nn]o .model. to sample from")
  # flat sampling needs both bounds
  mk_latent$model |>
    add_latent(hidden = list(min = 0)) |>
    expect_error("[Nn]o .model. to sample from")
})

test_that("Passing a prior raises a hint on model", {
  sms <- mk_latent$model |>
    add_latent(hidden = list(prior = "uniform", min = -1, max = 2)) |>
    tryCatch(error = conditionMessage)
  sms |> expect_match("[Ss]pare.*prior.*hidden")
  sms |> expect_match("[Dd]id you mean.*model")
})

test_that("A func definition is deterministic", {
  # a plain number is a constant value, whichever argument carried it
  by_func <- mk_latent$model |> add_latent(hidden = list(func = 21))
  by_value <- mk_latent$model |> add_latent(hidden = list(value = 21))
  get_latent(by_func) |> expect_identical(get_latent(by_value))

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
    expect_error("[Oo]nly one of .value. and .func.")
})

test_that("A name is resolved, not guessed", {
  # a name two holders answer to is never read
  sms <- mk_latent$model |>
    add_parameter(datum = list(prior = "gamma", shape = 2, scale = 1)) |>
    add_observation(record = list(model = "normal", mean = datum, sd = 1)) |>
    evaluate_log_lik() |>
    tryCatch(error = conditionMessage)
  sms |> expect_match("[Cc]ould not resolve the source of .datum.")
  sms |> expect_match("[Dd]eclared both as data column and parameter")
  sms |> expect_match("[Rr]ename one of them")

  # nothing of that name anywhere: the defect is named, with the near miss
  sms <- mk_latent$model |>
    add_observation(record = list(model = "normal", mean = datums, sd = 1)) |>
    evaluate_log_lik() |>
    tryCatch(error = conditionMessage)
  sms |> expect_match("[Cc]ould not resolve the source of .datums.")
  sms |> expect_match("[Dd]id you mean .datum.")
})

test_that("A forward consumer waits resolves after declaration", {
  obs <- mk_latent$model |>
    add_observation(record = list(model = "normal", mean = rho, sd = 1)) |>
    get_observations()
  obs[["args"]][[1]] |> expect_identical(c(mean = "rho|?", sd = "1"))

  # the consumer is encoded before the parameter exists, and resolves after
  obs <- mk_latent$model |>
    add_observation(record = list(model = "normal", mean = rho, sd = 1)) |>
    add_parameter(rho = list(prior = "gamma", shape = 2, scale = 1)) |>
    get_observations()
  obs[["args"]][[1]] |> expect_identical(c(mean = "rho|P", sd = "1"))
})

test_that("A status left without definitions names the status", {
  # definitions are left unresolved until a wall
  mk_latent$model |>
    add_latent(hidden_2 = list(func = min(datum_2), status = 2),
               hidden_2 = list(func = min(datum_2), status = 1)) |>
    get_latent() |>
    expect_nrow(2L)

  # every per-member definition is pinned to a status no individual carries
  mk_latent$model |>
    add_latent("hidden_2", func = min(datum_2), status = 2) |>
    evaluate_log_lik() |>
    expect_error("No.*latent.*definition.*status 1")

  # no status block to fall back on: the generic definition is the missing one
  data.frame(hh_id = "h1", status = 0, datum_2 = 2) |>
    new_housepi_model("test") |>
    add_latent("hidden_2", func = min(datum_2), status = c(2, 3)) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*generic.*latent.*definition")
})

test_that("Latent variable definitions are correct", {
  # basic constant definition
  augm <- mk_latent$model |>
    add_latent("hidden", value = 3) |>
    get_latent()
  augm$name |> expect_identical("hidden")
  augm$status |> expect_equal(NA_real_)
  augm$const |> expect_identical(TRUE)
  augm$discr |> expect_identical(FALSE)
  augm$value |> expect_equal(3)
  augm$proposal |> expect_equal(NA_real_)
  # a bare value is the same constant
  get_latent(add_latent(mk_latent$model, hidden = 3)) |> expect_identical(augm)

  # a model with no value is sampled by random walk
  augm <- mk_latent$model |>
    add_latent("hidden", model = "normal") |>
    get_latent()
  augm$const |> expect_identical(FALSE)
  augm$proposal |> expect_identical("rw")
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
    add_latent("hidden", value = datum, model = "uniform",
               min = layer_1, max = rho) |>
    get_latent("hidden")
  augm$status |> expect_equal(NA_real_)
  augm$value |> expect_identical("datum|D")
  augm$model |> expect_identical("uniform")
  augm$min |> expect_identical("layer_1|A")
  augm$max |> expect_identical("rho|P")
})

test_that("Positional arguments bind a variadic distribution only", {
  # categorical weights keep one empty key each
  augm <- mk_latent$model |>
    add_latent(hidden = list(model = "categorical", 0.2, 0.8)) |>
    get_latent()
  augm$model |> expect_identical("categorical")
  augm$args[[1]] |> expect_equal(c(0.2, 0.8) |> stats::setNames(c("", "")))

  # a fixed-arity distribution names its arguments
  mk_latent$model |>
    add_latent(hidden = list(model = "normal", 1, 2)) |>
    expect_error("[Uu]nnamed argument passed to distribution .normal.")
})

test_that("Latent independence proposal folds its arguments", {
  augm <- mk_latent$model |>
    add_latent("hidden", model = "uniform", min = 0, max = 5,
               proposal = "normal", proposal_mean = 0, proposal_sd = 1) |>
    get_latent()
  augm$proposal |> expect_identical("normal")
  augm$proposal_args[[1]] |> expect_equal(c(mean = 0, sd = 1))

  # random-walk scale is not accepted by an independence sampler
  mk_latent$model |>
    add_latent("hidden", model = "uniform", min = 0, max = 5,
               proposal = "normal", proposal_scale = 0.5) |>
    expect_error("random.?walk")
})

test_that("Independence proposal must cover the prior domain", {
  # a live minimum drawn in [0, 1] stays within the support of the proposal
  mk_latent$model |>
    add_parameter("rho", prior = "uniform") |>
    add_latent("hidden", model = "uniform", min = rho, max = 5,
               proposal = "gamma", proposal_shape = 2) |>
    helper_flat_sink(loc = hidden, column = "record") |>
    evaluate_log_lik() |>
    expect_all_finite()
  # a live minimum drawn in [-1, 1] reaches below the support of the proposal
  mk_latent$model |>
    add_parameter("rho", prior = "uniform", min = -1) |>
    add_latent("hidden", model = "uniform", min = rho, max = 5,
               proposal = "gamma", proposal_shape = 2) |>
    evaluate_log_lik() |>
    expect_error("support \\[0, Inf\\) that cannot cover")

  # a proposal bounded on both sides against an unbounded model
  mk_latent$model |>
    add_latent("hidden", model = "normal",
               proposal = "beta", proposal_shape1 = 2, proposal_shape2 = 2) |>
    evaluate_log_lik() |>
    expect_error("support \\[0, 1\\] that cannot cover")
  # a proposal minimum above the minimum of the model
  mk_latent$model |>
    add_latent("hidden", model = "beta", shape1 = 2, shape2 = 2, min = 0.2,
               proposal = "gamma", proposal_shape = 2, proposal_min = 0.5) |>
    evaluate_log_lik() |>
    expect_error("support \\[0.5, Inf\\) that cannot cover")
})

test_that("Latent function definition is correct", {
  hh_model <- mk_latent$model |>
    add_latent(layer_1 = list(model = "normal", mean = 0, sd = 1),
               layer_2 = list(model = "gamma", shape = 2, scale = 1),
               layer_3 = list(model = "lognormal", mean = 0, sd = 1))

  # deterministic definition: subtraction (as sum of minuses)
  augm <- hh_model |>
    add_latent("hidden", func = layer_1 - layer_2 - layer_3) |>
    get_latent("hidden")
  augm$value |> expect_equal(NA_real_)
  augm$func |> expect_identical("+")
  augm$func_args[[1]] |>
    expect_identical(c(x = "layer_1|A",
                       x2 = "minus|x|[layer_2|A]",
                       x3 = "minus|x|[layer_3|A]"))

  # alternative function call definitions
  by_args <- hh_model |>
    add_latent("hidden", func = "log", x = layer_1)
  by_expr <- hh_model |> add_latent(hidden = log(layer_1))
  get_latent(by_expr) |> expect_identical(get_latent(by_args))

  # an expression carries its own arguments: spare ones have no reading
  hh_model |>
    add_latent("hidden", func = log(layer_1), x = layer_1) |>
    expect_error("[Ee]xpression already carries its arguments")

  expect_warning(
    overwritten <- by_expr |> add_latent("hidden", func = log(layer_1)),
    "[Oo]verwriting"
  )
  get_latent(overwritten) |> expect_identical(get_latent(by_args))

  # multi-status definition: one row per status, identical otherwise
  multi <- hh_model |>
    add_latent("hidden", func = log(layer_1), status = c(1, 3)) |>
    get_latent("hidden")
  multi$status |> expect_identical(c("1", "3"))
  multi$status <- NULL
  multi[1, , drop = TRUE] |> expect_identical(multi[2, , drop = TRUE])
})

test_that("Symbolic definitions are equivalent", {
  hh_model <- mk_latent$model |>
    add_latent("layer_1", model = "uniform", min = 0, max = 2)

  # rlang double bang injects external values
  external_value <- 33
  hh_model_2 <- hh_model |>
    add_latent(
      hidden = !!external_value,
      hidden_2 = log(layer_1 + !!external_value),
      hidden_3 = c(
        model = "uniform",
        min = !!external_value, max = log(layer_1 + !!external_value)
      ),
      hidden_4 = list(
        model = "uniform",
        min = !!external_value, max = log(layer_1 + !!external_value)
      )
    )

  get_latent(hh_model_2, "hidden")$value |>
    expect_equal(external_value)
  get_latent(hh_model_2, "hidden_2")$func |>
    expect_identical("log")
  get_latent(hh_model_2, "hidden_2")$func_args[[1]] |>
    expect_equal(c(x = sprintf("+|x:x2|[layer_1|A]:[%s|N]", external_value)))

  # the two spellings of the same bounded definition agree
  by_vector <- get_latent(hh_model_2, "hidden_3")
  by_list <- get_latent(hh_model_2, "hidden_4")
  by_vector$name <- by_list$name <- NULL
  by_vector |> expect_equal(by_list, ignore_attr = TRUE)

  # a bare name, a bang and a formula encode the same read
  val <- 1
  augm <- hh_model |>
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
  unique(augm) |> expect_nrow(1L)
})

test_that("An external list of latent definitions can be imported", {
  # an external list is loaded through 'import', not as a bare definition
  external <- list(
    hidden = list(model = "gamma", shape = ~hidden_2, scale = ~hidden_3)
  )
  hh_model <- mk_latent$model |>
    add_latent(hidden_2 = list(model = "uniform", min = 1, max = 2),
               hidden_3 = list(model = "uniform", min = 1, max = 2))
  hh_model |>
    add_latent(external, append = FALSE) |>
    expect_error("non.?empty")

  augm <- hh_model |>
    add_latent(import = external) |>
    get_latent("hidden")
  augm$model |> expect_identical("gamma")
  augm$args[[1]] |>
    expect_identical(c(shape = "hidden_2|A", scale = "hidden_3|A"))

  # formulas import as the bare expressions they wrap
  mk_latent$model |>
    add_latent(import = list(layer_1 = ~ datum,
                             hidden = ~ sqrt(log(layer_1)))) |>
    get_latent() |>
    expect_equal(get_latent(add_latent(mk_latent$model, layer_1 = datum,
                                       hidden = sqrt(log(layer_1)))))
})

test_that("Latent template import loads", {
  hh_model <- mk_outbreak$model |>
    add_latent(hidden = end_date + 1, append = FALSE)

  # appending keeps the pre-existing definitions beside the template ones
  hh_model |>
    add_latent(import = "template") |>
    get_latent() |>
    getElement("name") |>
    expect_setequal(c("hidden", "infec_time"))

  # dropping resets the store first
  hh_model |>
    add_latent(import = "template", append = FALSE) |>
    get_latent() |>
    getElement("name") |>
    expect_setequal("infec_time")
})

test_that("Nested latent definitions raise errors", {
  hh_model <- data.frame(hh_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("test")

  hh_model |>
    add_latent(layer_2 = list(model = "uniform"),
               layer_1 = list(model = "uniform", status = 1),
               layer_1 = list(value = layer_2, status = c(2, 3, 4)),
               hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_all_finite()

  # simple nest
  hh_model |>
    add_latent(layer_1 = hidden, hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic latent definition")

  # triangular nest
  hh_model |>
    add_latent(hidden = layer_1, layer_1 = layer_2, layer_2 = hidden) |>
    evaluate_log_lik() |>
    expect_error("[Cc]yclic latent definition")

  # first generation of definitions: no layer_1 for status 3 and 4
  hh_model |>
    add_latent(layer_2 = list(model = "lognormal", status = 2),
               layer_1 = list(model = "normal", status = 1),
               layer_1 = list(value = layer_2, status = 2),
               hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")

  # second generation of definitions: no layer_2 for status 4
  hh_model |>
    add_latent(layer_2 = list(model = "lognormal", status = c(2, 3)),
               layer_1 = list(model = "normal", status = 1),
               layer_1 = list(value = layer_2, status = c(2, 3, 4)),
               hidden = layer_1,
               append = FALSE) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")
})

test_that("A discrete latent resolves 'uniform' to the discrete distribution", {
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"),
                         status = c(1, 1)) |>
    new_housepi_model("test")

  hh_model |>
    add_latent("hidden", model = "uniform", min = 0, max = 3, discr = TRUE) |>
    evaluate_log_lik() |>
    expect_all_finite()

  # the flat lattice law owns every integer value in the window, so a bound that
  # is not one has no reading, and the error names the latent that carries it
  sms <- hh_model |>
    add_latent("hidden", model = "uniform", min = 0.5, max = 3, discr = TRUE) |>
    evaluate_log_lik() |>
    tryCatch(error = conditionMessage)
  sms |> expect_match("min.*uniform.*discrete")
  sms |> expect_match("hidden")
})

test_that("A location off the integers is refused when discrete", {
  hh_model <- data.frame(hh_id = 1, id = c("m1", "m2"), status = c(1, 0)) |>
    new_housepi_model("test")

  # a sampled location would shift the lattice under the chain at every step:
  # the fix is on the declaration it reads, named in place of `loc`
  hh_model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 1) |>
    add_latent("hidden", model = "poisson", lambda = 2, loc = rho) |>
    evaluate_log_lik() |>
    expect_error("[Cc]hange parameter.*rho.*discrete")
})

test_that("A latent bound to a distribution returns its density, not a draw", {
  # a distribution called as a func gives back a mass: reading its support here
  # would type the column integer and truncate every mass to zero
  simulated <- data.frame(hh_id = rep(seq(3), each = 2), id = seq(6),
                          status = 1) |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "uniform", max = 5, discr = TRUE),
               hidden = list(func = "poisson", x = layer_1, lambda = 3)) |>
    run_simulation()

  simulated$layer_1 |> expect_type("integer")
  simulated$hidden |> expect_type("double")
  simulated$hidden |> expect_equal(dpois(simulated$layer_1, 3))
})

test_that("Column reads inherit the discretization of their column", {
  # sources resolve before `discr` does, so a read reaches the column and
  # inherits its storage type
  augm <- data.frame(hh_id = rep(c("h1", "h2"), each = 2),
                     datum = as.Date("2005-02-28") + rep(0:1, 2),
                     datum_2 = rep(c(3.5, 5.5), 2),
                     status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "normal", value = datum),
               hidden_2 = list(model = "normal", value = datum_2)) |>
    get_latent()

  augm$value |> expect_identical(c("datum|D", "datum_2|D"))
  # the date column is discrete, the double one is not
  augm$discr |> expect_identical(c(TRUE, FALSE))
})

###################################
###  OBSERVATION : DEFINITIONS  ###
###################################

test_that("Observation definitions are encoded", {
  obs <- mk_latent$model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(
      record = list(model = "normal", mean = hidden, sd = 1)
    ) |>
    get_observations()

  # marked as observation, data-bound and constant (never sampled)
  obs[["level"]] |> expect_identical("member")
  obs[["const"]] |> expect_true()
  obs[["value"]] |> expect_identical("record|D")
  obs[["proposal"]] |> expect_all_na()
  # distribution carries the observation arguments
  obs[["model"]] |> expect_identical("normal")
  obs[["args"]][[1]] |> expect_identical(c(mean = "hidden|A", sd = "1"))

  # a deterministic observation carries a function and no distribution
  obs <- mk_latent$model |>
    add_latent(layer_1 = list(model = "gamma", shape = 2, scale = 1),
               hidden = -layer_1) |>
    add_observation(record = floor(hidden + layer_1)) |>
    get_observations()

  obs[["level"]] |> expect_identical("member")
  obs[["model"]] |> expect_all_na()
  obs[["value"]] |> expect_all_na()
  obs[["func"]] |> expect_match("floor")
})

test_that("Observation node rejects conflicting options", {
  # observation nodes are data-bound, never sampled: the options are dropped
  ignored <- function(...) {
    expect_warning(hh_model <- add_observation(mk_latent$model, ...),
                   "[Ii]gnored.*const.*proposal")
    get_observations(hh_model)
  }
  obs <- ignored("record", model = "normal", const = FALSE,
                 proposal = "forward")
  obs[["const"]] |> expect_true()
  obs[["proposal"]] |> expect_all_na()

  # a deterministic observation drops them alike
  obs <- ignored(record_2 = list(func = min(datum_2), status = 1,
                                 const = TRUE, proposal = "rw"))
  obs[["func"]] |> expect_identical("min")
  obs[["const"]] |> expect_false()
  obs[["proposal"]] |> expect_all_na()

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
  hh_model <- mk_latent$model
  err <- "cannot be status specific"

  # guarded at the bulk default level and per definition
  add_observation(hh_model,
    record_2 = list(func = min(datum_2), status = 1), level = "house"
  ) |>
    expect_error(err)
  add_observation(hh_model,
    record_2 = list(func = min(datum_2), status = 1, level = "house")
  ) |>
    expect_error(err)
  add_observation(hh_model, record_2 = list(func = min(datum_2)),
                  status = 1, level = "house") |>
    expect_error(err)

  # a member-level declaration keeps its status, and the no-effect reduction
  # that would promote it to household level is not taken when one is set
  member <- hh_model |>
    add_observation(record_2 = list(func = min(datum_2), status = 1)) |>
    get_observations()
  member[["level"]] |> expect_identical("member")
  as.integer(member[["status"]]) |> expect_identical(1L)

  # a bulk 'status =' reaches every entry, single and bulk modes alike
  bulk <- hh_model |>
    add_observation(record_2 = list(func = min(datum_2)), status = 1) |>
    get_observations()
  as.integer(bulk[["status"]]) |> expect_identical(1L)
  imported <- hh_model |>
    add_observation(
      import = list(record_2 = list(func = quote(min(datum_2)))),
      status = 1
    ) |>
    get_observations()
  as.integer(imported[["status"]]) |> expect_identical(1L)
})

test_that("A name is declared under a single kind", {
  hh_model <- mk_latent$model
  err <- "already declared as"

  # latent and observations share one store: every pair collides, either way
  hh_model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(hidden = datum_2) |>
    expect_error(err)
  hh_model |>
    add_observation(hidden = datum_2) |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    expect_error(err)
  hh_model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(hidden = min(datum_2), level = "house") |>
    expect_error(err)
  hh_model |>
    add_observation(hidden = min(datum_2), level = "house") |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    expect_error(err)
  hh_model |>
    add_observation(hidden = datum_2) |>
    add_observation(hidden = min(datum_2), level = "house") |>
    expect_error(err)

  # a status-specific declaration is no loophole: the store is keyed by name
  hh_model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1, status = 1)) |>
    add_observation(hidden = min(datum_2), level = "house") |>
    expect_error(err)

  # a leading dot is reserved to the definitions housepi generates itself
  hh_model |>
    add_latent(.status = list(value = 1, const = TRUE)) |>
    expect_error("cannot start with a dot")

  # siblings of one kind stay legal: statuses of the same name, and overwrites
  hh_model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1, status = 1)) |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1, status = 2)) |>
    get_latent() |>
    expect_nrow(2L)
  hh_model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 2)) |>
    expect_warning("Overwriting")

  # an observation standing in for a like-named data column is the surrogate
  # design, not a clash: only declarations collide with each other
  hh_model |>
    add_observation(datum_2 = record) |>
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
      record_2 = list(func = min(datum_2), level = "house"),
      record_3 = list(func = 21, level = "member")
    ) |>
    get_observations()
  level_of(mixed) |> expect_identical(expected)
  ("args" %in% names(mixed)) |> expect_false()

  # the argument is the default: entries not setting their own inherit it
  mk_latent$model |>
    add_observation(
      record_2 = list(func = min(datum_2)),
      record_3 = list(func = 21, level = "member"),
      level = "house"
    ) |>
    get_observations() |>
    level_of() |>
    expect_identical(expected)

  # imported definitions honour it too
  mk_latent$model |>
    add_observation(import = list(
      record_2 = list(func = quote(min(datum_2)), level = "house"),
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
  hh_model <- mk_latent$model |>
    add_latent(hidden = list(model = "gamma", shape = 2, scale = 1)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1))

  # non-append household observation keeps latent definitions
  hh_model <- add_observation(
    hh_model, record_2 = min(datum_2), level = "house", append = FALSE
  )
  get_latent(hh_model)[["name"]] |> expect_identical("hidden")
  obs <- get_observations(hh_model)

  # observations are a single family: the member row went with the reset
  obs[["name"]] |> expect_setequal("record_2")
  obs[["level"]] |> expect_identical("house")
  obs[["func"]] |> expect_identical("min")
  obs[["func_args"]][[1]] |> expect_identical(c(x = "datum_2|D"))

  # appending keeps both levels side by side
  hh_model <- add_observation(
    hh_model, record_3 = max(datum_2), level = "house"
  )
  get_observations(hh_model)[["name"]] |>
    expect_setequal(c("record_2", "record_3"))
  get_latent(hh_model)[["name"]] |> expect_identical("hidden")

  # symmetrically, a non-append latent load keeps observations
  hh_model <- add_latent(hh_model, append = FALSE,
                         layer_1 = list(model = "gamma", shape = 2, scale = 1))
  get_latent(hh_model)[["name"]] |> expect_identical("layer_1")
  get_observations(hh_model)[["name"]] |>
    expect_setequal(c("record_2", "record_3"))
})

test_that("A degenerate reduction at member level is cast to house level", {
  # member-level reduction by mistake
  expect_warning(
    promoted <- mk_latent$model |>
      add_observation(record_2 = min(datum_2)) |>
      get_observations(),
    regexp = "no-effect.*reduction.*min"
  )
  promoted[["level"]] |> expect_identical("house")
  # explicit house-level reduction
  mk_latent$model |>
    add_observation(record_2 = min(datum_2), level = "house") |>
    get_observations() |>
    expect_identical(promoted)

  # redundant reductions raise an error and are not cast
  mk_latent$model |>
    add_observation(record_2 = min(datum_2), level = "house") |>
    add_observation(record_3 = min(record_2)) |>
    expect_error("redundant.*record_3")

  # genuine per-member expressions pass normally
  mk_latent$model |>
    add_observation(record_2 = min(datum_2, 0),
                    record_3 = datum_2,
                    record_4 = log(datum_2)) |>
    get_observations() |>
    getElement("level") |>
    expect_identical(rep("member", 3))
})

test_that("A missing observation leaves the likelihood", {
  # an observation value is a datum: where the cell is empty there is nothing
  # to condition on, so the term drops out instead of entering as NA
  mk <- function(x) {
    data.frame(id = c("m1", "m2"), hh_id = "h1", status = 1, record = x) |>
      new_housepi_model("test") |>
      add_latent("hidden", value = 0.5) |>
      add_observation("record", model = "uniform", loc = hidden,
                      min = -0.4, max = 0.4)
  }

  # both members observed: two terms of a uniform over a window of 0.8
  evaluate_log_lik(mk(c(0.5, 0.5))) |> expect_equal(2 * log(1 / 0.8))

  # one datum absent: one term left, and the run says which column and how many
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA))),
                 "[Mm]issing for 1 of 2 infected individuals")
  log_lik |> expect_equal(log(1 / 0.8))

  # the domain still binds where the datum is there
  evaluate_log_lik(mk(c(0.5, 9))) |> expect_equal(-Inf)

  # no datum anywhere: the model declares a density nothing feeds
  mk(c(NA_real_, NA_real_)) |>
    evaluate_log_lik() |>
    expect_error("missing for every infected individual")

  # a categorical column reaches the engine as integer codes: its empty cell
  # must code to no level, or the density would read the code of a non-value
  mk_categorical <- function(x) {
    data.frame(id = c("m1", "m2"), hh_id = "h1", status = 1, record = x) |>
      new_housepi_model("test") |>
      add_observation("record", model = "uniform", min = 0, max = 3)
  }
  expect_warning(
    log_lik <- evaluate_log_lik(mk_categorical(c("a", NA))),
    "[Mm]issing for 1 of 2 infected individuals"
  )
  log_lik |> expect_equal(log(1 / 3))
  mk_categorical(c(NA, NA)) |>
    evaluate_log_lik() |>
    expect_error("missing for every infected individual")

  # an observation is judged over the members carrying it: the status-1 member
  # never declares a status-2 observation, so its empty cell is not missing
  data.frame(id = c("m1", "m2"), hh_id = "h1",
             status = c(1, 2), record = c(NA, 0.5)) |>
    new_housepi_model("test") |>
    add_latent("hidden", value = 0.5) |>
    add_observation(record = list(model = "uniform", loc = hidden,
                                  min = -0.4, max = 0.4, status = 2)) |>
    evaluate_log_lik() |>
    expect_equal(log(1 / 0.8))
})

test_that("A latent goes with the observation that reads it", {
  # a latent nothing else reads is an unused parameter once its observation has
  # left the likelihood: marginalizing it integrates its model to one, so its
  # density goes with the datum. A flat model keeps the total free of the draw.
  mk <- function(x, x_2 = NULL, hop = FALSE) {
    hhdata <- data.frame(id = c("m1", "m2"), hh_id = "h1",
                         status = 1, record = x)
    if (!is.null(x_2)) hhdata$record_2 <- x_2
    hh_model <- hhdata |>
      new_housepi_model("test") |>
      add_latent(hidden = list(model = "uniform", min = 0, max = 2))
    hh_model <- if (hop) {
      hh_model |>
        add_latent(hidden_2 = hidden) |>
        add_observation("record", model = "uniform", loc = hidden_2,
                        min = -5, max = 5)
    } else {
      add_observation(hh_model, "record", model = "uniform", loc = hidden,
                      min = -5, max = 5)
    }
    if (is.null(x_2)) hh_model
    else add_observation(hh_model, "record_2", model = "uniform",
                         loc = hidden, min = -5, max = 5)
  }

  # both observed: two flat latent densities and two observation terms
  evaluate_log_lik(mk(c(0.5, 0.5))) |>
    expect_equal(2 * log(1 / 2) + 2 * log(1 / 10))

  # m2 carries no datum: its observation and the latent feeding it both drop
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA))),
                 "[Mm]issing for 1 of 2")
  log_lik |> expect_equal(log(1 / 2) + log(1 / 10))

  # the same one hop away, through a deterministic readout
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA), hop = TRUE)),
                 "[Mm]issing for 1 of 2")
  log_lik |> expect_equal(log(1 / 2) + log(1 / 10))

  # a second observation still reads it, so m2 keeps its draw
  expect_warning(log_lik <- evaluate_log_lik(mk(c(0.5, NA), x_2 = c(0.5, 0.5))),
                 "[Mm]issing for 1 of 2")
  log_lik |> expect_equal(2 * log(1 / 2) + 3 * log(1 / 10))

  # simulations generate every observation, so they prune nothing
  simulated <- run_simulation(mk(c(0.5, 0.5)))
  simulated$record |> expect_no_na()
})

test_that("Reading an observation in models and loops are rejected", {
  # a latent model must not read an observation (illegal in both
  # modes; caught at eager build, so the construction itself errors)
  data.frame(hh_id = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_observation("record", model = "uniform", min = 0, max = 1) |>
    add_latent("hidden", model = "normal", mean = record) |>
    expect_error("[Mm]odel.*reads.*observation")

  # a deterministic read that closes a loop (the observation's support
  # depends on the reader) is circular
  data.frame(hh_id = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_latent("hidden", value = record) |>
    add_observation("record", model = "uniform", min = hidden, max = 10) |>
    expect_error("[Cc]yclic")
})

test_that("A simulation override (value + model) rejects observation reads", {
  hh_model <- data.frame(hh_id = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_observation("record", model = "uniform", min = 0, max = 1)

  # the override is applied before the observation is generated
  hh_model |>
    add_latent("hidden", model = "normal", value = record) |>
    expect_error("[Ss]imulation override.*observation")

  # a deterministic readout does not make the observation available any earlier,
  # so the same read one hop away is rejected too
  hh_model |>
    add_latent(layer_1 = record,
               hidden = list(model = "normal", value = layer_1)) |>
    expect_error("[Ss]imulation override.*record.*through.*layer_1")
})

test_that("A surrogate column reading latent is simulation only", {
  hh_model <- data.frame(id = seq_len(10), hh_id = "h1",
                         status = c(rep(1, 9), 0), record = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record = hidden + 2)

  # the data view has no derived value to report: it keeps the imported datum
  get_data_wide(hh_model)[["record"]] |> expect_all_equal(0)

  # observation depends on sampled latent: not a valid surrogate definition
  evaluate_log_lik(hh_model) |>
    expect_error("must resolve from data alone")

  # simulation produces it for infected members, following the latent draw
  simulated <- run_simulation(hh_model)
  infected <- simulated$status > 0
  sum(names(simulated) == "record") |> expect_identical(1L)
  simulated$record[infected] |> expect_equal(simulated$hidden[infected] + 2)
  simulated$record[!infected] |> expect_all_na()
})

test_that("A surrogate observation overwrites the corresponding data column", {
  hh_model <- data.frame(hh_id = "h1", status = c(1, 0, 0),
                         datum_2 = c(2, 6, NA), record_2 = rep(99, 3)) |>
    new_housepi_model("test") |>
    add_observation(
      record_2 = list(value = min(datum_2), level = "house"),
      record_3 = datum_2 + 1
    ) |>
    add_latent(hidden = list(value = record_3 * 2))

  # the data view reports the derived value before evaluation
  unique(get_data_wide(hh_model)[["record_2"]]) |> expect_identical(2)
  ("record_3" %in% colnames(get_data_wide(hh_model))) |> expect_false()

  # and after evaluation
  evaluate_log_lik(hh_model) |> expect_all_finite()
  unique(get_data_wide(hh_model)[["record_2"]]) |> expect_identical(2)
  ("record_3" %in% colnames(get_data_wide(hh_model))) |> expect_false()

  # simulated observations are for all members, latent just for infected
  simulated <- run_simulation(hh_model)
  simulated[["record_3"]] |> expect_equal(c(3, 7, NA))
  simulated[["hidden"]] |> expect_equal(c(6, NA, NA))

  # missing datum is reported and attributed to the observation reading it
  data.frame(hh_id = "h1", status = 1, record = 5) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record_2 = list(model = "normal", mean = hidden, sd = 1)) |>
    add_observation(record = record_2 + 2) |>
    evaluate_log_lik() |>
    expect_error("missing data column.*record_2.*[Rr]equired by observation")
})

test_that("Surrogate columns of both levels build in dependency order", {
  hhdata <- data.frame(hh_id = "h1", status = c(1, 0, 0),
                       datum = c(2, 9, 4), record = 0, record_2 = 0)
  hh_model <- hhdata |>
    new_housepi_model("test") |>
    add_observation(record = datum + 3,
                    record_2 = list(value = max(record), level = "house"))

  # every member carries a value regardless of status
  hhdata_2 <- get_data_wide(hh_model)
  hhdata_2[["record"]] |> expect_identical(c(5, 12, 7))
  unique(hhdata_2[["record_2"]]) |> expect_identical(12)

  # pure surrogate observations with no raw data columns are not exported
  hh_model |>
    add_data_wide(hhdata[c("hh_id", "status", "datum")]) |>
    get_data_wide() |>
    colnames() |>
    intersect(c("record", "record_2")) |>
    expect_length(0L)
})

test_that("Reading only house observations casts member level to house", {
  hhdata <- data.frame(hh_id = "h1", status = c(1, 0, 0), datum = c(2, 9, 4),
                       record = 0, record_2 = 0, record_3 = 0,
                       record_4 = 0, record_5 = 0, record_6 = 0)
  hh_model <- hhdata |>
    new_housepi_model("test") |>
    add_observation(record = list(value = max(datum), level = "house"),
                    record_2 = record + 1,
                    record_3 = record_2 * 2,
                    record_4 = datum + 1, # <- member level
                    record_5 = record - 1,
                    record_6 = 18 / record)

  # a cast row is a household value to the rows reading it
  obs <- hh_model$observations
  names <- c("record_2", "record_3", "record_4", "record_5", "record_6")
  obs$level[match(names, obs$name)] |>
    expect_identical(c("house", "house", "member", "house", "house"))

  hhdata_2 <- get_data_wide(hh_model)
  hhdata_2[["record_2"]] |> expect_equal(c(10, 10, 10))
  hhdata_2[["record_3"]] |> expect_equal(c(20, 20, 20))
  hhdata_2[["record_4"]] |> expect_equal(c(3, 10, 5))
  hhdata_2[["record_5"]] |> expect_equal(c(8, 8, 8))
  hhdata_2[["record_6"]] |> expect_equal(c(2, 2, 2))
})

test_that("A household reduction folds its argument over the members", {
  # datum differs across members, the first one holds 2
  hhdata <- data.frame(hh_id = "h1", status = c(1, 0, 0),
                       datum = c(2, 9, 4), record = 0)
  reduce <- function(defin, data = hhdata) {
    data |>
      new_housepi_model("test") |>
      add_observation("record", value = !!defin, level = "house") |>
      get_data_wide() |>
      getElement("record") |>
      unique()
  }

  # a nested reduction reads every member
  reduce(quote(max(datum) + 1)) |> expect_equal(10)
  reduce(quote(log(max(datum)))) |> expect_equal(log(9))
  reduce(quote(max(datum) - min(datum))) |> expect_equal(7)
  # the argument is computed for each member, then folded
  reduce(quote(max(datum + 1))) |> expect_equal(10)
  reduce(quote(sum(datum + 1))) |> expect_equal(18)
  reduce(quote(max(log(datum)))) |> expect_equal(log(9))

  # a member with no datum drops out of the fold
  reduce(quote(max(datum + 1)), transform(hhdata, datum = c(2, NA, 4))) |>
    expect_equal(5)
  # with no datum at all there is nothing to fold
  reduce(quote(max(datum) + 1), transform(hhdata, datum = NA_real_)) |>
    is.na() |>
    expect_true()

  # a per-member column outside any reduction has no household value
  reduce(quote(datum + 1)) |>
    expect_error("[Pp]er-member data .datum. is read outside a reduction")
  reduce(quote(log(datum))) |>
    expect_error("[Pp]er-member data .datum. is read outside a reduction")
  reduce(quote(datum)) |>
    expect_error("[Pp]er-member data .datum. is read outside a reduction")
})

test_that("A household observation with no datum yet reads as NA", {
  # 'end_date' reduces over an observation the run has still to draw: there
  # is nothing to reduce, so the data view reports NA rather than reading the
  # empty reduction as a datum failing the discreteness 'case_ascertained' asks
  hhdata <- data.frame(id = c("m1", "m2", "m3", "m4"),
                       hh_id = c("h1", "h1", "h2", "h2"),
                       status = NA) |>
    new_housepi_model("case-ascertained") |>
    add_observation(
      onset_date = list(model = "normal", discr = TRUE, mean = 3, sd = 1),
      end_date = list(value = max(onset_date), level = "house")
    ) |>
    get_data_wide()

  hhdata$end_date |> expect_all_na()
})

test_that("A surrogate column proves the requirement its model asks", {
  # a surrogate filling 'end_date' answers the discreteness 'case_ascertained'
  # asks of the column, and answers for what it can return, not for the value it
  # happens to settle on. A household-level definition reduces over its members,
  # so its root is the reduction; the member level takes the algebra directly
  mk <- function(defin, level = "house") {
    data.frame(id = c("m1", "m2", "m3", "m4"),
               hh_id = c("h1", "h1", "h2", "h2"),
               status = c(1, 0, 1, 0), datum = 0.5, datum_2 = 16) |>
      new_housepi_model("case-ascertained") |>
      add_observation("end_date", value = !!defin, level = level) |>
      get_data_wide()
  }

  # integer-valued cells fold to the integer lattice through the algebra
  mk(quote(min(datum_2)))$end_date |> expect_equal(rep(16, 4))
  mk(quote(sum(datum_2)))$end_date |> expect_equal(rep(32, 4))
  mk(quote(datum_2 * 2 - 1), level = "member")$end_date |>
    expect_equal(rep(31, 4))

  # two halves sum to an integer, and are refused all the same: the column
  # answers for what the definition can return, not for the value it settles on
  mk(quote(sum(datum))) |>
    expect_error("not discrete-valued")
  mk(quote(datum_2 / 4), level = "member") |>
    expect_error("not discrete-valued")
  mk(quote(sqrt(datum_2)), level = "member") |>
    expect_error("not discrete-valued")

  # rounding restores the proof
  mk(quote(round(datum_2 / 4)), level = "member")$end_date |>
    expect_equal(rep(4, 4))
  mk(quote(ceil(sqrt(datum_2))), level = "member")$end_date |>
    expect_equal(rep(4, 4))

  # the refusal names the definition and the place it was built for
  sms <- tryCatch(mk(quote(sum(datum))), error = conditionMessage)
  sms |> expect_match("[Hh]ousehold observation")
  sms |> expect_match("h1")
  sms <- tryCatch(mk(quote(datum_2 / 4), level = "member"),
                  error = conditionMessage)
  sms |> expect_match("individual")
  sms |> expect_match("m1")
})

test_that("Data-bound household observation rejects latent dependence", {
  # a household observation standing in for data resolves from data alone
  mk_outbreak$model |>
    add_observation(start_date = min(test_date) + alpha,
                    level = "house") |>
    evaluate_log_lik() |>
    expect_error("must resolve from data alone")
})

############################
## STATUS : DEFINITIONS   ##
############################

test_that("Escaped label colliding with an infected status is rejected", {
  mk_status$data |>
    new_housepi_model("case-ascertained", escaped = "symptomatic") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
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

#############################
###  LATENT : END-TO-END  ###
#############################

test_that("Independence proposal inherits the model location shift", {
  results <- data.frame(id = c("m1", "m2"), hh_id = "h1",
                        status = 1, datum = c(4, 8), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", loc = datum,
               proposal = "gamma", proposal_shape = 2) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  draws <- as.data.frame(get_latent_draws(results))
  (draws$m1_hidden >= 4 & draws$m1_hidden <= 5) |> expect_all_true()
  (draws$m2_hidden >= 8 & draws$m2_hidden <= 9) |> expect_all_true()
})

test_that("Truncated latent models are normalized in log-likelihood traces", {
  # bounds read data, constants and a sampled parameter
  hhdata <- data.frame(id = c("m1", "m2"), hh_id = c("h1", "h2"),
                       status = 1, datum = c(3, 6), obs = 0)
  results <- hhdata |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform", min = 1, max = 3)) |>
    add_latent(
      layer_1 = list(model = "gamma", shape = 2, scale = 1.5, max = datum),
      hidden = list(model = "gamma", shape = rho, min = 0.5, max = 8),
      hidden_2 = list(model = "gamma", shape = 1.7, max = 5)
    ) |>
    helper_flat_sink(loc = layer_1 + hidden + hidden_2) |>
    run_inference(steps = 50L, adapt = 0L, thin = 1L, chains = 1L,
                  id_draws = c("m1", "m2")) |>
    suppressWarnings()

  draws <- latent_frame(results)
  rho <- draws_to_frame(get_parameter_draws(results))$rho
  log_lik <- draws_to_frame(get_log_lik_draws(results))
  log_dtrunc <- function(x, shape, scale, min, max) {
    dgamma(x, shape, scale = scale, log = TRUE) -
      log(pgamma(max, shape, scale = scale) - pgamma(min, shape, scale = scale))
  }
  # one member per household, the flat sink adds log(1 / 2e6)
  for (k in seq(2)) {
    col <- \(name) draws[[sprintf("m%d_%s", k, name)]]
    expected <- log_dtrunc(col("layer_1"), 2, 1.5, 0, hhdata$datum[k]) +
      log_dtrunc(col("hidden"), rho, 1, 0.5, 8) +
      log_dtrunc(col("hidden_2"), 1.7, 1, 0, 5) - log(2e6)
    log_lik[[sprintf("h%d", k)]] |> expect_equal(expected, tolerance = 1e-8)
  }
})

test_that("Latent draws stay inside their live support", {
  hhdata <- data.frame(id = seq(2), hh_id = "h1", status = 1, obs = 0)
  run_short <- function(hh_model) {
    hh_model |>
      run_inference(steps = 100L, adapt = 0L, thin = 1L, chains = 1L) |>
      suppressWarnings()
  }

  # other latent variable
  draws <- hhdata |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "uniform"),
               hidden = list(model = "uniform", loc = layer_1)) |>
    helper_flat_sink(loc = hidden) |>
    run_short() |>
    latent_frame()
  gap <- draws[["1_hidden"]] - draws[["1_layer_1"]]
  (gap >= 0 & gap <= 1) |> expect_all_true()

  # transitive deterministic latent variable
  draws <- hhdata |>
    new_housepi_model("test") |>
    add_latent(layer_2 = list(model = "uniform"),
               layer_1 = layer_2,
               hidden = list(model = "uniform", loc = layer_1)) |>
    helper_flat_sink(loc = hidden) |>
    run_short() |>
    latent_frame()
  gap <- draws[["1_hidden"]] - draws[["1_layer_2"]]
  (gap >= 0 & gap <= 1) |> expect_all_true()

  # parameter
  results <- hhdata |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform") |>
    add_latent("hidden", model = "uniform", min = rho, max = 5) |>
    helper_flat_sink(loc = hidden) |>
    run_short()
  hidden <- latent_frame(results)[["1_hidden"]]
  rho <- as.data.frame(get_parameter_draws(results))[["rho"]]
  (hidden >= rho) |> expect_all_true()
})

test_that("Draws requested where no latent is sampled are ignored", {
  # a member carrying no infected status draws no latent
  data.frame(id = c("m1", "m2"), hh_id = "h1", status = c(1, 0), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    helper_flat_sink(loc = hidden) |>
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
                          hh_id = rep(seq(10), each = 10),
                          status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform",
                             min = layer_1, max = layer_1 + 1),
               layer_1 = list(model = "gamma", shape = 2, scale = 1)) |>
    run_simulation()
  (simulated$layer_1 <= simulated$hidden &
     simulated$hidden <= simulated$layer_1 + 1) |>
    expect_all_true()
})

test_that("A value paired with model overrides the draw in simulations only", {
  hh_model <- data.frame(id = seq_len(100), hh_id = 1,
                         status = 1, datum = 7, obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", max = 100, value = datum) |>
    helper_flat_sink(loc = hidden)

  # simulation: pins the sampled variable to user value
  run_simulation(hh_model)$hidden |> expect_all_equal(7)

  # inference: the value does not pin the sampler and draws span the prior
  draws <- hh_model |>
    run_inference(steps = 40L, chains = 1L) |>
    get_latent_draws() |>
    suppressWarnings()
  all(draws[["1_hidden"]] == 7) |> expect_false()

  # a pin composed on another latent binds live: the relation holds member by
  # member, and follows the draw rather than collapsing to a single value
  simulated <- data.frame(id = seq_len(100), hh_id = 1,
                          status = 1, datum = 7) |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "gamma", shape = 2, scale = 1),
               hidden = list(model = "uniform", max = 100,
                             value = datum - layer_1)) |>
    run_simulation()
  simulated$hidden |> expect_equal(7 - simulated$layer_1)
  length(unique(simulated$layer_1)) |> expect_gt(1)
})

test_that("A latent without a model is one transmission generates", {
  # incubation periods are no transmission output: they need a model
  mk_status$model |>
    add_latent(incub_period = list(min = 0, max = 14)) |>
    expect_error("no model.*transmission model does not generate")
  # infection times are: transmission gives their density
  mk_status$model |>
    add_latent(infec_time = list(model = "uniform", min = -14, max = end_date),
               append = FALSE) |>
    expect_error("generated by the transmission model")
})

test_that("Bounds of a flat latent read neither draws nor free parameters", {
  # a free parameter moves the bound with each draw
  mk_outbreak$model |>
    add_parameter(window = list(prior = "uniform", min = 1, max = 30)) |>
    add_latent(infec_time = list(min = -window, max = end_date),
               append = FALSE) |>
    expect_error("[Bb]ound .min. of latent .infec_time. reads parameter")
  # so does another latent
  mk_outbreak$model |>
    add_latent(incub_period = list(model = "gamma", shape = 2, scale = 1),
               infec_time = list(min = -14, max = incub_period),
               append = FALSE) |>
    expect_error("[Bb]ound .max. of latent .infec_time. reads latent")

  # a constant parameter does not
  mk_outbreak$model |>
    add_parameter(window = 14) |>
    add_latent(infec_time = list(min = -window, max = end_date),
               append = FALSE) |>
    evaluate_log_lik() |>
    expect_all_finite()
  # nor does an observation, read as its datum even when the latent generates it
  results <- mk_onset$model |>
    add_latent(incub_period = list(model = "gamma", shape = 2, scale = 1),
               infec_time = list(min = -14, max = onset_date + 1),
               onset_time = infec_time + incub_period, append = FALSE) |>
    add_observation("onset_date", model = "normal", discr = TRUE,
                    mean = onset_time, sd = 1, min = floor(onset_time)) |>
    run_inference(steps = 10L, thin = 1L, adapt = 0L, chains = 1L,
                  progress = FALSE, warn = FALSE)
  get_log_lik_draws(results)$h1 |> expect_all_finite()
})

test_that("A flat latent is set by transmission, else by its value", {
  # transmission sets every infection but the first, which the value pins
  # a reader follows the time transmission sets, not the value pin
  simulated <- mk_outbreak$model |>
    add_latent(infec_time = list(value = 0, min = -14, max = end_date,
                                 proposal = "rw"),
               append = FALSE) |>
    add_observation(record = infec_time + 1) |>
    run_simulation()
  tapply(simulated$infec_time, simulated$hh_id, min, na.rm = TRUE) |>
    as.vector() |>
    expect_all_equal(0)
  any(simulated$infec_time > 0) |> expect_true()
  infected <- !is.na(simulated$infec_time)
  simulated$record[infected] |> expect_equal(simulated$infec_time[infected] + 1)

  # with no model and no value, the first infection has nothing to draw from
  mk_outbreak$model |>
    add_latent(infec_time = list(min = -14, max = end_date), append = FALSE) |>
    run_simulation() |>
    expect_error("[Nn]o .value. to simulate latent .infec_time.")
})

test_that("A simulation override orders its producers per run", {
  hh_model <- data.frame(id = sprintf("m%s", 1:4), hh_id = "h1",
                         status = 1, obs = 0) |>
    new_housepi_model("test")

  # simulation: hidden reads layer_1 through the value
  # inference: layer_1 reads hidden through the model
  # neither run closes a loop
  hh_model_2 <- hh_model |>
    add_latent(hidden = list(model = "normal", value = layer_1),
               layer_1 = list(model = "normal", mean = hidden, sd = 1,
                              value = 5)) |>
    helper_flat_sink(loc = hidden + layer_1)
  results <- run_simulation(hh_model_2)
  results$hidden |> expect_all_equal(5)
  results$layer_1 |> expect_all_equal(5)
  draws <- hh_model_2 |>
    run_inference(steps = 40L, chains = 1L) |>
    get_latent_draws() |>
    suppressWarnings()
  all(draws[["m1_layer_1"]] == 5) |> expect_false()

  # a support bound read by the model is an inference edge too
  results <- hh_model |>
    add_latent(hidden = list(model = "normal", value = layer_1),
               layer_1 = list(model = "uniform", min = hidden - 1,
                              max = hidden + 1, value = 2)) |>
    helper_flat_sink(loc = hidden + layer_1) |>
    run_simulation()
  results$hidden |> expect_all_equal(2)

  # simulation draws layer_2 before hidden reads it through the pin layer_1
  # inference orders hidden before layer_2
  results <- hh_model |>
    add_latent(layer_1 = list(model = "normal", value = layer_2),
               hidden = list(model = "normal", mean = layer_1, sd = 1e-3),
               layer_3 = list(model = "normal", mean = hidden, sd = 1,
                              value = 50),
               layer_2 = list(model = "normal", mean = layer_3, sd = 1)) |>
    helper_flat_sink(loc = hidden) |>
    run_simulation()
  (abs(results$hidden - results$layer_2) < 0.01) |> expect_all_true()

  # a loop closed in simulation only
  hh_model |>
    add_latent(hidden = list(model = "normal", value = hidden_2),
               hidden_2 = list(model = "normal", mean = hidden, sd = 1)) |>
    expect_error("[Cc]yclic.*hidden.*[Ss]imulation reads the value")

  # a loop closed in inference only
  hh_model |>
    add_latent(hidden = list(model = "normal", mean = hidden_2, value = 3),
               hidden_2 = list(model = "normal", mean = hidden, sd = 1)) |>
    expect_error("[Cc]yclic.*hidden")
})

test_that("Latent traces export only relevant variables", {
  # a constant carries no draw, an implicit function cannot be called on a
  # value, and what no sink reads is not sampled at all
  hh_model <- data.frame(id = c("m1", "m2"), hh_id = "h1",
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
      rho_3 = rho,
      rho_sum = rho + rho_2,
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
  observed <- add_observation(hh_model, "record", model = "normal",
                              mean = expl + rho_sum, sd = 1)

  traced <- function(hh_model) {
    hh_model |>
      run_inference(steps = 100L, burn = 0L, thin = 1L, adapt = 0L,
                    chains = 1L, id_draws = c("m1", "m2")) |>
      suppressWarnings() |>
      get_latent_draws() |>
      suppressWarnings() |>
      posterior::variables()
  }

  # the observation reads expl and rho_sum, so their producers are drawn too
  traced(observed) |>
    expect_setequal(
      paste(rep(c("m1", "m2"), each = 6),
            c("hidden", "hidden_2", "layer", "rho_sum", "alias", "expl"),
            sep = "_")
    )

  # with no observation nothing reads any of them: the model has no sink
  traced(hh_model) |> expect_length(0)
})

test_that("A latent degenerate domain warns and collapses", {
  # the bounds meet for m1 only: m1 is pinned, m2 is sampled
  # the pinned point is read in the located frame, loc + min
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                         datum = 1, datum_2 = c(1, 3), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", min = datum, max = datum_2,
               loc = 2) |>
    helper_flat_sink(loc = hidden)
  # the short run also skips the convergence diagnostics
  expect_warning(
    expect_warning(
      results <- hh_model |>
        run_inference(steps = 400L, chains = 1L, adapt = 50L),
      "[Cc]ollapsed .*latent .hidden. to constant"
    ),
    "[Ss]ample size was below"
  )
  draws <- as.data.frame(get_latent_draws(results))
  draws$m1_hidden |> expect_all_equal(3)
  length(unique(draws$m2_hidden)) |> expect_gt(1L)
  (draws$m2_hidden >= 3 & draws$m2_hidden <= 5) |> expect_all_true()

  # simulation draws the pinned member at its point too
  expect_warning(
    simulated <- run_simulation(hh_model, nrep = 3L),
    "[Cc]ollapsed .*latent .hidden. to constant"
  )
  vapply(simulated, \(x) x$hidden[1], numeric(1)) |> expect_all_equal(3)

  # both members collapse: one warning for the run
  hh_model_2 <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                           datum = 1) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", min = datum, max = datum,
               loc = 2)
  warnings <- capture_warnings(simulated <- run_simulation(hh_model_2,
                                                           nrep = 3L))
  warnings |> expect_length(1L)
  warnings |> expect_match("[Cc]ollapsed .*latent .hidden. to constant")
  unlist(lapply(simulated, `[[`, "hidden")) |> expect_all_equal(3)
})

test_that("Nested definitions works as expected", {
  hh_model <- new_housepi_model(model = "test") |>
    # define Box-Muller transformation
    add_latent(norm = list(model = "normal"),
               u1 = list(model = "uniform"),
               u2 = list(model = "uniform"),
               z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2))
  print(hh_model) |> expect_snapshot()

  # a flat observation on each chain leaves both posteriors equal to the model
  results <- hh_model |>
    add_data_wide(transform(hh_model$data_wide, obs = 0, obs_norm = 0)) |>
    helper_flat_sink(loc = z0) |>
    helper_flat_sink(loc = norm, column = "obs_norm") |>
    run_inference(steps = 1e5L, chains = 1L, adapt = 0L, thin = 1L)
  # no parameters nor priors, only latent
  posterior::nvariables(results$parameter_draws) |> expect_identical(0L)
  posterior::nvariables(results$latent_draws) |>
    expect_identical(dim(hh_model$latent)[1])

  # the S3 method of summarise_draws dispatches on the results
  out_vars <- posterior::summarise_draws(results)[["variable"]]
  startsWith(out_vars, "tommy") |> expect_all_true()

  # empirical tests for goodness of Gaussian distribution
  probs <- c(0.025, 0.5, 0.975)
  for (draws in list(results$latent_draws$tommy_norm,
                     results$latent_draws$tommy_z0)) {
    mean(draws) |> expect_equal(0, tolerance = 0.05)
    sd(draws) |> expect_equal(1, tolerance = 0.05)
    unname(stats::quantile(draws, probs)) |>
      expect_equal(qnorm(probs), tolerance = 0.05)
  }
})

test_that("Status-specific latent variables draw from per-status ranges", {
  hh_model <- data.frame(hh_id = 1, id = seq(4), status = seq(4), obs = 0) |>
    new_housepi_model("test") |>
    add_latent(
      hidden = list(status = 1, model = "uniform", min = 0, max = 1),
      hidden = list(status = 2, model = "uniform", min = 1, max = 2),
      hidden = list(status = 3, model = "uniform", min = 2, max = 3),
      hidden = list(status = 4, model = "uniform", min = 3, max = 4)
    ) |>
    helper_flat_sink(loc = hidden)

  # simulation routes draws into per-status supports
  simulated <- do.call(rbind, run_simulation(hh_model, nrep = 50))
  for (s in seq(4)) {
    draws <- simulated$hidden[simulated$status == s]
    (draws >= s - 1 & draws <= s) |> expect_all_true()
  }

  # inference traces every member within its status support
  results <- hh_model |>
    run_inference(steps = 5L, thin = 1L, adapt = 2L,
                  chains = 1L, id_draws = hh_model$data_wide$id) |>
    suppressWarnings()
  draws <- as.data.frame(get_latent_draws(results))
  for (s in seq(4)) {
    trace <- draws[[paste0(s, "_hidden")]]
    trace |> expect_length(5L)
    (trace >= s - 1 & trace <= s) |> expect_all_true()
  }
})

test_that("Evaluating a model allocates its latent storage", {
  hh_model <- data.frame(hh_id = 1, id = c("m1", "m2"), status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform")
  evaluate_log_lik(hh_model)
  augmented <- inspect_augmented_data(hh_model)
  is.data.frame(augmented) |> expect_true()
  ("hidden" %in% colnames(augmented)) |> expect_true()
})

##################################
###  OBSERVATION : END-TO-END  ###
##################################

test_that("Observation-aware initialization", {
  window <- 0.1
  # latent hidden ~ U(0, 1)
  # observation record ~ U(hidden - w, hidden + w)
  # if record = 0.5, valid for hidden in [0.4, 0.6]
  # without observation-aware restart, it would fail with 80% proba per-chain
  # probability of failing on K chains: 1 - 0.2^K
  hh_model <- data.frame(hh_id = "h1", status = 1, record = 0.5) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    add_observation("record", model = "uniform", loc = hidden,
                    min = -!!window, max = !!window)
  draws <- hh_model |>
    run_inference(steps = 10L, adapt = 0L, chains = 40L) |>
    suppressWarnings() |>
    get_latent_draws() |>
    as.data.frame()
  max(abs(draws[["1_hidden"]] - 0.5)) |> expect_lte(window)

  # datum record = 5 is unreachable
  hhdata <- data.frame(hh_id = "h1", status = 1, record = 5)
  hh_model |>
    add_data_wide(hhdata) |>
    run_inference() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  hhdata |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform") |>
    add_observation("record", model = "uniform", min = -rho, max = rho) |>
    run_inference() |>
    expect_error("[Ff]ree parameter dependencies")

  # generic error
  hhdata |>
    new_housepi_model("test") |>
    add_parameter("rho", value = 0.3) |>
    add_latent("hidden", model = "uniform") |>
    add_observation("record", model = "uniform", min = -rho, max = rho) |>
    run_inference() |>
    expect_error("[Oo]utside of respective domain")
})

test_that("An unbounded observation density stops the run", {
  # gamma with shape below 1 is unbounded at 0, reached at hidden = record
  hhdata <- data.frame(hh_id = "h1", status = 1, record = 5)
  mk <- function(...) {
    hhdata |>
      new_housepi_model("test") |>
      add_latent("hidden", ...) |>
      add_observation("record", model = "gamma", shape = 0.5, scale = 1,
                      loc = hidden)
  }

  mk(value = 5) |>
    evaluate_log_lik() |>
    expect_error("[Ll]og-likelihood is infinite")

  # a discrete latent reaches the bound within a few moves
  mk(model = "uniform", min = 0, max = 5, discr = TRUE) |>
    run_inference(steps = 200L, adapt = 0L) |>
    expect_error("[Ll]og-likelihood is infinite")
})

test_that("Members missing different observations still initialize", {
  # m1 carries a record and no record_2, m2 the reverse, m3 neither: the
  # ordinary censoring pattern must reach the sampler, not the
  # domain-feasibility wall of the observation-aware restart
  draws <- data.frame(
    id = c("m1", "m2", "m3"), hh_id = "h1", status = 1,
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
    get_latent_draws() |>
    as.data.frame()
  # each observed member stays inside its window
  max(abs(unlist(draws[c("m1_hidden", "m2_hidden")]) - 0.5)) |> expect_lte(0.1)
})

test_that("An observation infers discretization from the column it reads", {
  hh_model <- data.frame(hh_id = rep(c("h1", "h2"), each = 20),
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

  simulated <- run_simulation(hh_model)

  # dates and integers are discrete, doubles are not
  on_integer_lattice(simulated[["date_col"]]) |> expect_true()
  on_integer_lattice(simulated[["int_col"]]) |> expect_true()
  on_integer_lattice(simulated[["whole_col"]]) |> expect_false()
  on_integer_lattice(simulated[["kept_col"]]) |> expect_false()

  # a deterministic observation is not influenced
  hh_model |>
    add_observation(record = hidden + 1) |>
    run_simulation() |>
    getElement("record") |>
    on_integer_lattice() |>
    expect_false()
})

test_that("Simulated columns carry the discretization of their definition", {
  hh_model <- data.frame(hh_id = rep(c("h1", "h2"), each = 3),
                         id = seq(6), status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      record = list(model = "normal", discr = TRUE, mean = hidden, sd = 1),
      record_2 = list(model = "normal", mean = hidden, sd = 1),
      record_3 = list(value = min(record), level = "house")
    )

  simulated <- run_simulation(hh_model)

  # a discrete declaration, and the reduction folding it, export as integers
  simulated[["record"]] |> expect_type("integer")
  simulated[["record_3"]] |> expect_type("integer")
  # a continuous observation and a continuous latent stay double
  simulated[["record_2"]] |> expect_type("double")
  simulated[["hidden"]] |> expect_type("double")

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
  simulated <- data.frame(hh_id = rep(c("h1", "h2"), each = 2),
                          id = seq(4), status = rep(c(1, 2), 2)) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(
      record = list(status = 1, model = "normal", discr = TRUE, mean = hidden),
      record = list(status = 2, model = "normal", mean = hidden)
    ) |>
    run_simulation()

  # one column, two declarations: the continuous one has values to lose
  simulated[["record"]] |> expect_type("double")
  on_integer_lattice(simulated[["record"]][simulated$status == 2]) |>
    expect_false()
})

test_that("Discretization is resolved every time the data changes", {
  hhdata <- data.frame(hh_id = rep(c("h1", "h2"), 20), status = 1)
  hh_model <- new_housepi_model(model = "test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 2))

  # observations resolve to continuous in simulations and definitions
  hh_model <- hh_model |>
    add_data_wide(dplyr::mutate(hhdata, record = rep(c(2.4, 5.1), 20)))
  get_observations(hh_model)$discr |> expect_false()
  on_integer_lattice(run_simulation(hh_model)$record) |> expect_false()

  # observations resolve to discrete in simulations and definitions
  hh_model <- hh_model |>
    add_data_wide(dplyr::mutate(hhdata, record = rep(c(2L, 5L), 20)))
  get_observations(hh_model)$discr |> expect_true()
  on_integer_lattice(run_simulation(hh_model)$record) |> expect_true()
})

test_that("A bare observation name reads its column in a house reduction", {
  simulated <- mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -test_date)) |>
    add_observation(
      test_date = list(model = "normal", mean = infec_time + incub_period),
      start_date = list(value = min(test_date), level = "house")
    ) |>
    run_simulation()

  # the reduction runs over every infected member, not the first one
  infected <- simulated[simulated$status > 0, ]
  expected <- tapply(infected$test_date, infected$hh_id, min)
  per_house <- tapply(simulated$start_date, simulated$hh_id, unique)
  as.numeric(per_house[names(expected)]) |>
    expect_equal(as.numeric(expected))
})

test_that("A household-bound member observation bounds every member", {
  # the follow-up end reads the household start: a pinned index reads it too
  simulated <- mk_outbreak$data |>
    dplyr::select(-test_date, -end_date) |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-3, beta = 0.15, gamma = 1, hh_size_ref = 4) |>
    add_latent(
      incub_period = list(model = "gamma", mean = 2, sd = 1),
      infec_time = list(value = 0, min = start_date - 8, max = end_date,
                        proposal = "rw")
    ) |>
    add_observation(
      test_date = list(model = "normal", mean = infec_time + incub_period),
      start_date = list(value = min(test_date), level = "house"),
      end_date = start_date + 21,
      start_date_2 = list(value = min(test_date + 1), level = "house")
    ) |>
    run_simulation()

  simulated$end_date |> expect_no_na()
  simulated$end_date |> expect_equal(simulated$start_date + 21)
  # uninfected members leave the risk set at the end of follow-up
  mean(simulated$status > 0) |> expect_lt(1)
  # the argument is computed for each member: an uninfected one draws no test
  # date and drops out of the fold
  as.numeric(simulated$start_date_2) |>
    expect_equal(as.numeric(simulated$start_date) + 1)
})

test_that("An observation model reads an observation drawn before it", {
  # the reader is declared first: the draw order follows the dependency
  simulated <- mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -test_date)) |>
    add_observation(
      record = list(model = "normal", mean = test_date, sd = 1e-3),
      test_date = list(model = "normal", mean = infec_time + incub_period)
    ) |>
    run_simulation()

  infected <- simulated[simulated$status > 0, ]
  max(abs(infected$record - infected$test_date)) |> expect_lt(0.1)
})

test_that("Observation cyclic dependencies are captured", {
  data.frame(hh_id = "h1", status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "lognormal"),
               hidden_2 = list(model = "uniform", max = record_2)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1)) |>
    add_observation("record_2", value = min(record), level = "house") |>
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
  max(abs(recon - observed)) |> expect_lt(0.5)
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

  get_log_lik_draws(results)$h1 |> expect_all_finite()
})

test_that("A stochastic observation is simulated for infected members only", {
  simulated <- mk_outbreak$model |>
    add_latent(onset_time = infec_time + incub_period) |>
    add_observation("onset_date", model = "normal", mean = onset_time) |>
    run_simulation()

  # onset_date is generated only for infected members
  mask_infected <- simulated$status > 0
  simulated$onset_date[mask_infected] |> expect_no_na()
  simulated$onset_date[!mask_infected] |> expect_all_na()
})

test_that("Household observation derives from simulated observations", {
  # both the observation (test_date) and the surrogate column (start_date) are
  # produced internally: neither is provided in the input data
  reps <- mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -test_date)) |>
    add_observation("test_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1) |>
    add_observation(start_date = min(test_date), level = "house") |>
    run_simulation(nrep = 2)

  # each repetition is self-consistent: start_date is the per-household na.rm
  # min of that repetition's generated test_date, one value per household
  # replicated across rows, never NA (every house has a case)
  for (simulated in reps) {
    simulated$start_date |> expect_no_na()
    per_house <- tapply(simulated$start_date, simulated$hh_id, unique)
    lengths(per_house) |> as.vector() |> expect_all_equal(1)
    infected <- simulated[simulated$status > 0, ]
    expected <- tapply(infected$test_date, infected$hh_id, min)
    as.numeric(unlist(per_house)[names(expected)]) |>
      expect_equal(as.numeric(expected))
  }
})

test_that("A household reduction folds a latent over the members", {
  # the escaped member draws no 'hidden' and drops out of the fold
  hh_model <- data.frame(hh_id = "h1", status = c(1, 1, 0),
                         datum = c(2, 9, 4)) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1))

  results <- hh_model |>
    add_observation(record_2 = max(hidden), record_3 = sum(hidden + 1),
                    level = "house") |>
    run_simulation()
  unique(results$record_2) |>
    expect_equal(max(results$hidden, na.rm = TRUE))
  unique(results$record_3) |>
    expect_equal(sum(results$hidden + 1, na.rm = TRUE))

  # a latent outside any reduction has no household value
  hh_model |>
    add_observation("record_2", value = hidden + 0, level = "house") |>
    run_simulation() |>
    expect_error("[Pp]er-member latent .hidden. is read outside a reduction")
})

test_that("A column needed for every member is not drawn with infection", {
  hh_model <- mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -end_date)) |>
    add_latent(onset_time = infec_time + incub_period) |>
    add_observation("onset_date", model = "normal", mean = onset_time)

  # onset_date has no value for escaped members, nor has what reads it
  hh_model |>
    add_observation(end_date = onset_date + 21) |>
    run_simulation() |>
    expect_error("[Dd]ata .end_date. needs a value for every member")

  # a household value covers escaped members too
  results <- hh_model |>
    add_observation(end_date = min(onset_date) + 21, level = "house") |>
    run_simulation()
  results$end_date |> expect_no_na()
})

test_that("A model argument reading an empty cell names member and household", {
  sms <- data.frame(hh_id = "h1", status = 1, datum = c(2, NA, 4),
                    record = 0) |>
    new_housepi_model("test") |>
    add_observation(record = list(model = "normal", mean = datum + 1)) |>
    run_simulation() |>
    tryCatch(error = conditionMessage)
  sms |> expect_match("[Dd]ata .datum. has no value")
  sms |> expect_match("observation .record. for\\s+individual\\s+.2.")
  sms |> expect_match(".2.\\s+\\(status\\s+1\\)")
  sms |> expect_match("\\(status\\s+1\\)\\s+for\\s+household\\s+.h1.")
})

test_that("Household observations aggregate a data column across members", {
  # one household carries no test_date at all: its na.rm reduction is empty
  empty_house <- mk_outbreak$data$hh_id[[1]]
  hhdata <- dplyr::mutate(
    mk_outbreak$data,
    test_date = dplyr::if_else(hh_id == empty_house, as.Date(NA), test_date)
  )

  simulated <- mk_outbreak$model |>
    add_data_wide(hhdata) |>
    add_observation(
      first_test = min(test_date),
      last_test = max(test_date),
      level = "house"
    ) |>
    run_simulation()

  # one value per household, replicated across every member row
  per_house <- function(col) {
    values <- tapply(simulated[[col]], simulated$hh_id, unique)
    lengths(values) |> as.vector() |> expect_all_equal(1)
    unlist(values)[sort(names(values))]
  }
  # values equal the na.rm min/max of test_date over each household's members,
  # NA when every member is missing
  expected <- function(f) {
    reduce <- \(x) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE)
    values <- tapply(as.numeric(hhdata$test_date), hhdata$hh_id, reduce)
    values[sort(names(values))]
  }
  per_house("first_test") |> expect_equal(expected(min))
  per_house("last_test") |> expect_equal(expected(max))
  per_house("first_test")[[empty_house]] |> expect_all_na()
})

test_that("Simulation leaves overwritten data columns pristine", {
  # test_date is provided data and a generated observation: simulation
  # overwrites it with generated values, then must restore the input
  hh_model <- mk_outbreak$model |>
    add_observation("test_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1)

  before <- inspect_augmented_data(hh_model)
  invisible(run_simulation(hh_model, nrep = 2))
  after <- inspect_augmented_data(hh_model)

  before <- before[order(before$id), ]
  after <- after[order(after$id), ]
  # the overwritten column is restored, and no simulation-added column lingers
  as.numeric(after$test_date) |> expect_equal(as.numeric(before$test_date))
  names(after) |> expect_setequal(names(before))
})

test_that("A latent reads a static data-bound household observation", {
  # start_date stands in for provided (static) test_date, read by infec_time
  simulated <- mk_outbreak$model |>
    add_data_wide(mk_outbreak$data |>
                    dplyr::mutate(test_date = 5, end_date = 26)) |>
    add_latent(
      infec_time = list(value = start_date - 1, min = start_date - 5,
                        max = start_date + 5, proposal = "rw"),
      append = FALSE
    ) |>
    add_observation(start_date = min(test_date), level = "house") |>
    run_simulation()

  # start_date is materialized statically (= 5) before the index is pinned
  simulated$start_date |> expect_all_equal(5)
  tapply(simulated$infec_time, simulated$hh_id, min, na.rm = TRUE) |>
    as.vector() |>
    expect_all_equal(4)
})

test_that("Simulation rejects a latent reading a live household observation", {
  # simulation value of infec_time depends on start_date,
  # but start_date is consumer of test_date: cycle get rejected
  mk_outbreak$model |>
    add_data_wide(dplyr::select(mk_outbreak$data, -test_date)) |>
    add_latent(
      infec_time = list(value = start_date - 1, min = start_date - 5,
                        max = start_date + 5, proposal = "rw"),
      append = FALSE
    ) |>
    add_observation("test_date", model = "normal", mean = infec_time, sd = 1) |>
    add_observation(start_date = min(test_date), level = "house") |>
    run_simulation() |>
    expect_error("[Cc]yclic")
})

test_that("Data-bound household observation supplies a required data column", {
  # the dates of the dataset itself, where end_date is 21 days after the first
  # test date
  hhdata <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 40))
  hh_model <- new_housepi_model(model = "case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")

  # 'end_date' provided as an ordinary data column
  provided <- add_data_wide(hh_model, hhdata)

  # 'end_date' is passed as a surrogate observation
  surrogate <- hh_model |>
    add_data_wide(dplyr::select(hhdata, -end_date)) |>
    add_observation(end_date = min(test_date) + 21, level = "house")

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

test_that("Status process drives status-specific latent blocks", {
  simulated <- mk_outbreak$model |>
    # status-specific block: symptomatic incubation is longer
    add_latent(status = "symptomatic",
               incub_period = list(model = "gamma", shape = 8, scale = 1)) |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  by_status <- tapply(simulated$incub_period,
                      as.character(simulated$status), mean, na.rm = TRUE)
  # the resolved status selects the block: shape 8 symptomatic, 2 asymptomatic
  by_status[["symptomatic"]] |> expect_gt(5)
  by_status[["asymptomatic"]] |> expect_lt(4)
})

test_that("Escaped label is configurable", {
  simulated <- mk_outbreak$data |>
    new_housepi_model("case-ascertained", escaped = "healthy") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-12, beta = 0.001, gamma = 1, hh_size_ref = 4) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  # weak transmission leaves escapees, labelled per the construction argument
  ("healthy" %in% levels(simulated$status)) |> expect_true()
  any(as.character(simulated$status) == "healthy") |> expect_true()
})

test_that("Status-specific blocks alone cover every drawn status", {
  # simulation draws the status, so the data carries none: the blocks a drawn
  # status can land on are the declared ones, never the observed ones
  mk <- function(..., weights = c(0.5, 0.5)) {
    data.frame(
      id = c("m1", "m2", "m3", "m4"), hh_id = c("h1", "h1", "h2", "h2"),
      status = NA, end_date = 21
    ) |>
      new_housepi_model("case-ascertained") |>
      add_parameter(import = "template") |>
      add_infectivity_profile(import = "template") |>
      add_status_process(symptomatic = weights[1], asymptomatic = weights[2]) |>
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
  suppressWarnings(run_simulation(covered)) |> expect_s3_class("data.frame")

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
  sms |> expect_match("status 2")
  grepl("generic latent definition", sms) |> expect_false()

  # a generic block on its own still serves every drawn status
  generic <- mk(
    infec_time = 0,
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14),
    onset = list(func = infec_time + incub)
  )
  suppressWarnings(run_simulation(generic)) |> expect_s3_class("data.frame")
})
