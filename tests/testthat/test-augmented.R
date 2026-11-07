library(testthat)

####  Helpers  ####

# model carrying the data columns the definition blocks refer to
mk_latent <- local({
  hh_model <- data.frame(hh_id = "h1", status = 1, datum = 1,
                         record = 3, datum_2 = 2) |>
    new_housepi_model("test")
  list(model = hh_model)
})

# saturated outbreak: every household infects all of its members
mk_outbreak <- local({
  hh_data <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 50)) |>
    dplyr::mutate(end_date = 21)
  hh_model <- hh_data |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(alpha = 1e-12, beta = 100, gamma = 1, hh_size_ref = 4) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = list(min = -14, max = end_date)
    )
  list(data = hh_data, model = hh_model)
})

# case-ascertained household sampling incubation periods and infection times
mk_onset <- local({
  hh_data <- dplyr::select(helper_onset_house$data, -infec_time)
  hh_model <- hh_data |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1) |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1, max = 14),
      infec_time = list(min = -14, max = end_date)
    ) |>
    add_infectivity_profile(func = "gamma", shape = shape, scale = scale)
  list(data = hh_data, model = hh_model)
})

# two infected statuses under the template definitions
mk_status <- local({
  hh_data <- data.frame(id = c("m1", "m2", "m3", "m4"), hh_id = "h1",
                        status = c(1, 2, 0, 0), end_date = 21)
  hh_model <- hh_data |>
    new_housepi_model("case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")
  list(data = hh_data, model = hh_model)
})

# whether on an integer lattice
on_integer_lattice <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  length(x) > 0 && all(x == floor(x))
}

# draws of a run as a data frame of its latent traces
latent_frame <- function(x) {
  helper_draws_to_frame(get_latent_draws(x))
}

####  Latent definitions  ####

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
  # a uniform model spread to infinity holds no mass to draw from
  mk_latent$model |>
    add_latent(hidden = list(model = "uniform", min = -Inf, max = 1.1)) |>
    expect_error("min.*uniform.*must be finite.*Change .min. \\(-Inf\\)")
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
    add_latent("hidden", const = NA) |>
    expect_error(".const. of latent .hidden. must be the value")
  mk_latent$model |>
    add_latent("hidden", const = 3, proposal = "forward") |>
    expect_error("[Cc]onstant latent .hidden. cannot be sampled")
})

test_that("A latent without a model is sampled within its bounds", {
  # bounds alone sample, no proposal declared
  latent <- mk_status$model |>
    add_latent(infec_time = list(min = -1, max = end_date), append = FALSE) |>
    get_latent()
  latent$proposal |> expect_all_na()
  latent$const |> expect_all_na()
  latent$model |> expect_all_na()

  # a constant has no distribution to bound
  mk_status$model |>
    add_latent(infec_time = list(const = 3, min = 0, max = 5),
               append = FALSE) |>
    expect_error("[Cc]onstant latent .infec_time. has no distribution to bound")
})

test_that("A latent without a model rejects incomplete definitions", {
  # a forward proposal draws from the model
  mk_latent$model |>
    add_latent(hidden = list(min = 0, max = 5, proposal = "forward")) |>
    expect_error("[Nn]o .model. to sample from")
  # a bare latent needs both bounds
  mk_latent$model |>
    add_latent(hidden = list(min = 0)) |>
    expect_error("[Nn]o .model. to sample from")
})

test_that("Passing a prior raises a hint on model", {
  mk_latent$model |>
    add_latent(hidden = list(prior = "uniform", min = -1, max = 2)) |>
    expect_error("[Ss]pare.*prior.*hidden.*[Dd]id you mean.*model")
})

test_that("A func definition is deterministic", {
  # a plain number is a constant value, whichever argument carried it
  by_func <- mk_latent$model |> add_latent(hidden = list(func = 21))
  by_const <- mk_latent$model |> add_latent(hidden = list(const = 21))
  get_latent(by_func) |> expect_identical(get_latent(by_const))

  # deterministic definitions reject model or proposal
  mk_latent$model |>
    add_latent("hidden", func = "log", x = 2,
               model = "normal", mean = 0, sd = 1) |>
    expect_error("only one.*func.*model")
  mk_latent$model |>
    add_latent("hidden", func = "log", x = 2, proposal = "rw")   |>
    expect_error("[Dd]eterministic.*hidden.*cannot be sampled")
  # nor a constant alongside the function
  mk_latent$model |>
    add_latent("hidden", const = 2, func = "log", x = 1) |>
    expect_error("[Oo]nly one of .const. and .func.")
})

test_that("A name is resolved, not guessed", {
  # a name two holders answer to is never read
  mk_latent$model |>
    add_parameter(datum = list(prior = "gamma", shape = 2, scale = 1)) |>
    add_observation(record = list(model = "normal", mean = datum, sd = 1)) |>
    evaluate_log_lik() |>
    expect_error(paste0("[Cc]ould not resolve the source of .datum.",
                        ".*[Dd]eclared both as data column and parameter",
                        ".*[Rr]ename one of them"))

  # nothing of that name anywhere: the defect is named, with the near miss
  mk_latent$model |>
    add_observation(record = list(model = "normal", mean = datums, sd = 1)) |>
    evaluate_log_lik() |>
    expect_error(paste0("[Cc]ould not resolve the source of .datums.",
                        ".*[Dd]id you mean .datum."))
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
    add_latent("hidden", const = 3) |>
    get_latent()
  augm$name |> expect_identical("hidden")
  augm$status |> expect_equal(NA_real_)
  augm$const |> expect_equal(3)
  augm$discr |> expect_identical(FALSE)
  augm$proposal |> expect_equal(NA_real_)
  # a bare number is the same constant
  get_latent(add_latent(mk_latent$model, hidden = 3)) |> expect_identical(augm)

  # a model with no constant is sampled, no proposal declared
  augm <- mk_latent$model |>
    add_latent("hidden", model = "normal") |>
    get_latent()
  augm$const |> expect_all_na()
  augm$proposal |> expect_all_na()
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
    add_latent("hidden", const = datum, model = "uniform",
               min = layer_1, max = rho) |>
    get_latent("hidden")
  augm$status |> expect_equal(NA_real_)
  augm$const |> expect_identical("datum|D")
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
    expect_error("support \\[0, Inf\\).*cannot cover.*domain \\(-1, 5\\]")

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

test_that("A latent random walk adapts where its scale reads NA", {
  # the scale is read per member: m2 has none and adapts, m1 keeps its own
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                         datum = c(0.5, NA), record = 1) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", min = 0, max = 5,
               proposal = "rw", proposal_scale = datum) |>
    helper_flat_sink(loc = hidden, column = "record")
  hh_model |>
    run_inference(steps = 10L, thin = 1L, adapt = 200L, chains = 1L) |>
    suppressWarnings()

  moves <- inspect_state(hh_model, "h1")$.sampler$moves
  scale <- stats::setNames(moves$scale, moves$member)
  scale[["m1"]] |> expect_identical(0.5)
  # 1 is the placeholder a walk holds before its first reset
  (scale[["m2"]] == 1) |> expect_false()
})

test_that("Latent function definition is correct", {
  hh_model <- mk_latent$model |>
    add_latent(layer_1 = list(model = "normal", mean = 0, sd = 1),
               layer_2 = list(model = "gamma", shape = 2, scale = 1),
               layer_3 = list(model = "lognormal", mean = 1, sd = 1))

  # deterministic definition: subtraction (as sum of minuses)
  augm <- hh_model |>
    add_latent("hidden", func = layer_1 - layer_2 - layer_3) |>
    get_latent("hidden")
  augm$const |> expect_all_na()
  augm$func |> expect_identical("+")
  augm$func_args[[1]] |>
    expect_identical(c(x = "layer_1|A",
                       x2 = "-|x|[layer_2|A]",
                       x3 = "-|x|[layer_3|A]"))

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

  get_latent(hh_model_2, "hidden")$const |>
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
               layer_1 = list(func = layer_2, status = c(2, 3, 4)),
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
               layer_1 = list(func = layer_2, status = 2),
               hidden = layer_1) |>
    evaluate_log_lik() |>
    expect_error("[Nn]o.*definition.*status")

  # second generation of definitions: no layer_2 for status 4
  hh_model |>
    add_latent(layer_2 = list(model = "lognormal", status = c(2, 3)),
               layer_1 = list(model = "normal", status = 1),
               layer_1 = list(func = layer_2, status = c(2, 3, 4)),
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
  hh_model |>
    add_latent("hidden", model = "uniform", min = 0.5, max = 3, discr = TRUE) |>
    evaluate_log_lik() |>
    expect_error("min.*uniform.*discrete.*hidden")
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
  # an expression over a sampled declaration shifts it alike
  hh_model |>
    add_parameter("rho", prior = "uniform", min = 0, max = 1) |>
    add_latent("hidden", model = "poisson", lambda = 2, loc = rho * 2) |>
    evaluate_log_lik() |>
    expect_error("loc.*incompatible with discrete")
})

test_that("An argument read from a source is judged on every value it takes", {
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                         record = 1.5, datum = c(2, -1), datum_2 = c(0, 2),
                         obs = 0) |>
    new_housepi_model("test")

  # a prior reaching zero or below cannot set a standard deviation
  hh_model |>
    add_parameter(rho = list(prior = "normal")) |>
    add_observation("record", model = "normal", sd = rho) |>
    evaluate_log_lik() |>
    expect_error("sd.*normal.*must be greater than 0.*rho.*can take other")
  hh_model |>
    add_parameter(rho = list(prior = "exponential")) |>
    add_observation("record", model = "normal", sd = rho) |>
    evaluate_log_lik() |>
    expect_all_finite()

  # a data cell is judged on each member
  hh_model |>
    add_observation("record", model = "normal", sd = datum) |>
    evaluate_log_lik() |>
    expect_error("sd.*normal.*must be greater than 0.*datum.*-1.*m2")

  # an expression reaching zero or below is refused, live or constant
  hh_model |>
    add_parameter(rho = list(prior = "exponential")) |>
    add_observation("record", model = "normal", sd = rho - 1) |>
    evaluate_log_lik() |>
    expect_error("sd.*normal.*must be greater than 0")
  hh_model |>
    add_observation("record", model = "normal", sd = datum_2 - 1) |>
    evaluate_log_lik() |>
    expect_error("sd.*normal.*must be greater than 0")

  # an expression is judged on the values its operands take
  with_rho <- hh_model |> add_parameter(rho = list(prior = "exponential"))
  with_rho |>
    add_latent("hidden", model = "binomial", size = 3,
               prob = 1 - exp(-rho * (datum_2 + 1))) |>
    helper_flat_sink(loc = hidden) |>
    evaluate_log_lik() |>
    expect_all_finite()
  with_rho |>
    add_latent("hidden", model = "binomial", size = 3,
               prob = 1 / (1 + exp(-rho))) |>
    helper_flat_sink(loc = hidden) |>
    evaluate_log_lik() |>
    expect_all_finite()
  # a covariate at zero sets the probability to zero
  with_rho |>
    add_latent("hidden", model = "binomial", size = 3,
               prob = 1 - exp(-rho * datum_2)) |>
    helper_flat_sink(loc = hidden) |>
    evaluate_log_lik() |>
    expect_error("prob.*binomial.*must be greater than 0 and less than 1")
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
    add_latent(hidden = list(model = "normal", const = datum),
               hidden_2 = list(model = "normal", const = datum_2)) |>
    get_latent()

  augm$const |> expect_identical(c("datum|D", "datum_2|D"))
  # the date column is discrete, the double one is not
  augm$discr |> expect_identical(c(TRUE, FALSE))
})

####  Observation definitions  ####

test_that("Observation definitions are encoded", {
  obs <- mk_latent$model |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(
      record = list(model = "normal", mean = hidden, sd = 1)
    ) |>
    get_observations()

  # marked as observation, data-bound and never sampled
  obs[["level"]] |> expect_identical("member")
  obs[["const"]] |> expect_null()
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
  obs[["func"]] |> expect_match("floor")
})

test_that("Observation node rejects conflicting options", {
  # observation nodes are data-bound, never sampled: the options are dropped
  ignored <- function(...) {
    expect_warning(hh_model <- add_observation(mk_latent$model, ...),
                   "[Ii]gnored.*const.*proposal")
    get_observations(hh_model)
  }
  obs <- ignored("record", model = "normal", const = 3,
                 proposal = "forward")
  obs[["const"]] |> expect_null()
  obs[["proposal"]] |> expect_all_na()

  # a deterministic observation drops them alike
  obs <- ignored(record_2 = list(func = min(datum_2), status = 1,
                                 const = 3, proposal = "rw"))
  obs[["func"]] |> expect_identical("min")
  obs[["const"]] |> expect_null()
  obs[["proposal"]] |> expect_all_na()

  # a decorator carrying no definition at all
  mk_latent$model |>
    add_observation() |>
    expect_error("[Mm]issing or empty observation definitions")
  # an entry with neither model nor function
  mk_latent$model |>
    add_observation(record = list(status = 1)) |>
    expect_error("[Ii]ncomplete definition for observation .record.")
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
    add_latent(.status = list(const = 1)) |>
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

  # a deterministic observation standing in for a like-named data column is by
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
  names(mixed) |> expect_disjoint("args")

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
    expect_error("[Uu]nknown .level. .unknown")

  # 'level' is meaningless outside observations
  mk_latent$model |>
    add_latent(hidden = list(const = 1, level = "house")) |>
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
      add_latent("hidden", const = 0.5) |>
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
    add_latent("hidden", const = 0.5) |>
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
    hh_data <- data.frame(id = c("m1", "m2"), hh_id = "h1",
                          status = 1, record = x)
    if (!is.null(x_2)) hh_data$record_2 <- x_2
    hh_model <- hh_data |>
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
    add_latent("hidden", func = record) |>
    add_observation("record", model = "uniform", min = hidden, max = 10) |>
    expect_error("[Cc]yclic")
})

test_that("A deterministic observation reading a latent is simulation only", {
  hh_model <- data.frame(id = seq_len(10), hh_id = "h1",
                         status = c(rep(1, 9), 0), record = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record = hidden + 2)

  # the data view has no computed value to report: it keeps the imported datum
  get_dataset(hh_model)[["record"]] |> expect_all_equal(0)

  # observation reads a sampled latent: not a static deterministic observation
  evaluate_log_lik(hh_model) |>
    expect_error("must resolve from data alone")

  # simulation produces it for infected members, following the latent draw
  simulated <- run_simulation(hh_model)
  infected <- simulated$status > 0
  sum(names(simulated) == "record") |> expect_identical(1L)
  simulated$record[infected] |> expect_equal(simulated$hidden[infected] + 2)
  simulated$record[!infected] |> expect_all_na()
})

test_that("A deterministic observation overwrites the like-named data column", {
  hh_model <- data.frame(hh_id = "h1", status = c(1, 0, 0),
                         datum_2 = c(2, 6, NA), record_2 = rep(99, 3)) |>
    new_housepi_model("test") |>
    add_observation(
      record_2 = list(func = min(datum_2), level = "house"),
      record_3 = datum_2 + 1
    ) |>
    add_latent(hidden = list(func = record_3 * 2))

  # the data view reports the computed value before evaluation
  unique(get_dataset(hh_model)[["record_2"]]) |> expect_identical(2)
  colnames(get_dataset(hh_model)) |> expect_disjoint("record_3")

  # and after evaluation
  evaluate_log_lik(hh_model) |> expect_all_finite()
  unique(get_dataset(hh_model)[["record_2"]]) |> expect_identical(2)
  colnames(get_dataset(hh_model)) |> expect_disjoint("record_3")

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

test_that("Deterministic observations build producer before consumer", {
  hh_data <- data.frame(hh_id = "h1", status = c(1, 0, 0),
                        datum = c(2, 9, 4), record = 0, record_2 = 0)
  hh_model <- hh_data |>
    new_housepi_model("test") |>
    add_observation(record = datum + 3,
                    record_2 = list(func = max(record), level = "house"))

  # every member carries a value regardless of status
  hh_data_2 <- get_dataset(hh_model)
  hh_data_2[["record"]] |> expect_identical(c(5, 12, 7))
  unique(hh_data_2[["record_2"]]) |> expect_identical(12)

  # deterministic observations with no like-named data column are not exported
  hh_model |>
    add_dataset(hh_data[c("hh_id", "status", "datum")]) |>
    get_dataset() |>
    colnames() |>
    intersect(c("record", "record_2")) |>
    expect_length(0L)
})

test_that("Reading only house observations casts member level to house", {
  hh_data <- data.frame(hh_id = "h1", status = c(1, 0, 0), datum = c(2, 9, 4),
                        record = 0, record_2 = 0, record_3 = 0,
                        record_4 = 0, record_5 = 0, record_6 = 0)
  hh_model <- hh_data |>
    new_housepi_model("test") |>
    add_observation(record = list(func = max(datum), level = "house"),
                    record_2 = record + 1,
                    record_3 = record_2 * 2,
                    record_4 = datum + 1,
                    record_5 = record - 1,
                    record_6 = 18 / record)

  # a cast row is a household value to the rows reading it
  obs <- get_observations(hh_model)
  names <- c("record_2", "record_3", "record_4", "record_5", "record_6")
  obs$level[match(names, obs$name)] |>
    expect_identical(c("house", "house", "member", "house", "house"))

  hh_data_2 <- get_dataset(hh_model)
  hh_data_2[["record_2"]] |> expect_equal(c(10, 10, 10))
  hh_data_2[["record_3"]] |> expect_equal(c(20, 20, 20))
  hh_data_2[["record_4"]] |> expect_equal(c(3, 10, 5))
  hh_data_2[["record_5"]] |> expect_equal(c(8, 8, 8))
  hh_data_2[["record_6"]] |> expect_equal(c(2, 2, 2))
})

test_that("A household reduction folds its argument over the members", {
  # datum differs across members, the first one holds 2
  hh_data <- data.frame(hh_id = "h1", status = c(1, 0, 0),
                        datum = c(2, 9, 4), record = 0)
  reduce <- function(defin, data = hh_data) {
    data |>
      new_housepi_model("test") |>
      add_observation("record", func = !!defin, level = "house") |>
      get_dataset() |>
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
  reduce(quote(max(datum + 1)), transform(hh_data, datum = c(2, NA, 4))) |>
    expect_equal(5)
  # with no datum at all there is nothing to fold
  reduce(quote(max(datum) + 1), transform(hh_data, datum = NA_real_)) |>
    expect_all_na()

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
  hh_data <- data.frame(id = c("m1", "m2", "m3", "m4"),
                        hh_id = c("h1", "h1", "h2", "h2"),
                        status = NA) |>
    new_housepi_model("case-ascertained") |>
    add_observation(
      onset_date = list(model = "normal", discr = TRUE, mean = 3, sd = 1),
      end_date = list(func = max(onset_date), level = "house")
    ) |>
    get_dataset()

  hh_data$end_date |> expect_all_na()
})

test_that("A deterministic observation proves the requirement its model asks", {
  # a deterministic observation filling 'end_date' answers the discreteness
  # 'case_ascertained' asks of the column, and answers for what it can return,
  # not for the value it happens to settle on.
  # A household-level definition reduces over its members,
  # so its root is the reduction; the member level takes the algebra directly
  mk <- function(defin, level = "house") {
    data.frame(id = c("m1", "m2", "m3", "m4"),
               hh_id = c("h1", "h1", "h2", "h2"),
               status = c(1, 0, 1, 0), datum = 0.5, datum_2 = 16) |>
      new_housepi_model("case-ascertained") |>
      add_observation("end_date", func = !!defin, level = level) |>
      get_dataset()
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
  mk(quote(sum(datum))) |> expect_error("[Hh]ousehold observation.*h1")
  mk(quote(datum_2 / 4), level = "member") |> expect_error("individual.*m1")
})

test_that("Data-bound household observation rejects latent dependence", {
  # a household observation standing in for data resolves from data alone
  mk_outbreak$model |>
    add_latent("infec_time", min = -14, max = end_date, append = FALSE) |>
    add_observation(start_date = min(test_date) + alpha,
                    level = "house") |>
    evaluate_log_lik() |>
    expect_error("must resolve from data alone")
})

####  Status definitions  ####

test_that("Escaped label colliding with an infected status is rejected", {
  mk_status$data |>
    new_housepi_model("case-ascertained", escaped = "symptomatic") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    expect_error("collides with the escaped-status label")
})

test_that("Status declarations are all integers or all declared labels", {
  # one call carries both an integer and a character status: caught at the wall
  mk_status$model |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1,
                          status = "symptomatic"),
      onset = list(const = 1, status = 1)
    ) |>
    evaluate_log_lik() |>
    expect_error("[Mm]ix")
  # a label the status process does not declare points to the close one
  mk_status$model |>
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    add_latent(incub_period = list(model = "gamma", shape = 2, scale = 1,
                                   status = "symptomatc")) |>
    evaluate_log_lik() |>
    expect_error("[Uu]nknown status label .symptomatc.*[Dd]id you mean")
})

####  Latent end-to-end  ####

test_that("Independence proposal inherits the model location shift", {
  hh_fit <- data.frame(id = c("m1", "m2"), hh_id = "h1",
                       status = 1, datum = c(4, 8), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", loc = datum,
               proposal = "gamma", proposal_shape = 2) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  draws <- as.data.frame(get_latent_draws(hh_fit))
  (draws$m1_hidden >= 4 & draws$m1_hidden <= 5) |> expect_all_true()
  (draws$m2_hidden >= 8 & draws$m2_hidden <= 9) |> expect_all_true()
})

test_that("Latent variables moving together sample their models", {
  # no proposal declared: latents move jointly, `layer_1` on its log coordinate
  # `hidden` reads `layer_1` as its mean, a ridge of width 0.1
  # a lost Jacobian draws `layer_1` from a gamma of shape 3, mean 1.5
  hh_fit <- data.frame(id = "m1", hh_id = "h1", status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "gamma", shape = 4, scale = 0.5),
               hidden = list(model = "normal", mean = layer_1, sd = 0.1)) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 2e5L, adapt = 500L, thin = 5L, chains = 2L) |>
    suppressWarnings()

  draws <- latent_frame(hh_fit)
  mean(draws$m1_layer_1) |> expect_equal(2, tolerance = 0.1)
  sd(draws$m1_hidden - draws$m1_layer_1) |> expect_equal(0.1, tolerance = 0.05)
})

test_that("A latent bounded through a function of another moves alone", {
  # `hidden` reads `layer_1` through `shifted`:
  # moving them together biases `layer_1` to a mean above 3
  hh_fit <- data.frame(id = "m1", hh_id = "h1", status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent(
      layer_1 = list(model = "gamma", shape = 4, scale = 0.5),
      shifted = list(func = layer_1 + 0),
      hidden = list(model = "uniform", min = shifted, max = shifted + 1)
    ) |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 40000L, adapt = 500L, thin = 1L, chains = 2L) |>
    suppressWarnings()

  mean(latent_frame(hh_fit)$m1_layer_1) |> expect_equal(2, tolerance = 0.2)
})

test_that("Truncated latent models are normalized in log-likelihood traces", {
  # bounds read data, constants and a sampled parameter
  hh_data <- data.frame(id = c("m1", "m2"), hh_id = c("h1", "h2"),
                        status = 1, datum = c(3, 6), obs = 0)
  hh_fit <- hh_data |>
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

  draws <- latent_frame(hh_fit)
  rho <- helper_draws_to_frame(get_parameter_draws(hh_fit))$rho
  log_lik <- helper_draws_to_frame(get_log_lik_draws(hh_fit))
  log_dtrunc <- function(x, shape, scale, min, max) {
    dgamma(x, shape, scale = scale, log = TRUE) -
      log(pgamma(max, shape, scale = scale) - pgamma(min, shape, scale = scale))
  }
  # one member per household, the flat sink adds log(1 / 2e6)
  for (k in seq(2)) {
    col <- \(name) draws[[sprintf("m%d_%s", k, name)]]
    expected <- log_dtrunc(col("layer_1"), 2, 1.5, 0, hh_data$datum[k]) +
      log_dtrunc(col("hidden"), rho, 1, 0.5, 8) +
      log_dtrunc(col("hidden_2"), 1.7, 1, 0, 5) - log(2e6)
    log_lik[[sprintf("h%d", k)]] |> expect_equal(expected, tolerance = 1e-8)
  }
})

test_that("Latent draws stay inside their live support", {
  hh_data <- data.frame(id = seq(2), hh_id = "h1", status = 1, obs = 0)
  run_short <- function(hh_model) {
    hh_model |>
      run_inference(steps = 100L, adapt = 0L, thin = 1L, chains = 1L) |>
      suppressWarnings()
  }

  # other latent variable
  draws <- hh_data |>
    new_housepi_model("test") |>
    add_latent(layer_1 = list(model = "uniform"),
               hidden = list(model = "uniform", loc = layer_1)) |>
    helper_flat_sink(loc = hidden) |>
    run_short() |>
    latent_frame()
  gap <- draws[["1_hidden"]] - draws[["1_layer_1"]]
  (gap >= 0 & gap <= 1) |> expect_all_true()

  # transitive deterministic latent variable
  draws <- hh_data |>
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
  hh_fit <- hh_data |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform") |>
    add_latent("hidden", model = "uniform", min = rho, max = 5) |>
    helper_flat_sink(loc = hidden) |>
    run_short()
  hidden <- latent_frame(hh_fit)[["1_hidden"]]
  rho <- as.data.frame(get_parameter_draws(hh_fit))[["rho"]]
  (hidden >= rho) |> expect_all_true()
})

test_that("Draws requested where no latent is sampled are ignored", {
  # a member carrying no infected status draws no latent
  data.frame(id = c("m1", "m2"), hh_id = "h1", status = c(1, 0), obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    helper_flat_sink(loc = hidden) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  id_draws = c("m1", "m2")) |>
    expect_warning("[Ii]gnored draws of non infected individual .m2.") |>
    suppressWarnings()

  # a model sampling no latent at all has nothing to trace
  mk_latent$model |>
    add_parameter(rho = list(prior = "uniform")) |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  id_draws = "m1") |>
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

test_that("A latent reading a parameter follows each draw in simulation", {
  hh_model <- data.frame(id = seq(40), hh_id = rep(seq(10), each = 4),
                         status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform", max = 20)) |>
    add_latent(hidden = list(model = "gamma", mean = rho, sd = 0.1))

  # each replicate draws around its own value
  simulated <- run_simulation(hh_model, draws = data.frame(rho = c(1, 10)))
  (abs(simulated[[1]]$hidden - 1) < 1) |> expect_all_true()
  (abs(simulated[[2]]$hidden - 10) < 1) |> expect_all_true()
})

test_that("A constant with a model is pinned, its density kept", {
  mk <- function(max) {
    data.frame(id = c("m1", "m2"), hh_id = "h1", status = 1, datum = 7,
               obs = 0) |>
      new_housepi_model("test") |>
      add_latent("hidden", model = "uniform", max = !!max, const = datum) |>
      helper_flat_sink(loc = hidden)
  }

  # simulation reads the pin
  run_simulation(mk(100))$hidden |> expect_all_equal(7)
  # inference reads the density at the pin, one term per member
  (evaluate_log_lik(mk(50)) - evaluate_log_lik(mk(100))) |>
    expect_equal(2 * log(2))
  # the bounds truncate the density
  evaluate_log_lik(mk(5)) |> expect_equal(-Inf)
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

test_that("Bounds of a bare latent read neither draws nor free parameters", {
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
  hh_fit <- mk_onset$model |>
    add_latent(incub_period = list(model = "gamma", shape = 2, scale = 1),
               infec_time = list(min = -14, max = onset_date + 1),
               onset_time = infec_time + incub_period, append = FALSE) |>
    add_observation("onset_date", model = "normal", discr = TRUE,
                    mean = onset_time, sd = 1, min = floor(onset_time)) |>
    run_inference(steps = 10L, thin = 1L, adapt = 0L, chains = 1L,
                  progress = FALSE, warn = FALSE)
  get_log_lik_draws(hh_fit)$h1 |> expect_all_finite()
})

test_that("A bare latent is set by transmission, from the introduction on", {
  # transmission sets every infection, the introduction at the time origin
  # a reader follows the time transmission sets
  simulated <- mk_outbreak$model |>
    add_latent("infec_time", min = -14, max = end_date, append = FALSE) |>
    add_observation(record = infec_time + 1) |>
    run_simulation()
  tapply(simulated$infec_time, simulated$hh_id, min, na.rm = TRUE) |>
    as.vector() |>
    expect_all_equal(0)
  any(simulated$infec_time > 0) |> expect_true()
  infected <- !is.na(simulated$infec_time)
  simulated$record[infected] |> expect_equal(simulated$infec_time[infected] + 1)
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
      run_inference(steps = 100L, thin = 1L, adapt = 0L,
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
  # 20 draws cap the ESS below 100, so the diagnostics warn too
  expect_warning(
    expect_warning(
      hh_fit <- hh_model |>
        run_inference(steps = 20L, chains = 1L, adapt = 50L, thin = 1L),
      "[Cc]ollapsed .*latent .hidden. to constant"
    ),
    "ESS.*below 100 per chain"
  )
  draws <- as.data.frame(get_latent_draws(hh_fit))
  draws$m1_hidden |> expect_all_equal(3)
  length(unique(draws$m2_hidden)) |> expect_gt(1L)
  (draws$m2_hidden >= 3 & draws$m2_hidden <= 5) |> expect_all_true()

  # simulation draws the pinned member at its point too
  expect_warning(
    simulated <- run_simulation(hh_model, nrep = 3L),
    "[Cc]ollapsed .*latent .hidden. to constant"
  )
  vapply(simulated, \(x) x$hidden[1], numeric(1)) |> expect_all_equal(3)

  # both members collapse: one warning for the run, naming the bound sources
  hh_model_2 <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                           datum = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", min = datum, max = datum,
               loc = 2)
  collapsed <- paste0("[Cc]ollapsed .*latent .hidden. to constant",
                      ".*.min.*data .datum.*equals .max.*data .datum")
  expect_warning(simulated <- run_simulation(hh_model_2, nrep = 3L), collapsed)
  unlist(lapply(simulated, `[[`, "hidden")) |> expect_all_equal(3)
  expect_warning(
    log_lik <- hh_model_2 |>
      helper_flat_sink(loc = hidden) |>
      evaluate_log_lik(),
    collapsed
  )
  log_lik |> expect_all_finite()

  # an expression bound meeting the other collapses alike
  data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1, datum = 1,
             obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform", min = datum + 1, max = 2) |>
    helper_flat_sink(loc = hidden) |>
    evaluate_log_lik() |>
    expect_warning("[Cc]ollapsed .*latent .hidden. to constant")
})

test_that("A constant latent with a degenerate model is a point mass", {
  # the bounds meet on 2: the model puts all its mass there, the value stays
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                         record = 2) |>
    new_housepi_model("test") |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1))
  hh_model_2 <- hh_model |>
    add_latent("hidden", const = 3, model = "uniform", min = 2, max = 2)
  hh_model <- hh_model |>
    add_latent("hidden", const = 2, model = "uniform", min = 2, max = 2)
  point_mass <- paste0("[Ll]atent .hidden. has a single possible value",
                       ".*.min. \\(2\\) equals .max. \\(2\\)")

  # a value on the point contributes log 1, one off it is impossible
  expect_warning(log_lik <- evaluate_log_lik(hh_model), point_mass)
  log_lik |> expect_all_finite()
  expect_warning(log_lik_2 <- evaluate_log_lik(hh_model_2), point_mass)
  log_lik_2 |> expect_identical(-Inf)
})

test_that("Nested latents are sampled through their transformation", {
  hh_model <- new_housepi_model(model = "test") |>
    # define Box-Muller transformation
    add_latent(norm = list(model = "normal"),
               u1 = list(model = "uniform"),
               u2 = list(model = "uniform"),
               z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2))
  print(hh_model) |> expect_snapshot()

  # a flat observation on each chain leaves both posteriors equal to the model
  hh_fit <- hh_model |>
    add_dataset(transform(get_dataset(hh_model), obs = 0, obs_norm = 0)) |>
    helper_flat_sink(loc = z0) |>
    helper_flat_sink(loc = norm, column = "obs_norm") |>
    run_inference(steps = 2.5e5L, chains = 1L, adapt = 0L, thin = 1L)
  # no parameters nor priors, only latent
  get_parameter_draws(hh_fit) |> expect_warning("[Ee]mpty draws")
  posterior::nvariables(get_latent_draws(hh_fit)) |>
    expect_identical(nrow(get_latent(hh_model)))

  # the S3 method of summarise_draws dispatches on the fit
  out_vars <- posterior::summarise_draws(hh_fit)[["variable"]]
  out_vars |> expect_match("^tommy")

  # empirical tests for goodness of Gaussian distribution
  probs <- c(0.025, 0.5, 0.975)
  for (draws in list(get_latent_draws(hh_fit)$tommy_norm,
                     get_latent_draws(hh_fit)$tommy_z0)) {
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
  hh_fit <- hh_model |>
    run_inference(steps = 5L, thin = 1L, adapt = 2L,
                  chains = 1L, id_draws = get_dataset(hh_model)$id) |>
    suppressWarnings()
  draws <- as.data.frame(get_latent_draws(hh_fit))
  for (s in seq(4)) {
    trace <- draws[[paste0(s, "_hidden")]]
    trace |> expect_length(5L)
    (trace >= s - 1 & trace <= s) |> expect_all_true()
  }
})

test_that("The latents of a member that declare no proposal move in a block", {
  hh_model <- data.frame(id = c("m1", "m2", "m3"), hh_id = c("h1", "h1", "h2"),
                         status = c(1, 1, 0), obs = 0, obs_2 = 0) |>
    new_housepi_model("test")
  run <- function(...) {
    # each latent reads its own sink, flat over a window it cannot leave
    hh_model <- hh_model |>
      add_latent(..., hidden_2 = list(model = "uniform", min = 0, max = 5)) |>
      helper_flat_sink(loc = hidden) |>
      helper_flat_sink(column = "obs_2", loc = hidden_2)
    hh_model |>
      run_inference(steps = 100L, adapt = 40L, chains = 1L, id_draws = NULL,
                    progress = FALSE, warn = FALSE) |>
      suppressWarnings()
    inspect_state(hh_model, "h1")$.sampler
  }

  # one block per infected member
  sampler <- run(hidden = list(model = "uniform", min = 0, max = 5))
  sampler$moves$member |> expect_setequal(c("m1", "m2"))
  sampler$moves$name |> expect_setequal(c("hidden", "hidden_2"))
  sampler$moves$coord |> expect_all_equal("logit")
  sampler$moves$block |>
    as.vector() |>
    sort() |>
    expect_identical(c(1L, 1L, 2L, 2L))
  sampler$blocks$member |> expect_setequal(c("m1", "m2"))
  vapply(sampler$blocks$jump_cov, nrow, integer(1)) |> expect_all_equal(2L)
  sampler$blocks$jump_cov[[1]] |> expect_dim(c(2L, 2L))
  sampler$blocks$jump_cov[[1]] |>
    rownames() |>
    expect_setequal(c("hidden", "hidden_2"))

  # a declared proposal moves alone
  sampler <- run(hidden = list(model = "uniform", max = 5, proposal = "rw"))
  sampler$moves$block |> expect_all_na()
  sampler$moves$proposal[sampler$moves$name == "hidden"] |>
    expect_all_equal("rw")
  sampler$moves$proposal[sampler$moves$name == "hidden_2"] |> expect_all_na()
  sampler$blocks |> expect_nrow(0L)
})

test_that("Discrete latents join a block beside continuous ones", {
  # record pins hidden + hidden_2 to {11, 12, 13}
  # record_2 pins hidden_3 within 1.5 of 8 - hidden
  hh_model <- data.frame(id = c("m1", "m2"), hh_id = "h1", status = c(1, 0),
                         record = 12, record_2 = 8) |>
    new_housepi_model("test") |>
    add_latent(
      hidden = list(model = "poisson", lambda = 6, discr = TRUE),
      hidden_2 = list(model = "poisson", lambda = 4, discr = TRUE),
      hidden_3 = list(model = "gamma", shape = 3, scale = 1)
    ) |>
    add_observation("record", model = "uniform", min = -1.5, max = 1.5,
                    loc = hidden + hidden_2) |>
    add_observation("record_2", model = "uniform", min = -1.5, max = 1.5,
                    loc = hidden + hidden_3)
  hh_fit <- hh_model |>
    run_inference(steps = 50000L, adapt = 1000L, chains = 1L, id_draws = "m1",
                  progress = FALSE, warn = FALSE) |>
    suppressWarnings()

  # one block: the lattice walks in the identity coordinate
  moves <- inspect_state(hh_model, "h1")$.sampler$moves
  moves$block |> expect_all_equal(1L)
  moves$coord[moves$name != "hidden_3"] |> expect_all_equal("identity")

  # exact posterior over the pairs on the ridge, hidden_3 integrated per pair
  pairs <- expand.grid(hidden = 0:13, hidden_2 = 0:13)
  pairs <- pairs[abs(pairs$hidden + pairs$hidden_2 - 12) <= 1.5, ]
  lower <- pmax(0, 6.5 - pairs$hidden)
  upper <- pmax(0, 9.5 - pairs$hidden)
  prior <- dpois(pairs$hidden, 6) * dpois(pairs$hidden_2, 4)
  weight <- prior * (pgamma(upper, 3) - pgamma(lower, 3))
  weight <- weight / sum(weight)
  # E[hidden_3 1(window)] = 3 P_shape4(window) for a Gamma(3, 1)
  mean_3 <- sum(prior * 3 * (pgamma(upper, 4) - pgamma(lower, 4))) /
    sum(prior * (pgamma(upper, 3) - pgamma(lower, 3)))
  mean_1 <- sum(weight * pairs$hidden)

  draws <- get_latent_draws(hh_fit)
  mean(draws$m1_hidden) |> expect_equal(mean_1, tolerance = 0.03)
  var(draws$m1_hidden) |>
    expect_equal(sum(weight * (pairs$hidden - mean_1)^2), tolerance = 0.2)
  mean(draws$m1_hidden + draws$m1_hidden_2 == 12) |>
    expect_equal(sum(weight[pairs$hidden + pairs$hidden_2 == 12]),
                 tolerance = 0.15)
  mean(draws$m1_hidden_3) |> expect_equal(mean_3, tolerance = 0.07)
})

test_that("A block of discrete latents alone samples their exact posterior", {
  # record pins hidden + hidden_2 to {11, 12, 13}
  # a joint jump rounded to 0 on both is accepted without evaluation
  hh_model <- data.frame(id = c("m1", "m2"), hh_id = "h1", status = c(1, 0),
                         record = 12) |>
    new_housepi_model("test") |>
    add_latent(
      hidden = list(model = "poisson", lambda = 6, discr = TRUE),
      hidden_2 = list(model = "poisson", lambda = 4, discr = TRUE)
    ) |>
    add_observation("record", model = "uniform", min = -1.5, max = 1.5,
                    loc = hidden + hidden_2)
  hh_fit <- hh_model |>
    run_inference(steps = 20000L, adapt = 1000L, chains = 1L, id_draws = "m1",
                  progress = FALSE, warn = FALSE) |>
    suppressWarnings()
  inspect_state(hh_model, "h1")$.sampler$moves$block |> expect_all_equal(1L)

  # exact posterior over the pairs on the ridge
  # limits sit about 7 sd out, measured over 15 runs
  pairs <- expand.grid(hidden = 0:13, hidden_2 = 0:13)
  pairs <- pairs[abs(pairs$hidden + pairs$hidden_2 - 12) <= 1.5, ]
  weight <- dpois(pairs$hidden, 6) * dpois(pairs$hidden_2, 4)
  weight <- weight / sum(weight)
  mean_1 <- sum(weight * pairs$hidden)

  draws <- get_latent_draws(hh_fit)
  mean(draws$m1_hidden) |> expect_equal(mean_1, tolerance = 0.04)
  var(draws$m1_hidden) |>
    expect_equal(sum(weight * (pairs$hidden - mean_1)^2), tolerance = 0.19)
  mean(draws$m1_hidden + draws$m1_hidden_2 == 12) |>
    expect_equal(sum(weight[pairs$hidden + pairs$hidden_2 == 12]),
                 tolerance = 0.15)
})

test_that("Evaluating a model allocates its latent storage", {
  hh_model <- data.frame(hh_id = 1, id = c("m1", "m2"), status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform")
  evaluate_log_lik(hh_model)
  augmented <- inspect_augmented_data(hh_model)
  augmented |> expect_s3_class("data.frame")
  colnames(augmented) |> expect_contains("hidden")
})

####  Observation end-to-end  ####

test_that("A degenerate observation domain is a point mass", {
  # the bounds meet on 2: a continuous model puts all its mass there
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 1,
                         record = 2) |>
    new_housepi_model("test") |>
    add_observation(record = list(model = "normal", mean = 2, sd = 1,
                                  min = 2, max = 2))
  hh_model_2 <- hh_model |>
    add_dataset(transform(get_dataset(hh_model), record = c(2, 3)))
  point_mass <- paste0("[Oo]bservation .record. has a single possible value",
                       ".*.min. \\(2\\) equals .max. \\(2\\)")

  # a datum on the point contributes log 1, one warning for the run
  expect_warning(log_lik <- evaluate_log_lik(hh_model), point_mass)
  log_lik |> expect_equal(0)
  # a datum off the point is impossible
  expect_warning(log_lik_2 <- evaluate_log_lik(hh_model_2), point_mass)
  log_lik_2 |> expect_identical(-Inf)

  # simulation draws the point
  expect_warning(simulated <- run_simulation(hh_model, nrep = 3L), point_mass)
  unlist(lapply(simulated, `[[`, "record")) |> expect_all_equal(2)
})

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
  hh_data <- data.frame(hh_id = "h1", status = 1, record = 5)
  hh_model |>
    add_dataset(hh_data) |>
    run_inference() |>
    expect_error("[Ee]xceeded attempts")

  # sampled parameter in the domain still makes null observation
  hh_data |>
    new_housepi_model("test") |>
    add_parameter("rho", prior = "uniform") |>
    add_observation("record", model = "uniform", min = -rho, max = rho) |>
    run_inference() |>
    expect_error("[Ff]ree parameter dependencies")

  # generic error
  hh_data |>
    new_housepi_model("test") |>
    add_parameter("rho", const = 0.3) |>
    add_latent("hidden", model = "uniform") |>
    add_observation("record", model = "uniform", min = -rho, max = rho) |>
    run_inference() |>
    expect_error("[Oo]utside of respective domain")
})

test_that("An unbounded observation density stops the run", {
  # gamma with shape below 1 is unbounded at 0, reached at hidden = record
  hh_data <- data.frame(hh_id = "h1", status = 1, record = 5)
  mk <- function(...) {
    hh_data |>
      new_housepi_model("test") |>
      add_latent("hidden", ...) |>
      add_observation("record", model = "gamma", shape = 0.5, scale = 1,
                      loc = hidden)
  }

  mk(const = 5) |>
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
      record_3 = list(func = min(record), level = "house")
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
  hh_data <- data.frame(hh_id = rep(c("h1", "h2"), 20), status = 1)
  hh_model <- new_housepi_model(model = "test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 2))

  # observations resolve to continuous in simulations and definitions
  hh_model <- hh_model |>
    add_dataset(dplyr::mutate(hh_data, record = rep(c(2.4, 5.1), 20)))
  get_observations(hh_model)$discr |> expect_false()
  on_integer_lattice(run_simulation(hh_model)$record) |> expect_false()

  # observations resolve to discrete in simulations and definitions
  hh_model <- hh_model |>
    add_dataset(dplyr::mutate(hh_data, record = rep(c(2L, 5L), 20)))
  get_observations(hh_model)$discr |> expect_true()
  on_integer_lattice(run_simulation(hh_model)$record) |> expect_true()
})

test_that("A bare observation name reads its column in a house reduction", {
  simulated <- mk_outbreak$model |>
    add_dataset(dplyr::select(mk_outbreak$data, -test_date)) |>
    add_observation(
      test_date = list(model = "normal", mean = infec_time + incub_period),
      start_date = list(func = min(test_date), level = "house")
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
      infec_time = list(min = start_date - 8, max = end_date)
    ) |>
    add_observation(
      test_date = list(model = "normal", mean = infec_time + incub_period),
      start_date = list(func = min(test_date), level = "house"),
      end_date = start_date + 21,
      start_date_2 = list(func = min(test_date + 1), level = "house")
    ) |>
    run_simulation()

  simulated$end_date |> expect_no_na()
  simulated$end_date |> expect_equal(simulated$start_date + 21)
  # uninfected members leave the rate set at the end of follow-up
  mean(simulated$status > 0) |> expect_lt(1)
  # the argument is computed for each member: an uninfected one draws no test
  # date and drops out of the fold
  as.numeric(simulated$start_date_2) |>
    expect_equal(as.numeric(simulated$start_date) + 1)
})

test_that("An observation model reads an observation drawn before it", {
  # the reader is declared first: the draw order follows the dependency
  simulated <- mk_outbreak$model |>
    add_dataset(dplyr::select(mk_outbreak$data, -test_date)) |>
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
    add_observation("record_2", func = min(record), level = "house") |>
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
  hh_fit <- mk_onset$model |>
    add_observation("onset_date", model = "normal", discr = TRUE, sd = 1,
                    mean = infec_time + incub_period) |>
    run_inference(steps = 50L, thin = 1L, adapt = 0L, chains = 1L) |>
    suppressWarnings()

  get_log_lik_draws(hh_fit)$h1 |> expect_all_finite()
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
  # both the stochastic observation (test_date) and the deterministic one
  # (start_date) are produced internally: neither is provided in the input data
  reps <- mk_outbreak$model |>
    add_dataset(dplyr::select(mk_outbreak$data, -test_date)) |>
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

test_that("A member observation follows the drawn observation it reads", {
  reps <- data.frame(id = seq(40), hh_id = rep(seq(10), each = 4),
                     status = 1) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 2, max = 8)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1),
                    record_2 = record + 1) |>
    run_simulation(nrep = 2)

  # 'record' is drawn per replicate, 'record_2' reads each draw
  for (simulated in reps) {
    (simulated$record_2 - simulated$record) |> expect_equal(rep(1, 40))
  }
})

test_that("A household reduction folds a latent over the members", {
  # the escaped member draws no 'hidden' and drops out of the fold
  hh_model <- data.frame(hh_id = "h1", status = c(1, 1, 0),
                         datum = c(2, 9, 4)) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "normal", mean = 0, sd = 1)) |>
    add_observation(record = list(model = "normal", mean = hidden, sd = 1))

  simulated <- hh_model |>
    add_observation(record_2 = max(hidden), record_3 = sum(hidden + 1),
                    level = "house") |>
    run_simulation()
  unique(simulated$record_2) |>
    expect_equal(max(simulated$hidden, na.rm = TRUE))
  unique(simulated$record_3) |>
    expect_equal(sum(simulated$hidden + 1, na.rm = TRUE))

  # a latent outside any reduction has no household value
  hh_model |>
    add_observation("record_2", func = hidden + 0, level = "house") |>
    run_simulation() |>
    expect_error("[Pp]er-member latent .hidden. is read outside a reduction")
})

test_that("A column needed for every member is not drawn with infection", {
  hh_model <- mk_outbreak$model |>
    add_dataset(dplyr::select(mk_outbreak$data, -end_date)) |>
    add_latent(onset_time = infec_time + incub_period) |>
    add_observation("onset_date", model = "normal", mean = onset_time)

  # onset_date has no value for escaped members, nor has what reads it
  hh_model |>
    add_observation(end_date = onset_date + 21) |>
    run_simulation() |>
    expect_error("[Dd]ata .end_date. needs a value for every member")

  # a household value covers escaped members too
  simulated <- hh_model |>
    add_observation(end_date = min(onset_date) + 21, level = "house") |>
    run_simulation()
  simulated$end_date |> expect_no_na()
})

test_that("A model argument reading an empty cell names member and household", {
  data.frame(hh_id = "h1", status = 1, datum = c(2, NA, 4), record = 0) |>
    new_housepi_model("test") |>
    add_observation(record = list(model = "normal", mean = datum + 1)) |>
    run_simulation() |>
    expect_error(paste0("[Dd]ata .datum. has no value",
                        ".*for\\s+household\\s+.h1.",
                        ".*observation .record. for\\s+individual\\s+.2.",
                        "\\s+\\(status\\s+1\\)"))
})

test_that("A call reading an empty cell is named as written", {
  hh_data <- data.frame(hh_id = "h1", status = 1, datum = c(2, NA, 4),
                        record = 0)
  # `1 - datum` reads datum through its unary minus
  reads <- list("operator\\s+.[+]." = quote(datum + 1),
                "operator\\s+.-." = quote(1 - datum),
                "operator\\s+.[*]." = quote(datum * 2),
                "operator\\s+./." = quote(datum / 2),
                "operator\\s+.\\^." = quote(datum^2),
                "function\\s+.sum." = quote(sum(datum, 1)))
  for (op in names(reads)) {
    hh_data |>
      new_housepi_model("test") |>
      add_observation(record = list(model = "normal", mean = !!reads[[op]])) |>
      run_simulation() |>
      expect_error(paste0("[Dd]ata .datum. has no value.*", op))
  }
})

test_that("Household observations aggregate a data column across members", {
  # one household carries no test_date at all: its na.rm reduction is empty
  empty_house <- mk_outbreak$data$hh_id[[1]]
  hh_data <- dplyr::mutate(
    mk_outbreak$data,
    test_date = dplyr::if_else(hh_id == empty_house, as.Date(NA), test_date)
  )

  simulated <- mk_outbreak$model |>
    add_dataset(hh_data) |>
    add_latent(infec_time = list(min = -14, max = end_date), append = FALSE) |>
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
    values <- tapply(as.numeric(hh_data$test_date), hh_data$hh_id, reduce)
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
  # start_date stands in for provided (static) test_date, read by incub_period
  simulated <- mk_outbreak$model |>
    add_dataset(mk_outbreak$data |>
                  dplyr::mutate(test_date = 5, end_date = 26)) |>
    add_latent(
      incub_period = list(model = "uniform", min = start_date - 1,
                          max = start_date),
      infec_time = list(min = start_date - 5, max = start_date + 5),
      append = FALSE
    ) |>
    add_observation(start_date = min(test_date), level = "house") |>
    add_observation(record = incub_period + 0) |>
    run_simulation()

  # start_date is materialized statically (= 5) before incub_period is drawn
  simulated$start_date |> expect_all_equal(5)
  record <- simulated$record[!is.na(simulated$record)]
  (record >= 4 & record <= 5) |> expect_all_true()
})

test_that("Data-bound household observation supplies a required data column", {
  # the dates of the dataset itself, where end_date is 21 days after the first
  # test date
  hh_data <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 40))
  hh_model <- new_housepi_model(model = "case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(import = "template") |>
    add_latent(import = "template")

  # 'end_date' provided as an ordinary data column
  provided <- add_dataset(hh_model, hh_data)

  # 'end_date' is passed as a deterministic observation
  computed <- hh_model |>
    add_dataset(dplyr::select(hh_data, -end_date)) |>
    add_observation(end_date = min(test_date) + 21, level = "house")

  # the two agree on the pointwise sum: the unpinned template draws the same
  # values under the same seed
  set.seed(42)
  log_lik <- evaluate_log_lik(provided)
  set.seed(42)
  evaluate_log_lik(computed) |>
    expect_equal(log_lik)
})

####  Status end-to-end  ####

test_that("Status process drives status-specific latent blocks", {
  simulated <- mk_outbreak$model |>
    # status-specific block: symptomatic incubation is longer
    add_latent(status = "symptomatic",
               incub_period = list(model = "gamma", shape = 8, scale = 1)) |>
    add_observation(record = incub_period + 0) |>
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
    add_status_process(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  # weak transmission leaves escapees, labelled per the construction argument
  levels(simulated$status) |> expect_contains("healthy")
  as.character(simulated$status) |> expect_contains("healthy")
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
    infec_time = list(const = 0, status = 1),
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14, status = 1),
    onset = list(func = infec_time + incub, status = 1),
    infec_time = list(const = 0, status = 2),
    incub = list(model = "gamma", mean = 3, sd = 1, max = 14, status = 2),
    onset = list(func = infec_time + incub, status = 2)
  )
  suppressWarnings(run_simulation(covered)) |> expect_s3_class("data.frame")

  # a status left without a block is named itself, not reported as a missing
  # generic definition the model never had to declare. Every infected member
  # draws that status, so the report does not depend on how the draw falls
  uncovered <- mk(
    infec_time = list(const = 0, status = 1),
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14, status = 1),
    onset = list(func = infec_time + incub, status = 1),
    weights = c(0, 1)
  )
  cnd <- suppressWarnings(run_simulation(uncovered)) |> expect_error("status 2")
  conditionMessage(cnd) |> expect_no_match("generic latent definition")

  # a generic block on its own still serves every drawn status
  generic <- mk(
    infec_time = 0,
    incub = list(model = "gamma", mean = 2, sd = 1, max = 14),
    onset = list(func = infec_time + incub)
  )
  suppressWarnings(run_simulation(generic)) |> expect_s3_class("data.frame")
})
