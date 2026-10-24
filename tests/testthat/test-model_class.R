library(testthat)

####  Helpers  ####

# transmission model over real data, bare and with its template parameters
mk_frederick <- local({
  hh_model <- new_housepi_model(flu_1918_frederick, "final-size")
  list(model = hh_model,
       template = add_parameter(hh_model, import = "template"))
})

# model carrying the data columns the source lookup refers to
mk_sources <- local({
  hh_model <- data.frame(hh_id = "h1", status = 1, datum = 1,
                         record = 3, datum_2 = 2) |>
    new_housepi_model("test")
  list(model = hh_model)
})

# the dictionary and its assertion are private: reach them through the object
dig_ <- \(hh_model, what, ...) hh_model$.__enclos_env__$private[[what]](...)

####  Standard calls  ####

test_that("Active fields are read-only", {
  new_housepi_model(model = "test") |>
    expect_read_only_field("name", "settings", "escaped")
})

test_that("Escaped label is readable back from the model", {
  new_housepi_model(model = "test")$escaped |>
    expect_identical("escaped")
  new_housepi_model(model = "test", escaped = "healthy")$escaped |>
    expect_identical("healthy")
})

test_that("Model class inconsistent IDs throw errors", {
  # missing individual ID
  data.frame(id = NA, hh_id = "h1") |>
    new_housepi_model("test") |>
    expect_error("[Mm]issing.*id")

  # missing house ID
  data.frame(id = "m1", hh_id = NA) |>
    new_housepi_model("test") |>
    expect_error("[Mm]issing.*hh_id")

  # not unique individual IDs
  data.frame(id = c("m1", "m1"), hh_id = "h1") |>
    new_housepi_model("test") |>
    expect_error("not unique")
})

test_that("Data status is capped by the engine status storage", {
  # above the highest storable infected status
  data.frame(hh_id = "h1", status = 255) |>
    new_housepi_model("test") |>
    expect_error("maximum of 254")

  # a negative status has no meaning: NA declares a missing status
  data.frame(hh_id = "h1", status = -200) |>
    new_housepi_model("test") |>
    expect_error("[Nn]egative status")

  # a missing status round trips as NA, the top level is never exported
  data.frame(hh_id = "h1", status = c(NA, 254)) |>
    new_housepi_model("test") |>
    get_data_wide() |>
    getElement("status") |>
    expect_equal(c(NA, 254))
})

test_that("Household size is capped by the engine member storage", {
  data.frame(hh_id = "h1", status = 0, id = seq(256)) |>
    new_housepi_model("test") |>
    expect_error("[Ee]xceeded the maximum of 255 members")
})

test_that("A missing datum is not a requirement violation", {
  # 'case_ascertained' requires end_date discrete; an absent cell carries no
  # verdict, so it imports and reads back as NA -- the requirement itself is
  # asserted at the run wall
  hh_data <- data.frame(
    id = c("m1", "m2", "m3"), hh_id = "h1", status = c(1, 1, 0),
    end_date = c(21, 21, NA)
  ) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    get_data_wide()

  hh_data$end_date |> expect_equal(c(21, 21, NA))
})

test_that("A datum violating a model requirement is refused at import", {
  mk <- function(...) {
    data.frame(id = c("m1", "m2", "m3"), hh_id = "h1",
               status = c(1, 1, 0), ...) |>
      new_housepi_model("case-ascertained", quiet = TRUE)
  }

  # 'case_ascertained' requires end_date discrete, and the offending row is
  # named alongside the column
  mk(end_date = c(21, 21, 20.5)) |>
    expect_error("end_date.*not discrete-valued")

  # integer-valued cells clear the requirement, negative ones included
  mk(end_date = c(-3, 0, 1)) |>
    get_data_wide() |>
    getElement("end_date") |>
    expect_equal(c(-3, 0, 1))
})

test_that("Model-required data can be postponed to the run wall", {
  # a model-required data column may be supplied later: construction
  # tolerates the gap, the wall (evaluate_*/run_*) names it
  hh_model <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 20)) |>
    dplyr::select(-end_date) |>
    new_housepi_model("case-ascertained")

  hh_model |>
    evaluate_log_lik() |>
    expect_error("end_date")
  hh_model |>
    run_simulation() |>
    expect_error("end_date")
})

test_that("Wide data can be deferred and attached with add_data_wide", {
  # construct without data: no data exposed
  hh_model <- new_housepi_model(model = "final-size") |>
    add_parameter(import = "template")
  get_data_wide(hh_model) |> expect_null()

  # evaluate_/run_ refuse a data-less model
  hh_model |>
    evaluate_log_lik() |>
    expect_error("[Nn]o household data")

  # bad inputs rejected by the wrapper
  hh_model |>
    add_data_wide(NULL) |>
    expect_error("data")
  hh_model |>
    add_data_wide("not a frame") |>
    expect_error("data.frame")

  # attach data later: populated and runnable
  hh_model <- add_data_wide(hh_model, flu_1918_frederick)
  get_data_wide(hh_model) |> expect_s3_class("data.frame")
  get_data_wide(hh_model) |> expect_nrow(nrow(flu_1918_frederick))
  evaluate_log_lik(hh_model) |> expect_all_finite()
})

test_that("Deferred data answers the model requests", {
  # case-ascertained follow-up dates are requested discrete
  hh_model <- new_housepi_model(model = "case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    add_infectivity_profile(import = "template")

  hh_model |>
    add_data_wide(transform(helper_onset_house$data, end_date = 21.5)) |>
    evaluate_log_lik() |>
    expect_error("[Dd]ata.*end_date.*not discrete-valued")
})

test_that("Data without hh_id warns on every attach", {
  # no household structure: every row becomes a household of its own
  hh_data <- data.frame(house_id = "x", id = c("m1", "m2", "m3"))

  new_housepi_model(hh_data, "test") |>
    expect_warning("house_id.*household identifier") |>
    expect_warning("[Ii]gnored household structure")
  # a second model in the same session must warn again
  hh_data$house_id <- NULL
  expect_warning(
    hh_model <- new_housepi_model(hh_data, "test"),
    "[Ii]gnored household structure"
  )
  # and so must a later attach on an existing model
  add_data_wide(hh_model, hh_data) |> expect_warning("[Ii]gnored household")
})

test_that("Model class 'test' calls", {
  hh_model <- new_housepi_model(model = "test")

  # 'test' default data wide format
  get_data_wide(hh_model) |> expect_equal(dummy_household)

  # templates
  add_latent(hh_model, import = "template") |>
    expect_warning("test")

  # zero log-likelihood ('test' model evaluates silently)
  evaluate_log_lik(hh_model) |> expect_equal(0.)

  # factorization round-trips through a fresh 'test' model
  hh_data <- sarscov2_2021_icofs
  hh_data$datum <- as.Date(runif(nrow(hh_data), max = 100))
  get_data_wide(new_housepi_model(hh_data, "test")) |> expect_equal(hh_data)

  # not existing household id
  hh_model |>
    evaluate_log_lik("h9") |>
    expect_error("[Ii]ndex.*h9.*did not match")
})

test_that("Household snapshot reads before any definition is pushed", {
  # entries required by the model are kept even when still empty
  snapshot <- dummy_household |>
    new_housepi_model("final-size") |>
    inspect_state("house")

  field_names <- c(
    "data_wide", "parameters", "group_risks", "contact_risks", "latent",
    "latent_log_lik", "observation_log_lik", "transmission_log_lik"
  )
  names(snapshot) |> expect_identical(field_names)

  for (f in c("parameters", "group_risks", "latent",
              "latent_log_lik", "observation_log_lik")) {
    snapshot[[f]] |> expect_length(0)
  }
  snapshot$transmission_log_lik |> expect_all_na()

  # single household
  snapshot$data_wide |> expect_equal(dummy_household)
})

test_that("Household snapshot reports a household once it is evaluated", {
  # `h2` and `h3` are equivalent: one member escaping alone
  hh_model <- data.frame(id = c("m1", "m2", "m3"), hh_id = c("h1", "h2", "h3"),
                         status = c(1, 0, 0)) |>
    new_housepi_model("final-size") |>
    add_parameter(alpha = 0.3, beta = 0.2, gamma = 1, hh_size_ref = 4)
  expect_unevaluated <- function(hh_model, hh_id) {
    snapshot <- inspect_state(hh_model, hh_id)
    snapshot$transmission_log_lik |> expect_all_na()
    snapshot$model_quantities |> expect_null()
  }
  # the member escapes the community alone
  expect_escaped <- function(hh_model, hh_id) {
    snapshot <- inspect_state(hh_model, hh_id)
    snapshot$transmission_log_lik |> expect_equal(-0.3)
    names(snapshot$model_quantities) |>
      expect_identical(c("community_escape", "p2p_escape"))
  }

  # never evaluated
  expect_unevaluated(hh_model, "h2")

  # one household evaluated, its equivalent with it
  evaluate_log_lik(hh_model, "h2")
  summary <- summarise_house_cases(hh_model)
  summary$count[summary$metric == "classes"] |> expect_equal(2)
  expect_escaped(hh_model, "h2")
  expect_escaped(hh_model, "h3")
  expect_unevaluated(hh_model, "h1")

  # each evaluation starts afresh from the definitions
  evaluate_log_lik(hh_model, "h1")
  inspect_state(hh_model, "h1")$transmission_log_lik |>
    is.finite() |>
    expect_true()
  expect_unevaluated(hh_model, "h2")

  # every household evaluated
  evaluate_log_lik(hh_model)
  expect_escaped(hh_model, "h2")

  # a redefinition or a simulation discards the evaluation
  expect_warning(hh_model_2 <- add_parameter(hh_model, alpha = 0.5),
                 "[Oo]verwriting")
  expect_unevaluated(hh_model_2, "h2")
  run_simulation(hh_model)
  expect_unevaluated(hh_model, "h2")

  # a model exporting nothing drops the entry
  hh_model_3 <- new_housepi_model(dummy_household, "test")
  evaluate_log_lik(hh_model_3)
  snapshot <- inspect_state(hh_model_3, "house")
  snapshot$transmission_log_lik |> is.finite() |> expect_true()
  snapshot$model_quantities |> expect_null()
})

test_that("Data long format is correctly handled", {
  hh_model <- new_housepi_model(model = "test")

  # data_long errors surface via add_data_long
  hh_model |>
    add_data_long("zero") |>
    expect_error("data.frame")

  # every long row names an individual of the wide data
  hh_model |>
    add_data_long(data.frame(id = NA_character_, measure = 1)) |>
    expect_error("[Ii]nvalid id found in row n.1")
  hh_model |>
    add_data_long(data.frame(id = "unknown", measure = 1)) |>
    expect_error("[Ii]d .unknown. from data long-format had no match")
  hh_model |>
    add_data_long(data.frame(measure = 1)) |>
    expect_error("[Ii]dentifier column .id. was missing")
  # long data hangs off the wide data
  new_housepi_model(model = "final-size") |>
    add_data_long(data.frame(id = "m1", measure = 1)) |>
    expect_error("before wide format")

  # correct data long format
  hh_data <- data.frame(id = get_data_wide(hh_model)$id[1],
                        time = seq(3), measure = seq(3)**3)
  hh_model <- add_data_long(hh_model, hh_data)
  get_data_long(hh_model) |> expect_equal(hh_data)

  # inspect_state returns refactored data_long when present
  inspect_state(hh_model, get_data_wide(hh_model)$hh_id[1])$data_long |>
    expect_s3_class("data.frame")
})

test_that("inspect_state reads the parameters and the sampler without a household", {
  hh_model <- data.frame(hh_id = c(1, 1, 2, 2), status = c(1, 0, 1, 1)) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")

  # the sampler reads the default scales, then those adapted by the run
  inspect_state(hh_model)$.sampler$moves$scale |> expect_all_equal(1)
  hh_model |>
    run_inference(steps = 20L, adapt = 10L, chains = 1L, progress = FALSE,
                  warn = FALSE) |>
    suppressWarnings()
  state <- inspect_state(hh_model)
  state |> names() |> expect_identical(c("parameters", ".sampler"))
  state$.sampler$moves$name |> expect_setequal(c("alpha", "beta", "gamma"))
  state$.sampler$moves$block |> expect_all_na()
  (state$.sampler$moves$scale > 0) |> expect_all_true()
  state$.sampler$moves$member |> expect_all_na()
  state$.sampler$blocks |> nrow() |> expect_identical(0L)

  # a household adds its own entries, the sampler last
  inspect_state(hh_model, 1) |> names() |> tail(2L) |>
    expect_identical(c("model_quantities", ".sampler"))

  # a household is one identifier, or NULL
  hh_model |> inspect_state(NA) |> expect_error("single household ID")
  hh_model |> inspect_state(c(1, 2)) |> expect_error("single household ID")
})

test_that("Print behavior works as expected", {
  # 'test' model
  print(new_housepi_model(model = "test")) |>
    expect_snapshot()

  # without data
  print(new_housepi_model(model = "final-size")) |>
    expect_snapshot()

  # without parameters
  print(new_housepi_model(flu_1918_frederick, "final-size")) |>
    expect_snapshot()

  # with data and parameters
  hh_model <- mk_frederick$template
  print(hh_model) |>
    expect_snapshot()

  # help() opens the template page, mocked so no pager is spawned
  called_with <- NULL
  with_mocked_bindings(
    open_model_helper_ = \(name) called_with <<- name,
    code = hh_model$help()
  )
  called_with |> expect_identical("final_size")
})

test_that("Definitions are listed by status, then by name", {
  hh_model <- data.frame(hh_id = "h1", status = c(1, 2), datum = c(3, 4),
                         datum_2 = c(1, 2)) |>
    new_housepi_model("test") |>
    add_parameter(rho_2 = 1, rho = 1) |>
    add_latent(hidden_2 = list(func = min(datum_2), status = 10),
               hidden_2 = list(func = min(datum_2), status = 2),
               hidden = list(func = min(datum), status = 2),
               layer_1 = list(model = "uniform", min = 0, max = 3),
               hidden_3 = list(model = "uniform", min = 0, max = 3))

  # alphabetical by name
  get_parameters(hh_model)[["name"]] |> expect_identical(c("rho", "rho_2"))

  # generic definitions first, then numerically by status, then by name
  latent <- get_latent(hh_model)
  latent[["name"]] |>
    expect_identical(c("hidden_3", "layer_1", "hidden", "hidden_2", "hidden_2"))
  latent[["status"]] |>
    expect_identical(c(NA, NA, "2", "2", "10"))
})

test_that("A copy is deep", {
  house <- data.frame(
    hh_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  )
  hh_model <- house |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  invisible(suppressWarnings(
    run_inference(hh_model, steps = 2L, chains = 1L, adapt = 0L, thin = 1L)
  ))
  hh_model_2 <- hh_model$copy()

  # data
  get_data_wide(hh_model) |> expect_identical(get_data_wide(hh_model_2))
  get_data_long(hh_model) |> expect_identical(get_data_long(hh_model_2))

  # a copy's settings start from the default configuration
  hh_model_2$settings |>
    expect_identical(new_housepi_model(house, "final-size")$settings)

  # definitions
  get_parameters(hh_model) |> expect_identical(get_parameters(hh_model_2))
  get_latent(hh_model) |> expect_identical(get_latent(hh_model_2))
  get_group_risks(hh_model) |> expect_identical(get_group_risks(hh_model_2))

  # changes on the original do not reach the copy
  hh_model <- add_group_risk(hh_model, "susceptibility", by = "member") |>
    suppressWarnings() # suppress baseline collapse warning
  isTRUE(all.equal(get_group_risks(hh_model), get_group_risks(hh_model_2))) |>
    expect_false()
})

####  Source lookup  ####

test_that("Every declared name maps to the source declaring it", {
  # the reserved columns are not addressable as data, the others are
  mk_sources$model |>
    dig_("source_dictionary_") |>
    expect_identical(c(datum = "D", record = "D", datum_2 = "D"))

  dict <- mk_sources$model |>
    add_parameter("rho", prior = "gamma", shape = 2, scale = 1) |>
    add_latent("hidden", model = "gamma", shape = 2, scale = 1) |>
    # a latent shadows the like-named column, as a simulated output fed back
    add_latent("record", model = "gamma", shape = 2, scale = 1) |>
    add_observation("datum_2", model = "normal", mean = hidden, sd = 1) |>
    add_observation("datum_3", level = "house", func = min(datum_2)) |>
    dig_("source_dictionary_")
  unname(dict[c("rho", "hidden", "record", "datum_2", "datum_3")]) |>
    expect_identical(c("P", "A", "A", "D", "D"))

  # two holders under one name: the marker says no read settles it
  dict <- mk_sources$model |>
    add_parameter("datum", prior = "gamma", shape = 2, scale = 1) |>
    dig_("source_dictionary_")
  unname(dict["datum"]) |> expect_identical("")
})

test_that("The summary reports an attached reporting process", {
  mk_frederick$template |>
    add_reporting_process(identity) |>
    print() |>
    expect_message("Reporting process")
})

test_that("An unresolved source is reported softly and at the wall", {
  # nothing unresolved: both modes stay quiet
  dig_(mk_sources$model, "assert_sources_", wall = FALSE) |> expect_length(0)
  dig_(mk_sources$model, "assert_sources_", wall = TRUE) |> expect_length(0)

  hh_model <- mk_sources$model |>
    add_observation(record = list(model = "normal", mean = unknown, sd = 1))
  # soft: the messages are returned, and the model summary carries them
  length(dig_(hh_model, "assert_sources_", wall = FALSE)) |> expect_gt(0)
  print(hh_model) |> expect_message("Sources")
  # wall: the same messages are raised
  dig_(hh_model, "assert_sources_", wall = TRUE) |>
    expect_error("[Cc]ould not resolve the source of .unknown.")
})

####  Data summaries  ####

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hh_model <- new_housepi_model(flu_2014_matsumoto, "final-size")

  case_summary <- summarise_house_cases(hh_model)
  # expect same household count as in Endo, et al. ('19)
  subset(case_summary, level == "household" & metric == "total")$count |>
    expect_identical(10486L)

  # house_size vs status equivalence groups correspond to measure non-zeros
  cases <- measure_cases_by_hh_size(flu_2014_matsumoto)
  cases$hh_size <- NULL
  subset(case_summary, level == "household" & metric == "classes")$count |>
    expect_identical(sum(cases[!is.na(cases)] > 0))

  # fix expected equivalence classes by member composition
  hh_model |>
    add_group_risk("community", by = "member") |>
    summarise_house_cases() |>
    subset(level == "household" & metric == "classes") |>
    getElement("count") |>
    expect_identical(454L)
})

####  Posterior sampling  ####

test_that("Model class inconsistent inference throw errors", {
  # inference without definitions
  new_housepi_model(model = "test") |>
    run_inference() |>
    expect_error("[Nn]othing to sample")

  # missing status: all
  data.frame(id = "m1", hh_id = "h1", status = NA) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference() |>
    expect_error("missing status")
})

test_that("Posterior memory is checked against R's vector limit", {
  # R ignores a limit below its heap capacity, which a collection does not
  # shrink: set one above any heap and ask for draws that cannot fit under it
  limit <- mem.maxVSize()
  withr::defer(mem.maxVSize(limit))
  mem.maxVSize(1e9)
  expect_equal(mem.maxVSize(), 1e9)
  mk_frederick$template |>
    run_inference(steps = 1e15, chains = 1L) |>
    expect_error("[Pp]osterior draws would exceed R.s memory limit")
  # no limit, no check: the engine stops at a model with nothing to sample
  mem.maxVSize(Inf)
  new_housepi_model(model = "test") |>
    run_inference(steps = 1e15, chains = 1e3) |>
    expect_error("[Nn]othing to sample")
})

test_that("Inference rejects invalid settings", {
  mk_frederick$template |>
    run_inference(steps = 5L, thin = 10L) |>
    expect_error("[Tt]hin.*cannot be greater than.*steps")
  mk_frederick$template |>
    run_inference(steps = 5L, thin = 0L) |>
    expect_error("[Tt]hinning interval should be a positive integer")
  mk_frederick$template |>
    run_inference(steps = 5.5) |>
    expect_error("[Nn]umber of steps should be a positive integer")
  mk_frederick$template |>
    run_inference(steps = 5L, adapt = -1L) |>
    expect_error("[Nn]umber of adaptive steps should be a non-negative integer")
  mk_frederick$template |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  threads = 1000L) |>
    expect_warning("[Tt]hreads.*available") |>
    suppressWarnings()
})

test_that("A run stores the sample size its settings ask for", {
  hh_model <- new_housepi_model(model = "test") |>
    add_parameter(rho = list(prior = "uniform", proposal = "forward"))
  n_draws <- function(hh_model, ...) {
    run_inference(hh_model, chains = 1L, ...) |>
      suppressWarnings() |>
      get_parameter_draws() |>
      posterior::niterations()
  }
  # steps are extended to store the last thinned draw: 10 thinned by 3
  n_draws(hh_model, steps = 10L, adapt = 0L, thin = 3L) |> expect_identical(4L)
  # the adaptation budget of a non-adaptive sampler runs as burn-in
  n_draws(hh_model, steps = 5L, adapt = 3L, thin = 1L) |> expect_identical(5L)
  hh_model$settings[c("burn", "adapt")] |>
    expect_equal(list(burn = 3, adapt = 0))
  # an adaptive sampler adapts instead, with no burn-in
  hh_model_2 <- new_housepi_model(model = "test") |>
    add_parameter(rho = list(prior = "uniform"))
  n_draws(hh_model_2, steps = 5L, adapt = 3L, thin = 1L) |> expect_identical(5L)
  hh_model_2$settings[c("burn", "adapt")] |>
    expect_equal(list(burn = 0, adapt = 3))
})

test_that("id_draws selections are resolved before the run", {
  hh_model <- data.frame(id = c("m1", "m2"), hh_id = "h1",
                         status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    helper_flat_sink(loc = hidden)

  hh_model |>
    run_inference(steps = 5L, chains = 1L, id_draws = list(member = "m1")) |>
    expect_error("[Uu]nknown .id_draws. name.*member")

  # NULL stores no latent draw
  hh_model |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, id_draws = NULL) |>
    suppressWarnings() |>
    get_latent_draws() |>
    suppressWarnings() |>
    posterior::nvariables() |>
    expect_identical(0L)

  # NA, the default, follows a household
  hh_model |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, id_draws = NA) |>
    suppressWarnings() |>
    get_latent_draws() |>
    suppressWarnings() |>
    posterior::variables() |>
    expect_setequal(c("m1_hidden", "m2_hidden"))
})

test_that("A followed household without infected members warns", {
  hh_model <- data.frame(id = c("m1", "m2"), hh_id = c("h1", "h2"),
                         status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform")

  hh_model |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L,
                  id_draws = list(hh_id = "h2")) |>
    expect_warning("no infected members") |>
    suppressWarnings()
})

####  Simulation run  ####

test_that("Simulation rejects wrong arguments", {
  hh_model <- new_housepi_model(model = "test")
  # wrong number of replicates
  hh_model |>
    run_simulation(nrep = 0) |>
    expect_error("positive integer")
  # the integration step is a positive duration
  hh_model |>
    run_simulation(time_step = 0) |>
    expect_error("time_step. must be strictly positive")
  run_simulation(mk_frederick$template, time_step = 0.5) |>
    expect_s3_class("data.frame")
})

test_that("Param table with degenerate columns raises an ESS error", {
  new_housepi_model(model = "test") |>
    run_simulation(draws = posterior::as_draws_df(
      data.frame(rho = rep(1, 100))
    )) |>
    expect_error("ESS")
})

test_that("Thinning to ESS keeps from a few draws up to all draws", {
  hh_model <- data.frame(hh_id = "h1", id = c("m1", "m2"), status = 0) |>
    new_housepi_model("test")
  # two disjoint trending chains of 20 iterations: ESS near the chain count
  trend <- posterior::as_draws_df(data.frame(
    rho = rep(c(10, 20), each = 20) + rep(seq_len(20L), 2) / 20,
    .chain = rep(1:2, each = 20), .iteration = rep(seq_len(20L), 2)
  ))
  thinned <- run_simulation(hh_model, draws = trend)
  length(thinned) |> expect_gte(2L)
  length(thinned) |> expect_lt(5L)
  # without thinning, every draw is simulated, under a cap above their count
  hh_model |>
    run_simulation(draws = trend, thin_ess = FALSE, nrep = 100) |>
    expect_length(40L)
  # one chain of 20 iterations alternating low and high draws:
  # its ESS lies above the draw count, where `posterior` caps it
  alternate <- data.frame(rho = c(rbind(1:10, 20:11)))
  hh_model |>
    run_simulation(draws = posterior::as_draws_df(alternate)) |>
    expect_length(20L)
  # one trending chain: its ESS floors to a single draw, still in a list
  hh_model |>
    run_simulation(draws = posterior::as_draws_df(
      data.frame(rho = seq_len(100L))
    )) |>
    expect_length(1L)
})

test_that("Simulation reads undeclared parameters from draws", {
  hh_model <- mk_frederick$model$copy()
  param_table <- data.frame(alpha = c(0.3, 0.1), beta = c(0.2, 0.4),
                            gamma = 1, hh_size_ref = 4)
  hh_model |>
    run_simulation(draws = param_table) |>
    expect_length(2L)
  # the model holds the parameters of the last replicate
  state <- hh_model |> inspect_state(get_data_wide(hh_model)$hh_id[1])
  unlist(state$parameters[names(param_table)]) |>
    expect_equal(unlist(param_table[2, ]), ignore_attr = TRUE)
  # the table is not declared on the model
  get_parameters(hh_model) |> expect_nrow(0L)
  hh_model |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing parameter definition")

  # a parameter defined as a function takes no value from the table
  data.frame(hh_id = "h1", status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform"), rho_2 = exp(rho)) |>
    run_simulation(draws = data.frame(rho = 0.5, rho_2 = 3)) |>
    expect_error("[Nn]on numeric parameter .rho_2.")
})

test_that("Simulation returns combined wide/long output", {
  # a long reporting process forces the long branch into a non-empty list
  hh_model <- mk_frederick$template |>
    add_reporting_process(long = \(x) data.frame(id = "m1", measure = 1.0))
  simulated <- hh_model |> run_simulation()
  simulated |> expect_named(c("wide", "long"))
  simulated$wide |> expect_s3_class("data.frame")
  simulated$long |> expect_s3_class("data.frame")

  # multiple repetitions
  simulated <- hh_model |> run_simulation(nrep = 2)
  simulated$wide[[1]] |> expect_s3_class("data.frame")
  simulated$long[[2]] |> expect_s3_class("data.frame")
})

test_that("Each simulated repetition holds its own draws", {
  # a long frame passing a column through keeps its repetition's draws
  simulated <- data.frame(id = c("m1", "m2", "m3"), hh_id = "h1",
                          status = c(1, 1, 0), record = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record = hidden + 2) |>
    add_reporting_process(long = \(x) x[, c("id", "hidden")]) |>
    run_simulation(nrep = 2)
  simulated$long[[1]]$hidden |> expect_identical(simulated$wide[[1]]$hidden)
  simulated$long[[2]]$hidden |> expect_identical(simulated$wide[[2]]$hidden)
})

test_that("Simulation replicates are a list of data frames", {
  # even a single replicate, once `nrep` asks for replicates
  simulated <- run_simulation(mk_frederick$template, nrep = 1)
  simulated |> expect_type("list")
  simulated |> expect_length(1L)
  simulated[[1]] |> expect_s3_class("data.frame")

  # more draws than nrep are thinned down to nrep
  param_table <- data.frame(alpha = rep(0.15, 10L), beta = 0.2, gamma = 1)
  run_simulation(mk_frederick$template, draws = param_table, nrep = 4) |>
    expect_length(4L)
})

test_that("A simulation with a reporting process is reproducible", {
  hh_model <- mk_frederick$template |>
    add_reporting_process(
      \(x) dplyr::mutate(dplyr::rowwise(x), flag = runif(1))
    )
  simulated <- withr::with_seed(1, run_simulation(hh_model))
  ("flag" %in% colnames(simulated)) |> expect_true()
  withr::with_seed(1, run_simulation(hh_model)) |>
    expect_identical(simulated)
})

test_that("A household error in simulation names the household", {
  # a negative datum in h2 empties the latent domain
  hh_model <- data.frame(hh_id = c("h1", "h1", "h2"), id = c("m1", "m2", "m3"),
                         status = 1, datum = c(1, 1, -1), obs = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 0, max = datum)) |>
    helper_flat_sink(loc = hidden)
  err <- run_simulation(hh_model) |> expect_error("[Ii]nconsistent domain")
  # cli markup is rendered, not printed raw
  conditionMessage(err) |> expect_match("household \"h2\"")
  conditionMessage(err) |> expect_no_match("{.val", fixed = TRUE)
})

test_that("An error raised by a reporting process reaches R verbatim", {
  mk_frederick$template |>
    add_reporting_process(\(x) stop("a {b} c")) |>
    run_simulation() |>
    expect_error("{b} c", fixed = TRUE)
})

test_that("A failed simulation leaves the observed data in place", {
  hh_model <- mk_frederick$model |>
    add_parameter(alpha = 0.15, beta = 0.2, gamma = 1, hh_size_ref = 4)
  log_lik <- evaluate_log_lik(hh_model)
  # the second replicate fails after the first one is simulated
  calls <- 0L
  fail_second <- \(x) if ((calls <<- calls + 1L) == 2L) stop("stop") else x
  hh_model |>
    add_reporting_process(fail_second) |>
    run_simulation(nrep = 3L) |>
    expect_error("[Ss]top")
  hh_model |>
    evaluate_log_lik() |>
    expect_identical(log_lik)
})
