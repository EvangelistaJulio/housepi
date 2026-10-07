library(testthat)

# transmission model over real data, bare and with its template parameters
mk_frederick <- local({
  hhmodel <- new_housepi_model(flu_1918_frederick, "final-size")
  list(model = hhmodel, template = add_parameter(hhmodel, import = "template"))
})

# model carrying the data columns the source lookup refers to
mk_sources <- local({
  hhmodel <- data.frame(hhid = "h1", status = 1, datum = 1,
                        record = 3, datum_2 = 2) |>
    new_housepi_model("test")
  list(model = hhmodel)
})

# the dictionary and its assertion are private: reach them through the object
dig_ <- \(hhmodel, what, ...) hhmodel$.__enclos_env__$private[[what]](...)

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Active fields are read-only", {
  new_housepi_model(model = "test") |>
    expect_read_only_field(
      "data_wide", "data_long", "parameters", "latent", "observations",
      "infectivity_profiles", "group_risks", "contact_risks", "name", "setup",
      "escaped"
    )
})

test_that("Escaped label is readable back from the model", {
  new_housepi_model(model = "test")$escaped |>
    expect_identical("escaped")
  new_housepi_model(model = "test", escaped = "healthy")$escaped |>
    expect_identical("healthy")
})

test_that("Model class inconsistent IDs throw errors", {
  # missing individual ID
  data.frame(id = NA, hhid = "h1") |>
    new_housepi_model("test") |>
    expect_error("[Mm]issing.*id")

  # missing house ID
  data.frame(id = "m1", hhid = NA) |>
    new_housepi_model("test") |>
    expect_error("[Mm]issing.*hhid")

  # not unique individual IDs
  data.frame(id = c("m1", "m1"), hhid = "h1") |>
    new_housepi_model("test") |>
    expect_error("not unique")
})

test_that("Data status is capped by the engine status storage", {
  # above the highest storable infected status
  data.frame(hhid = "h1", status = 255) |>
    new_housepi_model("test") |>
    expect_error("maximum of 254")

  # a negative status has no meaning: NA declares a missing status
  data.frame(hhid = "h1", status = -200) |>
    new_housepi_model("test") |>
    expect_error("[Nn]egative status")

  # a missing status round trips as NA, the top level is never exported
  data.frame(hhid = "h1", status = c(NA, 254)) |>
    new_housepi_model("test") |>
    get_data_wide() |>
    getElement("status") |>
    expect_equal(c(NA, 254))
})

test_that("Household size is capped by the engine member storage", {
  data.frame(hhid = "h1", status = 0, id = seq(256)) |>
    new_housepi_model("test") |>
    expect_error("[Ee]xceeded the maximum of 255 members")
})

test_that("A missing datum is not a requirement violation", {
  # 'case_ascertained' requires start_date and end_date discrete; an absent
  # cell carries no verdict, so it imports and reads back as NA -- the
  # requirement itself is asserted at the run wall
  hhdata <- data.frame(
    id = c("m1", "m2", "m3"), hhid = "h1", status = c(1, 1, 0),
    start_date = c(1, NA, 1), end_date = c(21, 21, NA)
  ) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    get_data_wide()

  hhdata$start_date |> expect_equal(c(1, NA, 1))
  hhdata$end_date |> expect_equal(c(21, 21, NA))
})

test_that("A datum violating a model requirement is refused at import", {
  mk <- function(...) {
    data.frame(id = c("m1", "m2", "m3"), hhid = "h1",
               status = c(1, 1, 0), ...) |>
      new_housepi_model("case-ascertained", quiet = TRUE)
  }

  # 'case_ascertained' requires start_date and end_date discrete, and the
  # offending row is named alongside the column
  mk(start_date = c(1, 2.5, 1), end_date = 21) |>
    expect_error("start_date.*not discrete-valued")

  mk(start_date = 1, end_date = c(21, 21, 20.5)) |>
    expect_error("not discrete-valued")

  # integer-valued cells clear the requirement, negative ones included
  mk(start_date = c(-3, 0, 1), end_date = 21) |>
    get_data_wide() |>
    getElement("start_date") |>
    expect_equal(c(-3, 0, 1))
})

test_that("Model-required data can be postponed to the run wall", {
  # a model-required data column may be supplied later: construction
  # tolerates the gap, the wall (evaluate_*/run_*) names it
  hhmodel <- sarscov2_2021_icofs |>
    dplyr::filter(hhid %in% utils::head(unique(hhid), 20)) |>
    dplyr::select(-start_date) |>
    new_housepi_model("case-ascertained")

  hhmodel |>
    evaluate_log_lik() |>
    expect_error("start_date")
  hhmodel |>
    run_simulation() |>
    expect_error("start_date")
})

test_that("Wide data can be deferred and attached with add_data_wide", {
  # construct without data: no data exposed
  hhmodel <- new_housepi_model(model = "final-size") |>
    add_parameter(import = "template")
  hhmodel$data_wide |> expect_null()

  # evaluate_/run_ refuse a data-less model
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("[Nn]o household data")

  # bad inputs rejected by the wrapper
  hhmodel |>
    add_data_wide(NULL) |>
    expect_error("data")
  hhmodel |>
    add_data_wide("not a frame") |>
    expect_error("data.frame")

  # attach data later: populated and runnable
  hhmodel <- add_data_wide(hhmodel, flu_1918_frederick)
  hhmodel$data_wide |> expect_s3_class("data.frame")
  hhmodel$data_wide |> expect_nrow(nrow(flu_1918_frederick))
  is.finite(evaluate_log_lik(hhmodel)) |> expect_true()
})

test_that("Deferred data answers the model requests", {
  # case-ascertained follow-up dates are requested discrete
  hhmodel <- new_housepi_model(model = "case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    add_infectivity_profile(import = "template")

  hhmodel |>
    add_data_wide(transform(helper_onset_house$data, end_date = 21.5)) |>
    evaluate_log_lik() |>
    expect_error("[Dd]ata.*end_date.*not discrete-valued")
})

test_that("Data without hhid warns on every attach", {
  # no household structure: every row becomes a household of its own
  hhdata <- data.frame(id = c("m1", "m2", "m3"), status = c(1, 0, 1))

  new_housepi_model(hhdata, "test") |>
    expect_warning("[Ii]gnored household structure")
  # a second model in the same session must warn again
  expect_warning(
    hhmodel <- new_housepi_model(hhdata, "test"),
    "[Ii]gnored household structure"
  )
  # and so must a later attach on an existing model
  add_data_wide(hhmodel, hhdata) |> expect_warning("[Ii]gnored household")
})

test_that("Model class 'test' calls", {
  hhmodel <- new_housepi_model(model = "test")

  # 'test' default data wide format
  hhmodel$data_wide |> expect_equal(dummy_household)

  # templates
  add_latent(hhmodel, import = "template") |>
    expect_warning("test")

  # zero log-likelihood ('test' model evaluates silently)
  evaluate_log_lik(hhmodel) |> expect_equal(0.)

  # factorization round-trips through a fresh 'test' model
  hhdata <- sarscov2_2021_icofs
  hhdata$datum <- as.Date(runif(nrow(hhdata), max = 100))
  new_housepi_model(hhdata, "test")$data_wide |> expect_equal(hhdata)

  # not existing household id
  hhmodel |>
    evaluate_log_lik("h9") |>
    expect_error("[Ii]ndex.*h9.*did not match")
})

test_that("Household snapshot reads before any definition is pushed", {
  # entries required by the model are kept even when still empty
  snapshot <- dummy_household |>
    new_housepi_model("final-size") |>
    inspect_house_state("house")

  field_names <- c(
    "data_wide", "parameters", "group_risks", "contact_risks", "latent",
    "latent_log_lik", "observation_log_lik", "transmission_log_lik"
  )
  names(snapshot) |> expect_identical(field_names)

  for (f in c("parameters", "group_risks", "latent",
              "latent_log_lik", "observation_log_lik")) {
    snapshot[[f]] |> expect_length(0)
  }

  # single household
  snapshot$data_wide |> expect_equal(dummy_household)
})

test_that("Data long format is correctly handled", {
  hhmodel <- new_housepi_model(model = "test")

  # data_long errors surface via add_data_long
  hhmodel |>
    add_data_long("zero") |>
    expect_error("data.frame")

  # every long row names an individual of the wide data
  hhmodel |>
    add_data_long(data.frame(id = NA_character_, measure = 1)) |>
    expect_error("[Ii]nvalid id found in row n.1")
  hhmodel |>
    add_data_long(data.frame(id = "unknown", measure = 1)) |>
    expect_error("[Ii]d .unknown. from data long-format had no match")
  hhmodel |>
    add_data_long(data.frame(measure = 1)) |>
    expect_error("[Ii]dentifier column .id. was missing")
  # long data hangs off the wide data
  new_housepi_model(model = "final-size") |>
    add_data_long(data.frame(id = "m1", measure = 1)) |>
    expect_error("before wide format")

  # correct data long format
  hhdata <- data.frame(id = hhmodel$data_wide$id[1],
                       time = seq(3), measure = seq(3)**3)
  hhmodel <- add_data_long(hhmodel, hhdata)
  get_data_long(hhmodel) |> expect_equal(hhdata)

  # inspect_house_state returns refactored data_long when present
  inspect_house_state(hhmodel, hhmodel$data_wide$hhid[1])$data_long |>
    expect_s3_class("data.frame")
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
  hhmodel <- mk_frederick$template
  print(hhmodel) |>
    expect_snapshot()

  # help() opens the template page, mocked so no pager is spawned
  called_with <- NULL
  with_mocked_bindings(
    open_model_helper_ = \(name) called_with <<- name,
    code = hhmodel$help()
  )
  called_with |> expect_identical("final_size")
})

test_that("A copy is deep", {
  house <- data.frame(
    hhid = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  )
  hhmodel <- house |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  invisible(suppressWarnings(
    run_inference(hhmodel, steps = 2L, chains = 1L, adapt = 0L, thin = 1L)
  ))
  hhmodel_2 <- hhmodel$copy()

  # data
  hhmodel$data_wide |> expect_identical(hhmodel_2$data_wide)
  hhmodel$data_long |> expect_identical(hhmodel_2$data_long)

  # a copy's setup starts from the default configuration
  hhmodel_2$setup |>
    expect_identical(new_housepi_model(house, "final-size")$setup)

  # definitions
  hhmodel$parameters |> expect_identical(hhmodel_2$parameters)
  hhmodel$latent |> expect_identical(hhmodel_2$latent)
  hhmodel$group_risks |> expect_identical(hhmodel_2$group_risks)

  # changes on the original do not reach the copy
  hhmodel <- add_group_risk(hhmodel, "susceptibility", by = "member") |>
    suppressWarnings() # suppress baseline collapse warning
  isTRUE(all.equal(hhmodel$group_risks, hhmodel_2$group_risks)) |>
    expect_false()
})

#################################
## MODEL CLASS : SOURCE LOOKUP ##
#################################

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
    add_observation("datum_3", level = "house", value = min(datum_2)) |>
    dig_("source_dictionary_")
  unname(dict[c("rho", "hidden", "record", "datum_2", "datum_3")]) |>
    expect_identical(c("P", "A", "A", "D", "D"))

  # two holders under one name: the marker says no read settles it
  dict <- mk_sources$model |>
    add_parameter("datum", prior = "gamma", shape = 2, scale = 1) |>
    dig_("source_dictionary_")
  unname(dict["datum"]) |> expect_identical("")
})

test_that("An unresolved source is reported softly and at the wall", {
  # nothing unresolved: both modes stay quiet
  dig_(mk_sources$model, "assert_sources_", wall = FALSE) |> expect_length(0)
  dig_(mk_sources$model, "assert_sources_", wall = TRUE) |> expect_length(0)

  hhmodel <- mk_sources$model |>
    add_observation(record = list(model = "normal", mean = unknown, sd = 1))
  # soft: the messages are returned, and the model summary carries them
  length(dig_(hhmodel, "assert_sources_", wall = FALSE)) |> expect_gt(0)
  print(hhmodel) |> expect_output("Sources")
  # wall: the same messages are raised
  dig_(hhmodel, "assert_sources_", wall = TRUE) |>
    expect_error("[Cc]ould not resolve the source of .unknown.")
})

##################################
## MODEL CLASS : DATA SUMMARIES ##
##################################

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hhmodel <- new_housepi_model(flu_2014_matsumoto, "final-size")

  case_summary <- summarise_house_cases(hhmodel)
  # expect same household count as in Endo, et al. ('19)
  subset(case_summary, level == "household" & metric == "total")$count |>
    expect_identical(10486L)

  # house_size vs status equivalence groups correspond to measure non-zeros
  cases <- measure_cases_by_hhsize(flu_2014_matsumoto)
  cases$hhsize <- NULL
  subset(case_summary, level == "household" & metric == "classes")$count |>
    expect_identical(sum(cases[!is.na(cases)] > 0))

  # fix expected equivalence classes by member composition
  hhmodel |>
    add_group_risk("community", by = "member") |>
    summarise_house_cases() |>
    subset(level == "household" & metric == "classes") |>
    getElement("count") |>
    expect_identical(454L)
})

######################################
## MODEL CLASS : POSTERIOR SAMPLING ##
######################################

test_that("Model class inconsistent inference throw errors", {
  # inference without definitions
  new_housepi_model(model = "test") |>
    run_inference() |>
    expect_error("[Nn]othing to sample")

  # missing status: all
  data.frame(id = "m1", hhid = "h1", status = NA) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference() |>
    expect_error("missing status")
})

test_that("Posterior memory budget is asserted", {
  data.frame(hhid = 1, id = "m1", status = 0) |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference(steps = 1e15, chains = 1e3) |>
    expect_error("memory budget")
})

test_that("Inference rejects invalid setup arguments", {
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
    run_inference(steps = 5L, burn = -1L) |>
    expect_error("[Nn]umber of burn-in steps should be a non-negative integer")
  mk_frederick$template |>
    run_inference(steps = 5L, adapt = -1L) |>
    expect_error("[Nn]umber of adaptive steps should be a non-negative integer")
  mk_frederick$template |>
    run_inference(steps = 5L, adapt_window = 0L) |>
    expect_error("[Aa]daptation window should be a positive integer")
  mk_frederick$template |>
    run_inference(steps = 5L, adapt_windows = 10L) |>
    expect_error("[Uu]nexpected argument.*adapt_windows")
  mk_frederick$template |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  threads = 1000L) |>
    expect_warning("[Tt]hreads.*available") |>
    suppressWarnings()
})

test_that("A run stores the sample size its setup asks for", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(rho = list(prior = "uniform", proposal = "forward"))
  n_draws <- function(...) {
    run_inference(hhmodel, chains = 1L, ...) |>
      suppressWarnings() |>
      get_parameter_draws() |>
      posterior::niterations()
  }
  # burn-in runs on top of steps
  n_draws(steps = 5L, burn = 10L, adapt = 0L, thin = 1L) |> expect_identical(5L)
  # steps are extended to store the last thinned draw: 10 thinned by 3
  n_draws(steps = 10L, burn = 0L, adapt = 0L, thin = 3L) |> expect_identical(4L)
  # the adaptation budget of a non-adaptive sampler runs as burn-in
  n_draws(steps = 5L, burn = 0L, adapt = 3L, thin = 1L) |> expect_identical(5L)
  hhmodel$setup[c("burn", "adapt")] |> expect_equal(list(burn = 3, adapt = 0))
  # unless a burn-in is set: then it is dropped
  n_draws(steps = 5L, burn = 2L, adapt = 3L, thin = 1L) |> expect_identical(5L)
  hhmodel$setup[c("burn", "adapt")] |> expect_equal(list(burn = 2, adapt = 0))
})

test_that("id_draws selections are resolved before the run", {
  hhmodel <- data.frame(id = c("m1", "m2"), hhid = "h1",
                        status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    helper_flat_sink(loc = hidden)

  hhmodel |>
    run_inference(steps = 5L, chains = 1L, id_draws = list(member = "m1")) |>
    expect_error("[Uu]nknown .id_draws. name.*member")

  # the explicit opt-out stores no latent draw
  hhmodel |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, id_draws = ".none") |>
    suppressWarnings() |>
    get_latent_draws() |>
    suppressWarnings() |>
    posterior::nvariables() |>
    expect_identical(0L)
})

test_that("A followed household without infected members warns", {
  hhmodel <- data.frame(id = c("m1", "m2"), hhid = c("h1", "h2"),
                        status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform")

  hhmodel |>
    run_inference(steps = 20L, adapt = 0L, chains = 1L,
                  id_draws = list(hhid = "h2")) |>
    expect_warning("no infected members") |>
    suppressWarnings()
})

##################################
## MODEL CLASS : SIMULATION RUN ##
##################################

test_that("Simulation rejects wrong arguments", {
  hhmodel <- new_housepi_model(model = "test")
  # wrong observation process
  hhmodel |>
    run_simulation(obs_process = "not a function") |>
    expect_error("function")
  # wrong number of replicates
  hhmodel |>
    run_simulation(nrep = 0) |>
    expect_error("positive integer")
  # the integration step is a positive duration
  hhmodel |>
    run_simulation(time_step = 0) |>
    expect_error("time_step. must be strictly positive")
  run_simulation(mk_frederick$template, time_step = 0.5) |>
    expect_s3_class("data.frame")
})

test_that("Param table with degenerate columns raises an ESS error", {
  new_housepi_model(model = "test") |>
    run_simulation(param_table = posterior::as_draws_df(
      data.frame(rho = rep(1, 100))
    )) |>
    expect_error("ESS")
})

test_that("Thinning to ESS keeps from a few draws up to all draws", {
  hhmodel <- data.frame(hhid = "h1", id = c("m1", "m2"), status = 0) |>
    new_housepi_model("test")
  # two disjoint trending chains of 20 iterations: ESS near the chain count
  trend <- posterior::as_draws_df(data.frame(
    rho = rep(c(10, 20), each = 20) + rep(seq_len(20L), 2) / 20,
    .chain = rep(1:2, each = 20), .iteration = rep(seq_len(20L), 2)
  ))
  thinned <- run_simulation(hhmodel, param_table = trend)
  length(thinned) |> expect_gte(2L)
  length(thinned) |> expect_lt(5L)
  # without thinning, every draw is simulated, under a cap above their count
  hhmodel |>
    run_simulation(param_table = trend, thin_ess = FALSE, nrep = 100) |>
    expect_length(40L)
  # one chain of 20 iterations alternating low and high draws:
  # its ESS lies above the draw count, where `posterior` caps it
  alternate <- data.frame(rho = c(rbind(1:10, 20:11)))
  hhmodel |>
    run_simulation(param_table = posterior::as_draws_df(alternate)) |>
    expect_length(20L)
  # one trending chain: its ESS floors to a single draw, still in a list
  hhmodel |>
    run_simulation(param_table = posterior::as_draws_df(
      data.frame(rho = seq_len(100L))
    )) |>
    expect_length(1L)
})

test_that("Simulation reads undeclared parameters from param_table", {
  hhmodel <- mk_frederick$model$copy()
  param_table <- data.frame(alpha = c(0.3, 0.1), beta = c(0.2, 0.4),
                            gamma = 1, hhsize_ref = 4)
  hhmodel |>
    run_simulation(param_table = param_table) |>
    expect_length(2L)
  # the model holds the parameters of the last replicate
  state <- hhmodel |> inspect_house_state(hhmodel$data_wide$hhid[1])
  unlist(state$parameters[names(param_table)]) |>
    expect_equal(unlist(param_table[2, ]), ignore_attr = TRUE)
  # the table is not declared on the model
  get_parameters(hhmodel) |> expect_nrow(0L)
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing parameter definition")

  # a parameter defined as a function takes no value from the table
  data.frame(hhid = "h1", status = 1) |>
    new_housepi_model("test") |>
    add_parameter(rho = list(prior = "uniform"), rho_2 = exp(rho)) |>
    run_simulation(param_table = data.frame(rho = 0.5, rho_2 = 3)) |>
    expect_error("[Nn]on numeric parameter .rho_2.")
})

test_that("Simulation returns combined wide/long output", {
  # supplying obs_process_long forces the long branch into a non-empty list
  obs_long <- \(x) data.frame(id = "m1", measure = 1.0)
  simulated <- mk_frederick$template |>
    run_simulation(obs_process_long = obs_long)
  simulated |> expect_named(c("wide", "long"))
  simulated$wide |> expect_s3_class("data.frame")
  simulated$long |> expect_s3_class("data.frame")

  # multiple repetitions
  simulated <- mk_frederick$template |>
    run_simulation(obs_process_long = obs_long, nrep = 2)
  simulated$wide[[1]] |> expect_s3_class("data.frame")
  simulated$long[[2]] |> expect_s3_class("data.frame")
})

test_that("Each simulated repetition holds its own draws", {
  # a long frame passing a column through keeps its repetition's draws
  results <- data.frame(id = c("m1", "m2", "m3"), hhid = "h1",
                        status = c(1, 1, 0), record = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform")) |>
    add_observation(record = hidden + 2) |>
    run_simulation(nrep = 2, obs_process_long = \(x) x[, c("id", "hidden")])
  results$long[[1]]$hidden |> expect_identical(results$wide[[1]]$hidden)
  results$long[[2]]$hidden |> expect_identical(results$wide[[2]]$hidden)
})

test_that("Simulation replicates are a list of data frames", {
  # even a single replicate, once `nrep` asks for replicates
  simulated <- run_simulation(mk_frederick$template, nrep = 1)
  simulated |> expect_type("list")
  simulated |> expect_length(1L)
  simulated[[1]] |> expect_s3_class("data.frame")

  # a param_table larger than nrep is thinned down to it
  params <- get_parameters(mk_frederick$template)
  param_table <- as.data.frame(as.list(stats::setNames(params$value,
                                                       params$name)))
  param_table <- param_table[rep(1L, 10L), , drop = FALSE]
  run_simulation(mk_frederick$template, param_table = param_table, nrep = 4) |>
    expect_length(4L)
})

test_that("A simulation with an observation process is reproducible", {
  flag <- \(x) dplyr::mutate(dplyr::rowwise(x), flag = runif(1))
  simulated <- withr::with_seed(
    1, run_simulation(mk_frederick$template, obs_process = flag)
  )
  ("flag" %in% colnames(simulated)) |> expect_true()
  withr::with_seed(
    1, run_simulation(mk_frederick$template, obs_process = flag)
  ) |>
    expect_identical(simulated)
})

test_that("Measures on long data require obs_process_long", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(rho = 1)

  dig_(hhmodel, "run_measure_",
    param_table = data.frame(rho = c(0.1, 0.2, 0.3)),
    apply_wide = list,
    apply_long = list,
    obs_proc_wide = NULL,
    obs_proc_long = NULL
  ) |>
    expect_error("obs_process_long")
})

test_that("A household error in simulation names the household", {
  # a negative datum in h2 empties the latent domain
  hhmodel <- data.frame(hhid = c("h1", "h1", "h2"), id = c("m1", "m2", "m3"),
                        status = 1, datum = c(1, 1, -1), obs = 0) |>
    new_housepi_model("test") |>
    add_latent(hidden = list(model = "uniform", min = 0, max = datum)) |>
    helper_flat_sink(loc = hidden)
  err <- run_simulation(hhmodel) |> expect_error("[Ii]nconsistent domain")
  # cli markup is rendered, not printed raw
  conditionMessage(err) |> expect_match("household \"h2\"")
  conditionMessage(err) |> expect_no_match("{.val", fixed = TRUE)
})

test_that("An error raised by an observation process reaches R verbatim", {
  mk_frederick$template |>
    run_simulation(obs_process = \(x) stop("a {b} c")) |>
    expect_error("{b} c", fixed = TRUE)
})

test_that("A failed simulation leaves the observed data in place", {
  hhmodel <- mk_frederick$model |>
    add_parameter(alpha = 0.15, beta = 0.2, gamma = 1, hhsize_ref = 4)
  log_lik <- evaluate_log_lik(hhmodel)
  # the second replicate fails after the first one is simulated
  calls <- 0L
  fail_second <- \(x) if ((calls <<- calls + 1L) == 2L) stop("stop") else x
  hhmodel |>
    run_simulation(nrep = 3L, obs_process = fail_second) |>
    expect_error("[Ss]top")
  hhmodel |>
    evaluate_log_lik() |>
    expect_identical(log_lik)
})
