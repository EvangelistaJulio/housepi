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
dig_ <- function(hhmodel, what, ...) {
  hhmodel$.__enclos_env__$private[[what]](...)
}

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Active fields are read-only", {
  new_housepi_model(model = "test") |>
    helper_test_read_only(
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
    expect_error("missing.*id")

  # missing house ID
  data.frame(id = "m1", hhid = NA) |>
    new_housepi_model("test") |>
    expect_error("missing.*hhid")

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
    expect_error("negative status")

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

  expect_equal(hhdata$start_date, c(1, NA, 1))
  expect_equal(hhdata$end_date, c(21, 21, NA))
})

test_that("A datum violating a model requirement is refused at import", {
  mk <- function(...) {
    data.frame(id = c("m1", "m2", "m3"), hhid = "h1",
               status = c(1, 1, 0), ...) |>
      new_housepi_model("case-ascertained", quiet = TRUE)
  }

  # 'case_ascertained' requires start_date and end_date discrete, and the
  # offending row is named alongside the column
  sms <- tryCatch(mk(start_date = c(1, 2.5, 1), end_date = 21),
                  error = conditionMessage)
  expect_match(sms, "not discrete-valued")
  expect_match(sms, "start_date")
  expect_match(sms, "row n.2")

  mk(start_date = 1, end_date = c(21, 21, 20.5)) |>
    expect_error("not discrete-valued")

  # integer-valued cells clear the requirement, negative ones included
  expect_no_error(mk(start_date = c(-3, 0, 1), end_date = 21))
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
  expect_null(hhmodel$data_wide)

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
  expect_s3_class(hhmodel$data_wide, "data.frame")
  expect_equal(nrow(hhmodel$data_wide), nrow(flu_1918_frederick))
  expect_true(is.finite(evaluate_log_lik(hhmodel)))
})

test_that("Data without hhid warns on every attach", {
  # no household structure: every row becomes a household of its own
  hhdata <- data.frame(id = c("m1", "m2", "m3"), status = c(1, 0, 1))

  expect_warning(
    new_housepi_model(hhdata, "test"), "[Ii]gnored household structure"
  )
  # a second model in the same session must warn again
  expect_warning(
    hhmodel <- new_housepi_model(hhdata, "test"),
    "[Ii]gnored household structure"
  )
  # and so must a later attach on an existing model
  expect_warning(add_data_wide(hhmodel, hhdata), "[Ii]gnored household")
})

test_that("Model class 'test' calls", {
  hhmodel <- new_housepi_model(model = "test")

  # 'test' default data wide format
  expect_equal(hhmodel$data_wide, dummy_household)

  # templates
  add_latent(hhmodel, import = "template") |>
    expect_warning("test")

  # zero log-likelihood ('test' model evaluates silently)
  expect_equal(evaluate_log_lik(hhmodel), 0.)

  # factorization round-trips through a fresh 'test' model
  hhdata <- sarscov2_2021_icofs
  hhdata$datum <- as.Date(runif(nrow(hhdata), max = 100))
  expect_equal(new_housepi_model(hhdata, "test")$data_wide, hhdata)

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
  expect_equal(names(snapshot), field_names)

  for (f in c("parameters", "group_risks", "latent",
              "latent_log_lik", "observation_log_lik")) {
    expect_length(snapshot[[f]], 0)
  }

  # single household
  expect_equal(snapshot$data_wide, dummy_household)
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
  expect_equal(get_data_long(hhmodel), hhdata)

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

  # help() dispatches without error in non-interactive mode (emits warning)
  expect_warning(hhmodel$help())
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
  expect_equal(hhmodel$data_wide, hhmodel_2$data_wide)
  expect_equal(hhmodel$data_long, hhmodel_2$data_long)

  # a copy's setup starts from the default configuration
  expect_equal(hhmodel_2$setup, new_housepi_model(house, "final-size")$setup)

  # definitions
  expect_equal(hhmodel$parameters, hhmodel_2$parameters)
  expect_equal(hhmodel$latent, hhmodel_2$latent)
  expect_equal(hhmodel$group_risks, hhmodel_2$group_risks)

  # changes on the original do not reach the copy
  hhmodel <- add_group_risk(hhmodel, "susceptibility", by = "member") |>
    suppressWarnings() # suppress baseline collapse warning
  expect_false(isTRUE(all.equal(hhmodel$group_risks, hhmodel_2$group_risks)))
})

#################################
## MODEL CLASS : SOURCE LOOKUP ##
#################################

test_that("Every declared name maps to the source declaring it", {
  # the reserved columns are not addressable as data, the others are
  expect_equal(dig_(mk_sources$model, "source_dictionary_"),
               c(datum = "D", record = "D", datum_2 = "D"))

  dict <- mk_sources$model |>
    add_parameter("rho", prior = "gamma", shape = 2, scale = 1) |>
    add_latent("hidden", model = "gamma", shape = 2, scale = 1) |>
    dig_("source_dictionary_")
  expect_equal(unname(dict[c("rho", "hidden")]), c("P", "A"))

  # a member-level observation shadows the like-named column: the read-cell is
  # wired to the producer, so the producer is not an ambiguity
  dict <- mk_sources$model |>
    add_latent("hidden", model = "gamma", shape = 2, scale = 1) |>
    add_observation("datum_2", model = "normal", mean = hidden@A, sd = 1) |>
    dig_("source_dictionary_")
  expect_equal(unname(dict["datum_2"]), "A")

  # a household-level observation has no latent key: it is only ever a column
  dict <- mk_sources$model |>
    add_observation("record_2", level = "house", value = min(datum_2@D)) |>
    dig_("source_dictionary_")
  expect_equal(unname(dict["record_2"]), "D")

  # two holders under one name: the marker says no sigil-free read settles it
  dict <- mk_sources$model |>
    add_parameter("datum", prior = "gamma", shape = 2, scale = 1) |>
    dig_("source_dictionary_")
  expect_equal(unname(dict["datum"]), "")
})

test_that("An unresolved source is reported softly and at the wall", {
  # nothing unresolved: both modes stay quiet
  expect_length(dig_(mk_sources$model, "assert_sources_", wall = FALSE), 0)
  expect_no_error(dig_(mk_sources$model, "assert_sources_", wall = TRUE))

  hhmodel <- mk_sources$model |>
    add_observation(record = list(model = "normal", mean = unknown, sd = 1))
  # soft: the messages are returned, and the model summary carries them
  expect_gt(length(dig_(hhmodel, "assert_sources_", wall = FALSE)), 0)
  expect_output(print(hhmodel), "Sources")
  # wall: the same messages are raised
  dig_(hhmodel, "assert_sources_", wall = TRUE) |>
    expect_error("[Cc]ould not resolve the source of .unknown.")
})

##################################
## MODEL CLASS : DATA SUMMARIES ##
##################################

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hhmodel <- new_housepi_model(flu_2014_matsumoto, "final-size")

  hhmodel |>
    summarise_house_sizes() |>
    expect_s3_class("data.frame")

  case_summary <- summarise_house_cases(hhmodel)
  # expect same household count as in Endo, et al. ('19)
  subset(case_summary, level == "household" & metric == "total")$count |>
    expect_equal(10486L)

  # house_size vs status equivalence groups correspond to measure non-zeros
  cases <- measure_cases_by_hhsize(flu_2014_matsumoto)
  cases$hhsize <- NULL
  subset(case_summary, level == "household" & metric == "classes")$count |>
    expect_equal(sum(cases[!is.na(cases)] > 0))

  # fix expected equivalence classes by member composition
  hhmodel |>
    add_group_risk("community", by = "member") |>
    summarise_house_cases() |>
    subset(level == "household" & metric == "classes") |>
    getElement("count") |>
    expect_equal(454L)
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
    run_inference(steps = 1e20, chains = 1e3) |>
    expect_error("memory budget")
})

test_that("Inference rejects invalid setup arguments", {
  mk_frederick$template |>
    run_inference(steps = 5L, burn = 10L) |>
    expect_error("[Bb]urn")
  mk_frederick$template |>
    run_inference(steps = 5L, thin = 10L) |>
    expect_error("[Tt]hin")
  mk_frederick$template |>
    run_inference(steps = 5L, adapt_windows = 10L) |>
    expect_error("[Uu]nexpected argument.*adapt_windows")
  mk_frederick$template |>
    run_inference(steps = 5L, chains = 1L, adapt = 0L, thin = 1L,
                  threads = 1000L) |>
    expect_warning("[Tt]hreads.*available") |>
    suppressWarnings()
})

test_that("A run with nothing to adapt keeps its sample size", {
  # the adaptation budget of a non-adaptive sampler is spent as burn-in
  new_housepi_model(model = "test") |>
    add_parameter(rho = list(prior = "uniform", proposal = "forward")) |>
    run_inference(steps = 5L, burn = 0L, adapt = 3L, thin = 1L, chains = 1L) |>
    suppressWarnings() |>
    get_parameter_draws() |>
    posterior::niterations() |>
    expect_equal(5L)
})

test_that("id_draws selections are resolved before the run", {
  hhmodel <- data.frame(id = c("m1", "m2"), hhid = "h1",
                        status = 1, obs = 0) |>
    new_housepi_model("test") |>
    add_latent("hidden", model = "uniform") |>
    helper_flat_sink(loc = hidden@A)

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
    expect_equal(0)
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

test_that("Thinning to ESS keeps between one draw per chain and all draws", {
  hhmodel <- data.frame(hhid = "h1", id = c("m1", "m2"), status = 0) |>
    new_housepi_model("test")
  # two disjoint trending chains of 20 iterations: ESS below the chain count
  trend <- posterior::as_draws_df(data.frame(
    rho = rep(c(10, 20), each = 20) + rep(seq_len(20L), 2) / 20,
    .chain = rep(1:2, each = 20), .iteration = rep(seq_len(20L), 2)
  ))
  hhmodel |>
    run_simulation(param_table = trend) |>
    expect_length(2L)
  # without thinning, every draw is simulated
  hhmodel |>
    run_simulation(param_table = trend, thin_ess = FALSE) |>
    expect_length(40L)
  # one chain of 20 iterations alternating low and high draws:
  # its ESS lies above the draw count, where `posterior` caps it
  alternate <- data.frame(rho = c(rbind(1:10, 20:11)))
  hhmodel |>
    run_simulation(param_table = posterior::as_draws_df(alternate)) |>
    expect_length(20L)
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
  nrow(get_parameters(hhmodel)) |>
    expect_equal(0L)
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing parameter definition")
})

test_that("Simulation returns combined wide/long output", {
  # supplying obs_process_long forces the long branch into a non-empty list
  obs_long <- function(x) data.frame(id = "m1", measure = 1.0)
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

test_that("Simulation replicates are a list of data frames", {
  simulated <- run_simulation(mk_frederick$template, nrep = 3)
  expect_type(simulated, "list")
  expect_length(simulated, 3)
  expect_s3_class(simulated[[1]], "data.frame")

  # a param_table larger than nrep is thinned down to it
  params <- get_parameters(mk_frederick$template)
  param_table <- as.data.frame(as.list(stats::setNames(params$value,
                                                       params$name)))
  param_table <- param_table[rep(1L, 10L), , drop = FALSE]
  run_simulation(mk_frederick$template, param_table = param_table, nrep = 4) |>
    expect_length(4L)
})

test_that("A simulation with an observation process is reproducible", {
  flag <- function(x) dplyr::mutate(dplyr::rowwise(x), flag = runif(1))
  simulated <- withr::with_seed(
    1, run_simulation(mk_frederick$template, obs_process = flag)
  )
  expect_true("flag" %in% colnames(simulated))
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
    apply_wide = function(x) list(x),
    apply_long = function(x) list(x),
    obs_proc_wide = NULL,
    obs_proc_long = NULL
  ) |>
    expect_error("obs_process_long")
})
