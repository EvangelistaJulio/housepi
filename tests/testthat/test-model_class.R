library(testthat)

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Active fields are read-only", {
  new_housepi_model(model = "test") |>
    helper_test_read_only(
      "data_wide", "data_long", "parameters", "augmented", "observations",
      "infectivity_profiles", "group_risks", "contact_risks", "name", "setup"
    )
})

test_that("Model class inconsistent IDs throw errors", {
  # missing individual ID
  data.frame(id = NA, house_id = "h1") |>
    new_housepi_model("test") |>
    expect_error("missing.*id")

  # missing house ID
  data.frame(id = "m1", house_id = NA) |>
    new_housepi_model("test") |>
    expect_error("missing.*house_id")

  # not unique individual IDs
  data.frame(id = c("m1", "m1"), house_id = "h1") |>
    new_housepi_model("test") |>
    expect_error("not unique")
})

test_that("Data status is capped by the engine status storage", {
  # above the highest storable infected status
  data.frame(house_id = "h1", status = 255) |>
    new_housepi_model("test") |>
    expect_error("maximum of 254")

  # a negative status has no meaning: NA declares a missing status
  data.frame(house_id = "h1", status = -200) |>
    new_housepi_model("test") |>
    expect_error("negative status")

  # a missing status round trips as NA, the top level is never exported
  df <- data.frame(house_id = "h1", status = c(NA, 254)) |>
    new_housepi_model("test") |>
    get_data_wide()
  expect_equal(df$status, c(NA, 254))
})

test_that("data columns cannot claim the hidden-definition prefix", {
  df <- dummy_household
  df$.tx <- 55

  new_housepi_model(df, "test") |>
    expect_error("with a dot")
  new_housepi_model(model = "test") |>
    add_data_wide(df) |>
    expect_error("with a dot")

  new_housepi_model(model = "test") |>
    add_data_long(data.frame(id = "m1", .shagohod = 1)) |>
    expect_error("with a dot")
})

test_that("Wide data can be deferred and attached with add_data_wide", {
  # construct without data: silent, no data exposed
  hhmodel <- new_housepi_model(model = "final_size") |>
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
    add_data_wide("dead_cell") |>
    expect_error("data.frame")

  # attach data later: now populated and runnable
  hhmodel <- hhmodel |>
    add_data_wide(flu_1918_frederick)
  expect_s3_class(hhmodel$data_wide, "data.frame")
  expect_equal(nrow(hhmodel$data_wide), nrow(flu_1918_frederick))
  expect_true(is.finite(evaluate_log_lik(hhmodel)))
})

test_that("Data without house_id warns on every attach", {
  # no household structure: every row becomes a household of its own
  zanzibar <- data.frame(id = c("frank", "gray", "fox"), status = c(1, 0, 1))

  expect_warning(
    new_housepi_model(zanzibar, "test"), "[Ii]gnored household structure"
  )
  # a second model in the same session must warn again
  expect_warning(
    new_housepi_model(zanzibar, "test"), "[Ii]gnored household structure"
  )
  # and so must a later attach on an existing model
  hhmodel <- suppressWarnings(new_housepi_model(zanzibar, "test"))
  expect_warning(add_data_wide(hhmodel, zanzibar), "[Ii]gnored household")
})

test_that("Model class 'test' calls", {
  # test household model object
  hhmodel <- new_housepi_model(model = "test")

  # 'test' default data wide format
  expect_equal(hhmodel$data_wide, dummy_household)

  # templates
  add_parameter(hhmodel, import = "template") |>
    expect_warning("test")
  add_augmented(hhmodel, import = "template") |>
    expect_warning("test")

  # zero log-likelihood ('test' model evaluates silently)
  expect_equal(evaluate_log_lik(hhmodel), 0.)

  # factorization round-trips through a fresh 'test' model
  complex_data <- sarscov2_2021_icofs
  complex_data$dates <- as.Date(runif(dim(complex_data)[1], max = 100))
  complex_model <- new_housepi_model(complex_data, "test")
  expect_equal(complex_model$data_wide, complex_data)

  # not existing household id
  hhmodel |>
    evaluate_log_lik("mother-base") |>
    expect_error("[Ii]ndex.*mother-base.*did not match")
})

test_that("Household snapshot reads before any definition is pushed", {
  # entries required by the model are kept even when still empty
  hhmodel <- dummy_household |>
    new_housepi_model("final_size")
  snapshot <- hhmodel$export_house_snapshot("house")

  field_names <- c(
    "data_wide", "parameters", "group_risks", "contact_risks", "augmented",
    "augmented_log_lik", "observation_log_lik", "transmission_log_lik"
  )
  expect_equal(names(snapshot), field_names)

  for (f in c("parameters", "group_risks", "augmented",
              "augmented_log_lik", "observation_log_lik")) {
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

  # correct data long format
  data_long <- data.frame(id = hhmodel$data_wide$id[1],
                          time = seq(3), measure = seq(3)**3)

  hhmodel <- hhmodel |>
    add_data_long(data_long)
  expect_equal(hhmodel$data_long, data_long)

  # export_house_snapshot returns refactored data_long when present
  snap_long <- hhmodel$export_house_snapshot(hhmodel$data_wide$house_id[1])
  expect_true(is.data.frame(snap_long$data_long))
})

test_that("Print behavior works as expected", {
  # 'test' model
  print(new_housepi_model(model = "test")) |>
    expect_snapshot()

  # without data
  print(new_housepi_model(model = "final_size")) |>
    expect_snapshot()

  # without parameters
  print(new_housepi_model(flu_1918_frederick, "final_size")) |>
    expect_snapshot()

  # with data and parameters
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template")
  print(hhmodel) |>
    expect_snapshot()

  # help() dispatches without error in non-interactive mode (emits warning)
  expect_warning(hhmodel$help())
})

##################################
## MODEL CLASS : DATA SUMMARIES ##
##################################

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hhmodel <- flu_2014_matsumoto |>
    new_housepi_model("final_size")

  hhmodel |>
    summarise_house_sizes() |>
    is.data.frame() |>
    expect_true()

  case_summary <- summarise_house_cases(hhmodel)
  expect_true(is.data.frame(case_summary))
  # expect same household count as in Endo, et al. ('19)
  subset(case_summary, level == "household" & metric == "total")$count |>
    expect_equal(10486L)

  # house_size vs status equivalence groups correspond to measure non-zeros
  aux <- measure_cases_by_hhsize(flu_2014_matsumoto)
  aux$hhsize <- NULL
  subset(case_summary, level == "household" & metric == "distinct")$count |>
    expect_equal(sum(aux[!is.na(aux)] > 0))

  # fix expected equivalence classes by member composition
  case_summary <- hhmodel |>
    add_group_risk("community", by = "member") |>
    summarise_house_cases()
  subset(case_summary, level == "household" & metric == "distinct")$count |>
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
  data.frame(id = "m1", house_id = "h1", status = NA) |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template") |>
    run_inference() |>
    expect_error("missing status")
})

test_that("Posterior memory budget is asserted", {
  house <- data.frame(house_id = 1, id = "a", status = 0)
  hhmodel <- new_housepi_model(house, "final_size")
  hhmodel <- add_parameter(hhmodel, import = "template")
  hhmodel |>
    run_inference(steps = 1e20, chains = 1e3) |>
    expect_error("memory budget")
})

test_that("A followed household without infected members warns", {
  hhmodel <- data.frame(id = c("chico", "amanda"),
                        house_id = c("ground-zero", "mother-base"),
                        status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_augmented("test_date", model = "uniform")

  expect_warning(
    expect_warning(
      hhmodel |>
        run_inference(steps = 20L, adapt = 0L, chains = 1L,
                      id_draws = list(house_id = "mother-base")),
      regexp = "no infected members"
    ),
    regexp = "below 400"
  )
})

##################################
## MODEL CLASS : SIMULATION RUN ##
##################################

test_that("run_simulation rejects wrong arguments", {
  hhmodel <- new_housepi_model(model = "test")
  # wrong observation process
  hhmodel |>
    run_simulation(obs_process = "not a function") |>
    expect_error("function")
  # wrong number of replicates
  hhmodel |>
    run_simulation(nrep = 0) |>
    expect_error("positive integer")
  hhmodel |>
    run_simulation(nrep = -3) |>
    expect_error("positive integer")
})

test_that("Param table with degenerate columns raises an ESS error", {
  pt <- posterior::as_draws_df(data.frame(p = rep(1, 100)))

  hhmodel <- new_housepi_model(model = "test") |>
    run_simulation(param_table = pt) |>
    expect_error("ESS")
})

test_that("run_simulation auto-imports parameters from param_table", {
  hhmodel <- data.frame(house_id = 1, id = letters[1:3], status = 0) |>
    new_housepi_model("test")
  expect_equal(nrow(hhmodel$parameters), 0)
  pt <- data.frame(adamska = c(0.1, 0.2, 0.3))
  invisible(run_simulation(hhmodel, param_table = pt))
  expect_true("adamska" %in% hhmodel$parameters$name)
})

test_that("run_simulation returns combined wide/long output", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(import = "template")
  # supplying obs_process_long forces the long branch into a non-empty list
  # single simulation by default
  obs_long <- function(x) data.frame(id = "vic-boss", measure = 1.0)
  out <- hhmodel |>
    run_simulation(
      obs_process = function(x) x,
      obs_process_long = obs_long
    )
  expect_true(is.list(out))
  expect_equal(names(out), c("wide", "long"))
  expect_true(is.data.frame(out$wide))
  expect_true(is.data.frame(out$long))

  # multiple repetitions
  out <- hhmodel |>
    run_simulation(
      obs_process = function(x) x,
      obs_process_long = obs_long,
      nrep = 2
    )
  expect_true(is.data.frame(out$wide[[1]]))
  expect_true(is.data.frame(out$long[[2]]))
})

test_that("run_measure_ rejects apply_long without obs_process_long", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(alpha = 1)
  pvt <- hhmodel$.__enclos_env__$private

  pvt$run_measure_(
    param_table = data.frame(alpha = c(0.1, 0.2, 0.3)),
    apply_wide = function(x) list(x),
    apply_long = function(x) list(x),
    obs_proc_wide = NULL,
    obs_proc_long = NULL
  ) |>
    expect_error("obs_process_long")
})
