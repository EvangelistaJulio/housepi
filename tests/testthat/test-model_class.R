library(testthat)

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Model class inconsistent calls", {
  # non-test models may construct without data (held until the wall)
  new_housepi_model(model = "final_size") |>
    expect_silent()

  # test household model object
  hhmodel <- new_housepi_model(model = "test")

  # inference without definitions
  hhmodel |>
    run_inference() |>
    expect_error("[Nn]othing to sample")

  # data_long errors still surface via add_data_long
  hhmodel |>
    add_data_long("zero") |>
    expect_error("data.frame")

  # missing individual ID
  data.frame(id = NA, house_id = "h1") |>
    new_housepi_model("final_size") |>
    expect_error("missing.*id")

  # missing house ID
  data.frame(id = "m1", house_id = NA) |>
    new_housepi_model("final_size") |>
    expect_error("missing.*house_id")

  # not unique individual IDs
  data.frame(id = c("m1", "m1"), house_id = "h1") |>
    new_housepi_model("final_size") |>
    expect_error("not unique")

  # missing status: all
  data.frame(id = "m1", house_id = "h1", status = NA) |>
    new_housepi_model("final_size") |>
    add_parameter_template() |>
    run_inference() |>
    expect_error("missing status")
})

test_that("Wide data can be deferred and attached with add_data_wide", {
  # construct without data: silent, no data exposed
  hhmodel <- new_housepi_model(model = "final_size") |>
    add_parameter_template() |>
    expect_silent()
  expect_null(hhmodel$data_wide)
  expect_output(print(hhmodel))

  # the wall (evaluate/simulate) refuses a data-less model
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("[Nn]o household data")

  # attach data later: now populated and runnable
  snake <- add_data_wide(hhmodel, flu_2014_matsumoto)
  expect_s3_class(snake$data_wide, "data.frame")
  expect_equal(nrow(snake$data_wide), nrow(flu_2014_matsumoto))
  expect_true(is.finite(evaluate_log_lik(snake)))

  # bad inputs rejected by the wrapper
  add_data_wide(hhmodel, NULL) |>
    expect_error("data")
  add_data_wide(hhmodel, "ocelot") |>
    expect_error("data.frame")
})

test_that("Active fields are read-only", {
  new_housepi_model(model = "test") |>
    helper_test_read_only(
      "data_wide", "data_long", "parameters", "augmented", "infec_profile",
      "group_rates", "contact_rates", "name", "setup"
    )
})

test_that("Model class 'test' calls", {
  # test household model object
  hhmodel <- new_housepi_model(model = "test")
  snapshot <- hhmodel$export_house_snapshot("house")
  expect_equal(length(snapshot$param), 0)
  expect_equal(length(snapshot$augm), 0)
  expect_equal(length(snapshot$group_rate), 0)

  # 'test' default data wide format
  expect_equal(hhmodel$data_wide, dummy_household)
  expect_equal(snapshot$data_wide, hhmodel$data_wide)

  # correct data long format
  data_long <- data.frame(
    id = snapshot$data_wide$id[1],
    time = seq(3), measure = seq(3)**3
  )
  hhmodel <- add_data_long(hhmodel, data_long)
  expect_equal(hhmodel$data_long, data_long)
  # export_house_snapshot returns refactored data_long when present
  snap_long <- hhmodel$export_house_snapshot(snapshot$data_wide$house_id[1])
  expect_true(is.data.frame(snap_long$data_long))

  # factorization round-trips through a fresh 'test' model
  complex_data <- sarscov2_2021_icofs
  complex_data$dates <- as.Date(runif(dim(complex_data)[1], max = 100))
  complex_model <- new_housepi_model(complex_data, "test")
  expect_equal(complex_model$data_wide, complex_data)

  # templates
  add_parameter_template(hhmodel) |>
    expect_warning("test")
  add_augmented_template(hhmodel) |>
    expect_warning("test")

  # zero log-likelihood ('test' model evaluates silently)
  expect_equal(evaluate_log_lik(hhmodel), 0.)

  # not existing household id
  hhmodel |>
    evaluate_log_lik("mother-base") |>
    expect_error("[Ii]ndex.*mother-base.*did not match")

  # 'test' model simulation (runs silently)
  hhmodel <- new_housepi_model(model = "test")
  hhmodel <- add_augmented(hhmodel, name = "dummy", value = 3)
  out <- run_simulation(hhmodel)

  all(sapply(out$dummy[out$status == 0], function(x) identical(x, NA_real_))) |>
    expect_true()
})

test_that("Print behavior works as expected", {
  # with data and parameters
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter_template()
  expect_snapshot(print(hhmodel))

  # without parameters
  hhmodel2 <- flu_1918_frederick |>
    new_housepi_model("final_size")
  expect_snapshot(print(hhmodel2))

  # 'test' model
  expect_snapshot(print(new_housepi_model(model = "test")))

  # discrete-time flag branch and infectivity-profile section in print
  hhmodel_d <- sarscov2_2021_icofs |>
    dplyr::select(-start_date) |>
    new_housepi_model("case_ascertained", discr = TRUE, quiet = TRUE) |>
    add_house_observation(start_date = min(test_date@D))
  expect_snapshot(print(hhmodel_d))

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
    add_group_rate("community", by = "member") |>
    summarise_house_cases()
  subset(case_summary, level == "household" & metric == "distinct")$count |>
    expect_equal(454L)
})

######################################
## MODEL CLASS : POSTERIOR SAMPLING ##
######################################

test_that("Posterior memory budget is asserted", {
  house <- data.frame(house_id = 1, id = "a", status = 0)
  hhmodel <- new_housepi_model(house, "final_size")
  hhmodel <- add_parameter_template(hhmodel)
  hhmodel |>
    run_inference(steps = 1e20, chains = 1e3, progress = FALSE) |>
    expect_error("memory budget")
})

test_that("A followed household without infected members warns", {
  hhmodel <- data.frame(id = c("chico", "amanda"),
                        house_id = c("ground-zero", "mother-base"),
                        status = c(1, 0)) |>
    new_housepi_model("test") |>
    add_augmented("test_date", prior = "uniform")

  expect_warning(
    expect_warning(
      hhmodel |>
        run_inference(steps = 20, adapt = 0, chains = 1,
                      id_draws = list(house_id = "mother-base")),
      regexp = "no infected members"
    ),
    regexp = "below 400"
  )
})

test_that("Rlang calls work as expected", {
  # declare indiv_shape as an augmented variable so references resolve
  hhmodel <- new_housepi_model(model = "test") |>
    add_augmented("indiv_shape", prior = "normal")

  hhmodel |>
    add_parameter("kaz", prior = "normal", mean = indiv_shape + 1) |>
    expect_error("[Ii]ndividual specific")

  val <- 1
  # UPDATE: for the moment I am forcing parameter arguments to be co
  hhmodel <- hhmodel |>
    add_parameter(
      param_shape = 5,
      var1 = list(prior = "normal", mean = param_shape@p + 1),
      var2 = list(prior = "normal", mean = param_shape@p + !!val),
      var3 = list(prior = "normal", mean = ~ param_shape@P + 1),
      var4 = list(prior = "normal", mean = ~ param_shape@P + !!val),
      append = FALSE
    )
  aux <- hhmodel$parameters
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 2L)

  val <- 3
  hhmodel <- hhmodel |>
    add_augmented(
      var1 = list(prior = "normal", mean = indiv_shape + 3),
      var2 = list(prior = "normal", mean = indiv_shape + !!val),
      var3 = list(prior = "normal", mean = ~ indiv_shape + 3),
      var4 = list(prior = "normal", mean = ~ indiv_shape + !!val)
    )
  aux <- hhmodel$augmented
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 2L)
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

test_that("run_measure_ rejects apply_long without obs_process_long", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(alpha = 1)
  priv <- hhmodel$.__enclos_env__$private
  pt <- data.frame(alpha = c(0.1, 0.2, 0.3))

  priv$run_measure_(param_table = pt,
                    apply_wide = function(x) list(x),
                    apply_long = function(x) list(x),
                    obs_proc_wide = NULL,
                    obs_proc_long = NULL,
                    progress = FALSE) |>
    expect_error("obs_process_long")
})

test_that("run_simulation returns combined wide/long output", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter_template()
  # supplying obs_process_long forces the long branch into a non-empty list
  # single simulation by default
  obs_long <- function(x) data.frame(id = "vic-boss", measure = 1.0)
  out <- hhmodel |>
    run_simulation(obs_process = function(x) x,
                   obs_process_long = obs_long)
  expect_true(is.list(out))
  expect_equal(names(out), c("wide", "long"))
  expect_true(is.data.frame(out$wide))
  expect_true(is.data.frame(out$long))

  # multiple repetitions
  out <- hhmodel |>
    run_simulation(obs_process = function(x) x,
                   obs_process_long = obs_long,
                   nrep = 2)
  expect_true(is.data.frame(out$wide[[1]]))
  expect_true(is.data.frame(out$long[[2]]))
})
