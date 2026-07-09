library(testthat)

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Model class inconsistent calls", {
  # non-test models require data_wide at construction
  new_housepi_model(model = "final_size") |>
    expect_error("data")

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

  # discrete-time flag branch in print
  hhmodel_d <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained", discr = TRUE, quiet = TRUE)
  expect_snapshot(print(hhmodel_d))

  # infectivity-profile section in print
  hhmodel_ip <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |> 
    add_house_observ(start_date = min(test_date@D))
  expect_snapshot(print(hhmodel_ip))

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

test_that("Wrong definitions throw errors", {
  hhmodel <- new_housepi_model(model = "test")

  # eager rebuild: errors surface at definition, not at run_inference
  add_parameter(hhmodel, "a", prior = "uniform", min = ~ min_val@D) |>
    expect_error("not supported")
  # hint is thrown
  # FIXME: no-append exposes bug: faulty definition "a" stays stored!
  # add_parameter(hhmodel, b = list(prior = "uniform", min = min_val)) |>
  #   expect_error("rlang injection")

  # memory budget is exceeded
  house <- data.frame(house_id = 1, id = "a", status = 0)
  hhmodel <- new_housepi_model(house, "final_size")
  hhmodel <- add_parameter_template(hhmodel)
  hhmodel |>
    run_inference(steps = 1e20, chains = 1e3, progress = FALSE) |>
    expect_error("memory budget")
})

test_that("Discretization is enabled as expected", {
  discr <- TRUE
  val <- list(min = 1, max = 5)
  hhmodel <- new_housepi_model(model = "test")
  hhmodel <- add_parameter(hhmodel,
    sample = list(
      prior = "gamma", shape = 2, scale = 1, discr = TRUE,
      proposal = "sample", min = !!val$min, max = !!val$max
    ),
    uniform = list(
      prior = "gamma", shape = 2, scale = 1, discr = TRUE,
      proposal = "uniform", min = !!val$min, max = !!val$max
    ),
    normal = list(
      prior = "gamma", shape = 2, scale = 1, discr = TRUE,
      proposal = "normal", min = !!val$min, max = !!val$max
    ),
    lognormal = list(
      prior = "gamma", shape = 2, scale = 1, discr = TRUE,
      proposal = "lognormal", min = !!val$min, max = !!val$max
    )
  )

  results <- suppressWarnings(run_inference(hhmodel,
    steps = 50, burn = 0, thin = 1, adapt = 50, chains = 1
  ))
  y <- posterior::as_draws_matrix(get_parameter_draws(results))
  expect_true(all(y == floor(y)))
  expect_true(all(y >= val$min))
  expect_true(all(y <= val$max))
})

test_that("Invalid proposal name is rejected", {
  new_housepi_model(model = "test") |>
    add_parameter("snake", prior = "uniform", proposal = "venom") |>
    run_inference(steps = 1, chains = 1, progress = FALSE) |>
    expect_error("Unknown proposal.*venom")
})

test_that("Degenerate range raises an error for parameters", {
  new_housepi_model(model = "test") |>
    add_parameter("venom", prior = "uniform", min = 3, max = 3) |>
    expect_error("degenerate range")
})

test_that("Lognormal proposal non-positive domain", {
  # prior with unbounded lower bound
  new_housepi_model(model = "test") |>
    add_parameter("huey", prior = "normal", proposal = "lognormal") |>
    expect_error("non-positive.*huey")

  # prior with negative lower bound
  new_housepi_model(model = "test") |>
    add_parameter("miller", prior = "uniform",
                  min = -1, proposal = "lognormal") |>
    expect_error("non-positive")

  # location parameter
  suppressWarnings(
    param_draws <- new_housepi_model(model = "test") |>
      add_parameter("amanda", prior = "gamma", shape = 2,
                    max = 5, loc = -1, proposal = "lognormal") |>
      run_inference(steps = 20, chains = 1, adapt = 0, thin = 1) |>
      get_parameter_draws()
  )
  expect_true(all(param_draws[["amanda"]] > -1))
  expect_true(all(param_draws[["amanda"]] < 4))
})

#########################################
## MODEL CLASS : AUGMENTED DEFINITIONS ##
#########################################

test_that("Modules do not get overwritten", {
  # NOTE: I use profile as a  status flag for models
  out_1 <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    get_infec_profile()
  expect_gt(length(out_1), 0L)

  flu_1918_frederick |>
    new_housepi_model("final_size") |>
    get_infec_profile() |>
    expect_length(0L)

  # reload first model
  sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    get_infec_profile() |>
    expect_identical(out_1)
})

test_that("Custom models shares function registries", {
  hhmodel_0 <- new_housepi_model(model = "test")
  expose_function("growth") |>
    expect_error("Unknown function")

  # compile example
  custom <- system.file(file.path("examples", "custom"), package = "housepi")
  hhmodel_1 <- flu_1918_frederick |>
    new_housepi_model(custom, quiet = TRUE)

  # no templates available for custom example
  hhmodel_1 |>
    add_augmented_template() |>
    expect_error("file.*not found")

  # Growth_ overrides only raw_function_ and raw_primitive_
  # other raw_ virtuals fall back to raise_undefined_call_err_
  x <- seq(0, 3, 0.5)
  initial <- 2
  rate <- 0.5
  fn <- expose_function("growth", initial = initial, rate = rate)
  expect_equal(fn$func(x), initial * exp(rate * x))
  expect_equal(fn$prim(x), initial * exp(rate * x) / rate)
  fn$log_func(x) |>
    expect_error("undefined")
  # default raw_log_ratio_ returns 0
  fn$log_ratio(x, x + 0.1) |>
    expect_equal(rep(0L, length(x)))
})

#########################################
## MODEL CLASS : AUGMENTED DEFINITIONS ##
#########################################

test_that("Rlang calls work as expected", {
  # declare indiv_shape as an augmented variable so references resolve
  hhmodel <- new_housepi_model(model = "test") |>
    add_augmented("indiv_shape", prior = "normal")

  add_parameter(hhmodel, "kaz", prior = "normal", mean = indiv_shape + 1) |>
    expect_error("[Ii]ndividual specific")

  val <- 1
  # UPDATE: for the moment I am forcing parameter arguments to be co
  hhmodel <- add_parameter(hhmodel,
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
  hhmodel <- add_augmented(hhmodel,
    var1 = list(prior = "normal", mean = indiv_shape + 3),
    var2 = list(prior = "normal", mean = indiv_shape + !!val),
    var3 = list(prior = "normal", mean = ~ indiv_shape + 3),
    var4 = list(prior = "normal", mean = ~ indiv_shape + !!val),
    append = TRUE
  )
  aux <- hhmodel$augmented
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 2L)
})

test_that("Augmented aliasing works as expected", {
  # check that model arguments are correctly integrated
  simulated <- data.frame(id = seq(1, 1000), house_id = 1, status = 1) |>
    new_housepi_model("test") |>
    add_augmented(
      hour = list(prior = "uniform", min = day, max = day + 1),
      day = list(prior = "gamma", shape = 2, scale = 1)
    ) |>
    run_simulation()
  expect_true(all(
    simulated$day <= simulated$hour & simulated$hour <= simulated$day + 1
  ))
})

test_that("Nested definitions works as expected", {
  hhmodel <- new_housepi_model(model = "test") |>
    # define Box-Muller transformation
    add_augmented(norm = list(prior = "normal"),
                  u1 = list(prior = "uniform"),
                  u2 = list(prior = "uniform"),
                  z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2))
  expect_snapshot(print(hhmodel))

  results <- run_inference(hhmodel, steps = 1e5, chains = 1, adapt = 0,
                           thin = 1, progress = 0)
  # no parameters nor priors, only augmented
  expect_equal(posterior::nvariables(results$parameter_draws), 0L)
  expect_null(results$.__enclos_env__$private$param_prior_draws_)
  expect_equal(posterior::nvariables(results$augmented_draws),
               dim(hhmodel$augmented)[1])

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

test_that("Load augmented list works fine for external lists", {

  hhmodel <- data.frame(id = "id1", house_id = "h1",
                        status = 1, date = as.Date("2021-01-01")) |>
    new_housepi_model("test")

  # an external list of definitions can be imported
  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  tmp <- add_augmented(hhmodel, import = external_list)$augmented

  hhmodel <- add_augmented(hhmodel, z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})

test_that("Augmented traces export only relevant variables", {
  hhmodel <- new_housepi_model(model = "test") |>
    add_parameter(alpha = 1, beta = list(prior = "uniform")) |>
    add_augmented(
      # constants
      blocked = 42,
      # private variables
      .hidden = list(prior = "uniform"),
      .secret = list(prior = "normal"),
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
    add_augmented("d_horse", prior = "uniform", min = 5, max = 5)
  expect_warning(
    # FIXME double check this
    tryCatch(
      run_inference(hhmodel, steps = 10, chains = 1, adapt = 0, thin = 1),
      error = function(e) NULL
    ),
    regexp = "[Dd]egenerate range|[Cc]ollapsed"
  )
})

############################################
## MODEL CLASS : CONTACT RATE DEFINITIONS ##
############################################

test_that("Contact rate definition work as expected", {
  hhid <- "hh2"
  hhmodel <- data.frame(
    house_id = c(rep("hh1", 3), rep(hhid, 6)),
    status = c(1, 0, 0, 1, 0, 0, 0, 0, 0),
    member = c("student", "other", "other",
               "student", "sibling", "mother", "father", "other", "student")
  ) |>
    new_housepi_model("final_size") |>
    add_parameter_template()

  # no contact rate defs -> empty matrix
  hhsnapshot <- hhmodel$export_house_snapshot(hhid)
  cr_mtx <- hhsnapshot$contact_rate
  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), 0)
  expect_equal(ncol(cr_mtx), 0)

  # mutually exclusive strata
  strata <- data.frame(
    from_member = c(
      "student", "student", "sibling", "mother",
      "mother", "father", "other", NA
    ),
    to_member = c(
      "student", "sibling", "sibling", "student",
      "sibling", "mother", "other", NA
    ),
    param_name = c(
      "child_child", "child_child", "child_child", "mother_child",
      "mother_child", ".reference", "other_other", "cross_gen"
    )
  )
  # prime numbers as param values to ensure non duplicate pair productes
  params <- c(2, 3, 5, 7, 11, 13, 17, 19)
  names(params) <- strata$param_name
  hhmodel <- hhmodel |>
    add_contact_rate(strata = strata, params = params)

  # matrix is n x n
  hhsnapshot <- hhmodel$export_house_snapshot(hhid)
  cr_mtx <- hhsnapshot$contact_rate
  n <- nrow(hhsnapshot$data_wide)

  expect_true(is.matrix(cr_mtx))
  expect_equal(nrow(cr_mtx), n)
  expect_equal(ncol(cr_mtx), n)
  # diagonal is 0
  expect_equal(unname(diag(cr_mtx)), rep(0., n))

  # every off-diagonal cell resolves to the reference (1)
  # or to exactly one param value: no multiplications
  expect_true(all(cr_mtx %in% c(0, 1, params)))

  # drop_contact_rates_ cleans up the associated parameter definitions
  param_names_before <- hhmodel$parameters$name
  hhmodel$.__enclos_env__$private$drop_contact_rates_()
  any(setdiff(names(params), ".reference") %in% hhmodel$parameters$name) |>
    expect_false()
  # the template parameters added previously must survive
  template_pars <- c("alpha", "beta", "gamma")
  expect_true(all(template_pars %in% param_names_before))
})

test_that("Contact rates raise an error for incompatible models", {
  data.frame(
    house_id = 1, id = seq(4), status = seq(4),
    cohort = c("a", "a", "b", "b")
  ) |>
    new_housepi_model("test") |>
    add_contact_rate(by = "cohort") |>
    expect_error("not supported")
})

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

  hhmodel <- new_housepi_model(model = "test")
  hhmodel |>
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
  hhmodel <- new_housepi_model(model = "test")
  hhmodel <- add_parameter(hhmodel, alpha = 1)
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
  out <- run_simulation(hhmodel, obs_process = function(x) x,
                        obs_process_long = obs_long)
  expect_true(is.list(out))
  expect_equal(names(out), c("wide", "long"))
  expect_true(is.data.frame(out$wide))
  expect_true(is.data.frame(out$long))

  # multiple repetitions
  out <- run_simulation(hhmodel, obs_process = function(x) x,
                        obs_process_long = obs_long,
                        nrep = 2)
  expect_true(is.data.frame(out$wide[[1]]))
  expect_true(is.data.frame(out$long[[2]]))
})

test_that("add_augmented_template loads", {
  sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_augmented_template() |>
    get_augmented() |>
    nrow() |>
    expect_gt(0L)
})

test_that("add_infec_profile and list update the model directly", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("case_ascertained") |>
    add_infec_profile(func = "gamma", shape = 2, scale = 1)
  expect_true("gamma" %in% hhmodel$infec_profile$func)
  expect_warning(
    hhmodel <- add_infec_profile(hhmodel, func = "lognormal"),
    regexp = "Overwriting"
  )
  expect_true("lognormal" %in% hhmodel$infec_profile$func)
})

test_that("Parameter defined as a function tracks its dependencies", {
  set.seed(1)
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(
      alpha = 1e-4,
      beta = hhsize_ref@P^-gamma@P,
      gamma = list(prior = "uniform"),
      hhsize_ref = 4
    )
  res <- suppressWarnings(run_inference(hhmodel,
    steps = 100, chains = 1, adapt = 0, thin = 1
  ))
  vars <- posterior::variables(res$parameter_draws)
  # both sampled and function-derived parameters appear in the trace
  expect_true(all(c("gamma", "beta") %in% vars))
  # the function-derived parameter has no prior contribution
  # UPDATE: this will not require private access in the future
  expect_false("lprior_beta" %in% posterior::variables(
    res$.__enclos_env__$private$param_prior_draws_)
  )

  df <- as.data.frame(posterior::as_draws_df(res$parameter_draws))
  # every draw of beta equals hhsize_ref^(-gamma) exactly
  expect_equal(df$beta, 4^(-df$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size")

  # bare names default to augmented source! here @P is missing
  hhmodel |>
    add_parameter(alpha = 1e-4, gamma = 1, hhsize_ref = 4,
                  beta = hhsize_ref^-gamma) |>
    expect_error("cannot depend on augmented variables")

  # explicit reference to a data column is also rejected
  # FIXME: roll-back missing -> append = FALSE should not be necessary
  hhmodel |>
    add_parameter(alpha = 1e-4, gamma = 1, hhsize_ref = 4,
                  beta = status@D + 1, append = FALSE) |>
    expect_error("cannot depend on data variables")
})

test_that("Parameter defined as a function tracks its dependencies", {
  set.seed(1)
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter(alpha = 1e-4,
                  beta = hhsize_ref@P^-gamma@P,
                  gamma = list(prior = "uniform"),
                  hhsize_ref = 4)
  res <- run_inference(hhmodel, steps = 2000, progress = 0)
  vars <- posterior::variables(res$parameter_draws)
  # both sampled and function-derived parameters appear in the trace
  expect_true(all(c("gamma", "beta") %in% vars))
  # the function-derived parameter has no prior contribution
  expect_false("lprior_beta" %in% posterior::variables(
    res$.__enclos_env__$private$param_prior_draws_))

  df <- as.data.frame(posterior::as_draws_df(res$parameter_draws))
  # snake on the trail: every draw of beta equals hhsize_ref^(-gamma) exactly
  expect_equal(df$beta, 4^(-df$gamma))
})
