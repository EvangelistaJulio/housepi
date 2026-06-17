library(testthat)

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Model class inconsistent calls", {
  # non-blank models require data_wide at construction
  new_housepi_model(model = "chain_binomial") |>
    expect_error("data")

  # blank household model object
  hhmodel <- new_housepi_model(model = "blank")

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
    new_housepi_model("chain_binomial") |>
    expect_error("missing.*id")

  # missing house ID
  data.frame(id = "m1", house_id = NA) |>
    new_housepi_model("chain_binomial") |>
    expect_error("missing.*house_id")

  # not unique individual IDs
  data.frame(id = c("m1", "m1"), house_id = "h1") |>
    new_housepi_model("chain_binomial") |>
    expect_error("not unique")

  # missing status: all
  data.frame(id = "m1", house_id = "h1", status = NA) |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template() |>
    run_inference() |>
    expect_error("missing status")
})

test_that("Active fields are read-only", {
  new_housepi_model(model = "blank") |>
    helper_test_read_only(
      "data_wide", "data_long", "parameters", "augmented", "infec_profile",
      "group_rates", "contact_rates", "name", "setup"
    )
})

test_that("Model class 'blank' calls", {
  # blank household model object
  hhmodel <- new_housepi_model(model = "blank")
  snapshot <- hhmodel$export_house_snapshot("house")
  expect_equal(length(snapshot$param), 0)
  expect_equal(length(snapshot$augm), 0)
  expect_equal(length(snapshot$group_rate), 0)

  # blank default data wide format
  expect_equal(hhmodel$data_wide, blank_house)
  expect_equal(snapshot$data_wide, hhmodel$data_wide)

  # correct data long format
  data_long <- data.frame(
    id = snapshot$data_wide$id[1],
    time = seq(3), measure = seq(3)**3
  )
  hhmodel$add_data_long(data_long)
  expect_equal(hhmodel$data_long, data_long)
  # export_house_snapshot returns refactored data_long when present
  snap_long <- hhmodel$export_house_snapshot(snapshot$data_wide$house_id[1])
  expect_true(is.data.frame(snap_long$data_long))

  # factorization round-trips through a fresh blank model
  complex_data <- sarscov2_2021_icofs
  complex_data$dates <- as.Date(runif(dim(complex_data)[1], max = 100))
  complex_model <- new_housepi_model(complex_data, "blank")
  expect_equal(complex_model$data_wide, complex_data)

  # templates
  hhmodel$add_parameter_template() |>
    expect_warning("blank")
  hhmodel$add_augmented_template() |>
    expect_warning("blank")

  # zero log-likelihood
  expect_warning(x <- eval_log_density(hhmodel), regexp = "blank")
  expect_equal(x, 0)

  # not existing household id
  hhmodel$evaluate_log_lik("mother-base") |>
    expect_error("[Ii]ndex.*mother-base.*did not match")

  # blank simulation
  hhmodel <- new_housepi_model(model = "blank")
  hhmodel$add_augmented(name = "dummy", value = 3)
  expect_warning(out <- hhmodel$run_simulation(), regexp = "blank")
  expect_true(
    all(sapply(out$dummy[out$status == 0], function(x) identical(x, NA_real_)))
  )
})

test_that("Print behavior works as expected", {
  # with data and parameters
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  expect_snapshot(print(hhmodel))

  # without parameters
  hhmodel2 <- flu_1918_frederick |>
    new_housepi_model("chain_binomial")
  expect_snapshot(print(hhmodel2))

  # blank model
  expect_snapshot(print(new_housepi_model(model = "blank")))

  # discrete-time flag branch in print
  hhmodel_d <- sarscov2_2021_icofs |>
    new_housepi_model("augmented_infection_time", discr = TRUE, quiet = TRUE)
  expect_snapshot(print(hhmodel_d))

  # infectivity-profile section in print
  hhmodel_ip <- sarscov2_2021_icofs |>
    new_housepi_model("augmented_infection_time")
  expect_snapshot(print(hhmodel_ip))

  # help() dispatches without error in non-interactive mode (emits warning)
  expect_warning(hhmodel$help())
})

##################################
## MODEL CLASS : DATA SUMMARIES ##
##################################

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hhmodel <- flu_2014_matsumoto |>
    new_housepi_model("chain_binomial")

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
  aux <- measure_cases_by_house_size(flu_2014_matsumoto)
  aux$house_size <- NULL
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
  hhmodel <- new_housepi_model(model = "blank")

  # eager rebuild: errors surface at definition, not at run_inference
  hhmodel$add_parameter("a", prior = "uniform", min = ~ min_val@D) |>
    expect_error("not supported")
  # hint
  hhmodel$add_parameter_list(b = list(prior = "uniform", min = min_val)) |>
    expect_error("rlang injection")
  # FIXME : this also should fire the same:
  # hhmodel$add_parameter("b", prior = "uniform", min = min_val) |>
  #   expect_error("rlang injection")

  # memory budget is exceeded
  house <- data.frame(house_id = 1, id = "a", status = 0)
  hhmodel <- new_housepi_model(house, "chain_binomial")
  hhmodel$add_parameter_template()
  hhmodel$run_inference(steps = 1e20, chains = 1e3, progress = FALSE) |>
    expect_error("memory budget")
})

test_that("Discretization is enabled as expected", {
  discr <- TRUE
  val <- list(min = 1, max = 5)
  hhmodel <- new_housepi_model(model = "blank")
  hhmodel$add_parameter_list(
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

  results <- suppressWarnings(hhmodel$run_inference(
    steps = 50, burn = 0, thin = 1, adapt = 50, chains = 1
  ))
  y <- posterior::as_draws_matrix(get_parameter_draws(results))
  expect_true(all(y == floor(y)))
  expect_true(all(y >= val$min))
  expect_true(all(y <= val$max))
})

test_that("Invalid proposal name is rejected", {
  new_housepi_model(model = "blank") |>
    add_parameter("snake", prior = "uniform", proposal = "venom") |>
    run_inference(steps = 1, chains = 1, progress = FALSE) |>
    expect_error("Unknown proposal.*venom")
})

test_that("Degenerate range raises an error for parameters", {
  new_housepi_model(model = "blank") |>
    add_parameter("venom", prior = "uniform", min = 3, max = 3) |>
    expect_error("degenerate range")
})

test_that("Lognormal proposal non-positive domain", {
  # prior with unbounded lower bound
  new_housepi_model(model = "blank") |>
    add_parameter("huey", prior = "normal", proposal = "lognormal") |>
    run_inference() |>
    expect_error("non-positive.*huey")

  # prior with negative lower bound
  new_housepi_model(model = "blank") |>
    add_parameter("miller", prior = "uniform",
                  min = -1, proposal = "lognormal") |>
    run_inference() |>
    expect_error("non-positive")

  # location parameter
  suppressWarnings(
    param_draws <- new_housepi_model(model = "blank") |>
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
    new_housepi_model("augmented_infection_time") |>
    get_infec_profile()
  expect_gt(length(out_1), 0L)

  flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    get_infec_profile() |>
    expect_length(0L)

  # reload first model
  sarscov2_2021_icofs |>
    new_housepi_model("augmented_infection_time") |>
    get_infec_profile() |>
    expect_identical(out_1)
})

test_that("Custom models shares function registries", {
  hhmodel_0 <- new_housepi_model(model = "blank")
  expose_function("growth") |>
    expect_error("Unknown function")

  # compile example
  custom <- system.file(file.path("examples", "custom"), package = "housepi")
  hhmodel_1 <- flu_1918_frederick |>
    new_housepi_model(custom, quiet = TRUE)

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
  hhmodel <- new_housepi_model(model = "blank") |>
    add_augmented("indiv_shape", prior = "normal")

  hhmodel$add_parameter("kaz", prior = "normal", mean = indiv_shape + 1) |>
    expect_error("[Ii]ndividual specific")

  val <- 1
  # UPDATE: for the moment I am forcing parameter arguments to be co
  hhmodel$add_parameter_list(
    param_shape = 5,
    var1 = list(prior = "normal", mean = param_shape@p + 1),
    var2 = list(prior = "normal", mean = param_shape@p + !!val),
    var3 = list(prior = "normal", mean = ~ param_shape@P + 1),
    var4 = list(prior = "normal", mean = ~ param_shape@P + !!val)
  )
  aux <- hhmodel$parameters
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 2L)

  val <- 3
  hhmodel$add_augmented_list(
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
    new_housepi_model("blank") |>
    add_augmented_list(
      hour = list(prior = "uniform", min = day, max = day + 1),
      day = list(prior = "gamma", shape = 2, scale = 1)
    ) |>
    run_simulation()
  expect_true(all(
    simulated$day <= simulated$hour & simulated$hour <= simulated$day + 1
  ))
})

test_that("Nested definitions works as expected", {
  hhmodel <- new_housepi_model(model = "blank") |>
    # define Box-Muller transformation
    add_augmented_list(norm = list(prior = "normal"),
                       u1 = list(prior = "uniform"),
                       u2 = list(prior = "uniform"),
                       z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2))
  expect_snapshot(print(hhmodel))

  results <- hhmodel$run_inference(steps = 1e5, chains = 1, adapt = 0,
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
    new_housepi_model("blank")

  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  # raise error for on injection
  hhmodel$add_augmented_list(external_list) |>
    expect_error("injection operator")

  tmp <- hhmodel$add_augmented_list(!!!external_list)$augmented
  hhmodel$add_augmented_list(z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})

test_that("Augmented traces export only relevant variables", {
  hhmodel <- new_housepi_model(model = "blank") |>
    add_parameter_list(alpha = 1, beta = list(prior = "uniform")) |>
    add_augmented_list(
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
  result <- suppressWarnings(hhmodel$run_inference(
    steps = 200, burn = 0, thin = 1, adapt = 0, chains = 1,
    id_draws = infec_df$id
  ))
  expect_equal(
    posterior::variables(result$augmented_draws),
    paste(infec_df$id, c("layer", "alias", "alpha_beta", "expl"), sep = "_")
  )
})

test_that("Augmented degenerate range warns and collapses", {
  hhmodel <- new_housepi_model(model = "blank") |>
    add_parameter("d_dog", prior = "uniform") |>
    add_augmented("d_horse", prior = "uniform", min = 5, max = 5)
  expect_warning(
    tryCatch(
      hhmodel$run_inference(steps = 10, chains = 1, adapt = 0, thin = 1),
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
  house <- data.frame(
    house_id = c(rep("hh1", 3), rep(hhid, 6)),
    status = c(1, 0, 0, 1, 0, 0, 0, 0, 0),
    member = c("student", "other", "other",
               "student", "sibling", "mother", "father", "other", "student")
  )
  hhmodel <- new_housepi_model(house, "chain_binomial")
  hhmodel$add_parameter_template()

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
  hhmodel$add_contact_rate(strata = strata, params = params)

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

  # drop_contact_rates cleans up the associated parameter definitions
  param_names_before <- hhmodel$parameters$name
  hhmodel$drop_contact_rates()
  any(setdiff(names(params), ".reference") %in% hhmodel$parameters$name) |>
    expect_false()
  # the template parameters added previously must survive
  template_pars <- c("alpha", "beta", "gamma")
  expect_true(all(template_pars %in% param_names_before))
})

test_that("Contact rates raise an error for incompatible models", {
  hhmodel <- data.frame(
    house_id = 1, id = seq(4), status = seq(4),
    cohort = c("a", "a", "b", "b")
  ) |>
    new_housepi_model("augmented_example")
  hhmodel$add_contact_rate(by = "cohort") |>
    expect_error("not supported")
})

test_that("run_simulation rejects wrong arguments", {
  hhmodel <- new_housepi_model(model = "blank")
  # wrong observation process
  hhmodel$run_simulation(obs_process = "not a function") |>
    expect_error("function")
  # wrong number of replicates
  hhmodel$run_simulation(nrep = 0) |>
    expect_error("positive integer")
  hhmodel$run_simulation(nrep = -3) |>
    expect_error("positive integer")
})

test_that("Param table with degenerate columns raises an ESS error", {
  hhmodel <- new_housepi_model(model = "blank")
  pt <- posterior::as_draws_df(data.frame(p = rep(1, 100)))
  hhmodel$run_simulation(param_table = pt) |>
    expect_error("ESS")
})

test_that("run_simulation auto-imports parameters from param_table", {
  hhmodel <- data.frame(house_id = 1, id = letters[1:3], status = 0) |>
    new_housepi_model("blank")
  expect_equal(nrow(hhmodel$parameters), 0)
  pt <- data.frame(adamska = c(0.1, 0.2, 0.3))
  invisible(hhmodel$run_simulation(param_table = pt))
  expect_true("adamska" %in% hhmodel$parameters$name)
})

test_that("run_measure_ rejects apply_long without obs_process_long", {
  hhmodel <- new_housepi_model(model = "blank")
  hhmodel$add_parameter_list(alpha = 1)
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
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  hhmodel$add_parameter_template()
  # supplying obs_process_long forces the long branch into a non-empty list
  # single simulation by default
  obs_long <- function(x) data.frame(id = "vic-boss", measure = 1.0)
  out <- hhmodel$run_simulation(obs_process = function(x) x,
                                obs_process_long = obs_long)
  expect_true(is.list(out))
  expect_equal(names(out), c("wide", "long"))
  expect_true(is.data.frame(out$wide))
  expect_true(is.data.frame(out$long))

  # multiple repetitions
  out <- hhmodel$run_simulation(obs_process = function(x) x,
                                obs_process_long = obs_long,
                                nrep = 2)
  expect_true(is.data.frame(out$wide[[1]]))
  expect_true(is.data.frame(out$long[[2]]))
})

test_that("add_augmented_template loads from a non-blank model template", {
  hhmodel <- new_housepi_model(
    sarscov2_2021_icofs,
    "augmented_infection_time"
  )
  hhmodel$add_augmented_template()
  expect_gt(nrow(hhmodel$augmented), 0L)
})

test_that("add_infec_profile and list update the model directly", {
  hhmodel <- sarscov2_2021_icofs |>
    new_housepi_model("augmented_infection_time") |>
    add_infec_profile(func = "gamma", shape = 2, scale = 1)
  expect_true("gamma" %in% hhmodel$infec_profile$func)
  expect_warning(hhmodel$add_infec_profile(func = "lognormal"),
                 regexp = "Overwriting")
  expect_true("lognormal" %in% hhmodel$infec_profile$func)
})

test_that("Parameter defined as a function tracks its dependencies", {
  set.seed(1)
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_list(
      alpha = 1e-4,
      beta = size_ref@P^-gamma@P,
      gamma = list(prior = "uniform"),
      size_ref = 4
    )
  res <- suppressWarnings(hhmodel$run_inference(
    steps = 100, chains = 1, adapt = 0, thin = 1
  ))
  vars <- posterior::variables(res$parameter_draws)
  # both sampled and function-derived parameters appear in the trace
  expect_true(all(c("gamma", "beta") %in% vars))
  # the function-derived parameter has no prior contribution
  expect_false("lprior_beta" %in% posterior::variables(
    res$.__enclos_env__$private$param_prior_draws_))

  df <- as.data.frame(posterior::as_draws_df(res$parameter_draws))
  # every draw of beta equals size_ref^(-gamma) exactly
  expect_equal(df$beta, 4^(-df$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")

  # bare names default to augmented source! here @P is missing
  hhmodel$add_parameter_list(alpha = 1e-4, gamma = 1, size_ref = 4,
                             beta = size_ref^-gamma) |>
    expect_error("cannot depend on augmented variables")

  # explicit reference to a data column is also rejected
  hhmodel$add_parameter_list(alpha = 1e-4, gamma = 1, size_ref = 4,
                             beta = status@D + 1) |>
    expect_error("cannot depend on data variables")
})

test_that("Parameter defined as a function tracks its dependencies", {
  set.seed(1)
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  hhmodel$add_parameter_list(alpha = 1e-4,
                             beta = size_ref@P^-gamma@P,
                             gamma = list(prior = "uniform"),
                             size_ref = 4)
  res <- hhmodel$run_inference(steps = 2000, progress = 0)
  vars <- posterior::variables(res$parameter_draws)
  # both sampled and function-derived parameters appear in the trace
  expect_true(all(c("gamma", "beta") %in% vars))
  # the function-derived parameter has no prior contribution
  expect_false("lprior_beta" %in% posterior::variables(
    res$.__enclos_env__$private$param_prior_draws_))

  df <- as.data.frame(posterior::as_draws_df(res$parameter_draws))
  # snake on the trail: every draw of beta equals size_ref^(-gamma) exactly
  expect_equal(df$beta, 4^(-df$gamma))
})

test_that("Parameter function rejects non-parameter dependencies", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial")

  # bare names default to augmented source - rejected by C++ validator
  hhmodel$add_parameter_list(alpha = 1e-4,
                             beta = size_ref^-gamma,
                             gamma = list(prior = "uniform"),
                             size_ref = 4) |>
    expect_error("cannot depend on augmented variables")

  # explicit reference to a data column is also rejected
  hhmodel$add_parameter_list(alpha = 1e-4,
                             beta = list(prior = "uniform"),
                             gamma = status@D + 1) |>
    expect_error("cannot depend on data variables")
})

