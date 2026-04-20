library(testthat)

##################################
## MODEL CLASS : STANDARD CALLS ##
##################################

test_that("Model class inconsistent calls", {
  # Default blank
  expect_warning(new_housepi_model(), regexp = "Falling back to")

  # blank household model object
  hhmodel <- new_housepi_model(model = "blank")

  # inference without definitions
  expect_error(hhmodel$run_inference(), regexp = "Nothing to sample")

  # data_wide type errors surface at construction
  expect_error(new_housepi_model(data.frame(), "blank"),
               regexp = "empty")
  expect_error(new_housepi_model("ground", "blank"),
               regexp = "data.frame")

  # data_long errors still surface via add_data_long
  expect_error(hhmodel$add_data_long("zero"), regexp = "data.frame")

  # non-blank models require data_wide at construction
  expect_error(new_housepi_model(model = "chain_binomial"),
               regexp = "data")
})

test_that("Active fields are read-only", {
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  expect_error(hhmodel$data_wide <- "La", regexp = "read only")
  expect_error(hhmodel$data_long <- "Li", regexp = "read only")
  expect_error(hhmodel$parameters <- "Lu", regexp = "read only")
  expect_error(hhmodel$augmented <- "Le", regexp = "read only")
  expect_error(hhmodel$infec_profile <- "Lo", regexp = "read only")
  expect_error(hhmodel$group_rates <- "The", regexp = "read only")
  expect_error(hhmodel$name <- "Patriots", regexp = "read only")
  expect_error(hhmodel$setup <- "Cipher", regexp = "read only")
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
  data_long <- data.frame(id = snapshot$data_wide$id[1],
                          time = seq(3), measure = seq(3)**3)
  hhmodel$add_data_long(data_long)
  expect_equal(hhmodel$data_long, data_long)

  # factorization round-trips through a fresh blank model
  complex_data <- sarscov2_2021_icofs
  complex_data$dates <- as.Date(runif(dim(complex_data)[1], max = 100))
  complex_model <- new_housepi_model(complex_data, "blank")
  expect_equal(complex_model$data_wide, complex_data)

  # templates
  expect_warning(hhmodel$load_parameter_template(), regexp = "blank")
  expect_warning(hhmodel$load_augmented_template(), regexp = "blank")

  # zero log-likelihood
  expect_warning(x <- hhmodel$evaluate_log_lik(), regexp = "blank")
  expect_equal(x, 0)

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
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  hhmodel$load_parameter_template()
  expect_snapshot(print(hhmodel))

  # without parameters
  hhmodel2 <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  expect_snapshot(print(hhmodel2))

  # blank model
  expect_snapshot(print(new_housepi_model(model = "blank")))
})

##################################
## MODEL CLASS : DATA SUMMARIES ##
##################################

test_that("Flu 2014 Matsumoto is correctly pre-processed", {
  hhmodel <- flu_2014_matsumoto |>
    new_housepi_model("chain_binomial_grouped")

  case_summary <- summarise_cases(hhmodel)
  expect_true(is.data.frame(case_summary))
  # expect same household count as in Endo, et al. ('19)
  expect_equal(
    subset(case_summary, type == "household" & comment == "total")$counts,
    10486
  )
  # house_size vs status equivalence groups correspond to measure non-zeros
  aux <- measure_cases_by_house_size(flu_2014_matsumoto)
  aux$house_size <- NULL
  expect_equal(
    subset(case_summary, type == "household" & comment == "distinct")$counts,
    sum(aux[!is.na(aux)] > 0)
  )

  # fix expected equivalence classes by member composition
  case_summary <- hhmodel |>
    add_group_rate(by = "member") |>
    summarise_cases()
  expect_equal(
    subset(case_summary, type == "household" & comment == "distinct")$counts,
    454
  )

  expect_true(is.data.frame(summarise_house_sizes(hhmodel)))
})

######################################
## MODEL CLASS : POSTERIOR SAMPLING ##
######################################

test_that("Wrong definitions throw errors", {
  hhmodel <- new_housepi_model(model = "blank")

  # eager rebuild: errors surface at definition, not at run_inference
  expect_error(
    hhmodel$load_parameter_list(a = list(prior = "uniform", min = ~min_val@D)),
    regexp = "not supported"
  )
  # hint
  expect_error(
    hhmodel$load_parameter_list(b = list(prior = "uniform", min = min_val)),
    regexp = "rlang injection"
  )
})

test_that("Discretization is enabled as expected", {
  discr <- TRUE
  val <- list(min = 1, max = 5)
  hhmodel <- new_housepi_model(model = "blank")
  hhmodel$load_parameter_list(
    sample = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                  proposal = "sample", min = !!val$min, max = !!val$max),
    uniform = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                   proposal = "uniform", min = !!val$min, max = !!val$max),
    normal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                  proposal = "normal", min = !!val$min, max = !!val$max),
    lognormal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                     proposal = "lognormal", min = !!val$min, max = !!val$max)
  )

  results <- hhmodel$run_inference()
  y <- posterior::as_draws_matrix(get_param_draws(results))
  expect_true(all(y == floor(y)))
  expect_true(all(y >= val$min))
  expect_true(all(y <= val$max))
})

#########################################
## MODEL CLASS : AUGMENTED DEFINITIONS ##
#########################################

test_that("Modules do not get overwritten", {
  # NOTE: I use profile as a  status flag for models
  hhmodel_1 <- new_housepi_model(sarscov2_2021_icofs,
                                 "augmented_infection_time")
  out_1 <- hhmodel_1$infec_profile
  expect_gt(length(out_1), 0)

  hhmodel_2 <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  out_2 <- hhmodel_2$infec_profile
  expect_length(out_2, 0)

  # reload first model
  hhmodel_3 <- new_housepi_model(sarscov2_2021_icofs,
                                 "augmented_infection_time")
  out_3 <- hhmodel_3$infec_profile
  expect_identical(out_1, out_3)
})

test_that("Custom models shares function registries", {
  hhmodel_0 <- new_housepi_model(model = "blank")
  expect_error(expose_function("growth"), regexp = "Unknown function")

  # compile example
  custom <- system.file(file.path("examples", "custom"), package = "housepi")
  hhmodel_1 <- flu_1918_frederick |>
    new_housepi_model(custom, quiet = TRUE)
  expect_no_error(expose_function("growth")$func(seq(3)))
})

#########################################
## MODEL CLASS : AUGMENTED DEFINITIONS ##
#########################################

test_that("Rlang calls work as expected", {
  # declare indiv_shape as an augmented variable so references resolve
  hhmodel <- new_housepi_model(model = "blank") |>
    add_augmented("indiv_shape", model = "normal")

  expect_error(
    hhmodel$load_parameter_list(
      var1 = list(prior = "normal", mean = indiv_shape + 1)
    ),
    regexp = "[Ii]ndividual specific"
  )

  val <- 1
  # UPDATE: for the moment I am forcing parameter arguments to be co
  hhmodel$load_parameter_list(
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
  hhmodel$load_augmented_list(
    var1 = list(model = "normal", mean = indiv_shape + 3),
    var2 = list(model = "normal", mean = indiv_shape + !!val),
    var3 = list(model = "normal", mean = ~ indiv_shape + 3),
    var4 = list(model = "normal", mean = ~ indiv_shape + !!val),
    append = TRUE
  )
  aux <- hhmodel$augmented
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 2L)
})

test_that("Augmented aliasing works as expected", {
  # check that model arguments are correctly integrated
  hhmodel <- data.frame(id = seq(1, 1000), house_id = 1, status = 1) |>
    new_housepi_model("blank") |>
    load_augmented_list(
      hour = list(model = "uniform", min = day, max = day + 1),
      day = list(model = "gamma", shape = 2, scale = 1)
    )
  simulated <- run_simulation(hhmodel)
  expect_true(all(simulated$day <= simulated$hour &
                    simulated$hour <= simulated$day + 1))
})

test_that("Nested definitions works as expected", {
  hhmodel <- new_housepi_model(model = "blank")

  # define Box-Muller transformation
  hhmodel$load_augmented_list(
    norm = list(model = "normal"),
    u1 = list(model = "uniform"),
    u2 = list(model = "uniform"),
    z0 = sqrt(- 2 * log(u1)) * cos(2 * pi * u2)
  )
  results <- hhmodel$run_inference(steps = 1e5, progress = 0)
  # no parameters nor priors
  empty_draws <- posterior::as_draws_df(data.frame())
  expect_identical(results$param_draws, empty_draws)
  expect_identical(results$prior_draws, empty_draws)
  expect_snapshot(print(hhmodel))

  sample_norm <- results$augm_draws$tommy_norm
  sample_boxmul <- results$augm_draws$tommy_z0

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
  hhmodel <- data.frame(id = "id1", house_id = "h1", status = 1,
                        date = as.Date("2021-01-01")) |>
    new_housepi_model("blank")

  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  # raise error for on injection
  expect_error(hhmodel$load_augmented_list(external_list),
               regexp = "injection operator")

  tmp <- hhmodel$load_augmented_list(!!!external_list)$augmented
  hhmodel$load_augmented_list(z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})

test_that("Augmented traces export only relevant variables", {
  hhmodel <- new_housepi_model(model = "blank")
  infec_df <- subset(hhmodel$data_wide, status > 0)
  hhmodel$load_parameter_list(alpha = 1, beta = list(prior = "uniform"))
  hhmodel$load_augmented_list(
    # constants
    blocked = 42,
    # private variables
    .hidden = list(model = "uniform"),
    .secret = list(model = "normal"),
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
  result <- hhmodel$run_inference(id_draws = infec_df$id)
  expect_equal(
    posterior::variables(result$augm_draws),
    paste(infec_df$id, c("layer", "alias", "alpha_beta", "expl"), sep = "_")
  )
})

############################################
## MODEL CLASS : CONTACT RATE DEFINITIONS ##
############################################

test_that("Contact rate definition work as expected", {
  hhmodel <- new_housepi_model(flu_2014_matsumoto, "chain_binomial_grouped")
  hhmodel$load_parameter_template()

  # NOTE: larger households at the end
  hhid <- flu_2014_matsumoto$house_id[length(flu_2014_matsumoto$house_id)]

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
})
