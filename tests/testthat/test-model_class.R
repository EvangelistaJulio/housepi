library(testthat)

test_that("Model class inconsistent calls", {
  # Default blank
  expect_warning(new_housepi_model(), regexp = "Falling back to")

  # blank household model object
  hhmodel <- new_housepi_model("blank")

  # inference without definitions
  expect_error(hhmodel$run_inference(), regexp = "Nothing to sample")
})

test_that("Model class blank calls", {
  # blank household model object
  hhmodel <- new_housepi_model("blank")

  # blank default data wide format
  data_wide <- housepi:::template_blank_data_wide_()
  expect_equal(hhmodel$data_wide, data_wide)

  # correct data long format
  data_long <- data.frame(id = data_wide$id[1],
                          time = seq(3), measure = seq(3)**3)
  hhmodel$data_long <- data_long
  expect_equal(hhmodel$data_long, data_long)

  # factorization
  complex_data <- housepi::survey
  complex_data$dates <- as.Date(runif(dim(complex_data)[1], max = 100))
  hhmodel$data_wide <- complex_data
  expect_equal(hhmodel$data_wide, complex_data)

  # templates
  expect_warning(hhmodel$load_parameter_template(), regexp = "blank")
  expect_warning(hhmodel$load_augmented_template(), regexp = "blank")

  # zero log-likelihood
  expect_warning(x <- hhmodel$evaluate_log_lik(), regexp = "blank")
  expect_equal(x, 0)

  # blank simulation
  hhmodel <- new_housepi_model("blank")
  hhmodel$add_augmented(name = "dummy", value = 3)
  expect_warning(out <- hhmodel$run_simulation(), regexp = "blank")
  expect_true(
    all(sapply(out$dummy[out$status == 0], function(x) identical(x, NA_real_)))
  )
})

test_that("Active bindings work as expected", {
  hhmodel <- new_housepi_model("blank")
  
  # parameters
  hhmodel$parameters <- list(sample = list(prior = "uniform"),
                             uniform = list(prior = "uniform"))
  defs_a <- hhmodel$parameters
  hhmodel$load_parameter_list(sample = list(prior = "uniform"),
                              uniform = list(prior = "uniform"))
  expect_equal(defs_a, hhmodel$parameters)
  
  val <- 1
  hhmodel$load_parameter_list(
   var1 = list(prior = "normal", mean = indiv_shape + 1),
   var2 = list(prior = "normal", mean = indiv_shape + !!val),
   var3 = list(prior = "normal", mean = ~ indiv_shape + 1),
   var4 = list(prior = "normal", mean = ~ indiv_shape + !!val)
  )
  aux <- hhmodel$parameters
  aux$name <- NULL  
  expect_equal(nrow(unique(aux)), 1)
  
  # augmented
  hhmodel$augmented <- list(sample = list(model = "uniform"),
                            uniform = list(model = "uniform"))
  defs_a <- hhmodel$augmented
  hhmodel$load_augmented_list(sample = list(model = "uniform"),
                              uniform = list(model = "uniform"))
  expect_equal(defs_a, hhmodel$augmented)
  
  val <- 3
  hhmodel$load_augmented_list(
   var1 = list(model = "normal", mean = indiv_shape + 3),
   var2 = list(model = "normal", mean = indiv_shape + !!val),
   var3 = list(model = "normal", mean = ~ indiv_shape + 3),
   var4 = list(model = "normal", mean = ~ indiv_shape + !!val)
  )
  aux <- hhmodel$augmented
  aux$name <- NULL
  expect_equal(nrow(unique(aux)), 1)
})

test_that("Discretization is enabled as expected", {
  # blank household model object
  hhmodel <- new_housepi_model("blank")
  hhmodel$load_parameter_list(
    sample = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE),
    uniform = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                   proposal = "uniform"),
    normal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                  proposal = "normal"),
    lognormal = list(prior = "gamma", shape = 2, scale = 1, discr = TRUE,
                     proposal = "lognormal")
  )

  results <- hhmodel$run_inference()
  y <- posterior::as_draws_matrix(get_param_draws(results))
  expect_true(all(y == floor(y)))
})

test_that("Augmented aliasing works as expected", {
  # check that model arguments are correctly integrated
  hhmodel <- new_housepi_model("blank") |>
    add_data_wide(data.frame(id = seq(1, 1000), house_id = 1, status = 1)) |>
    load_augmented_list(
      hour = list(model = "uniform", min = day, max = day + 1),
      day = list(model = "gamma", shape = 2, scale = 1)
    )
  simulated <- run_simulation(hhmodel)
  expect_true(all(simulated$day <= simulated$hour &
                    simulated$hour <= simulated$day + 1))
})

test_that("Nested definittions works as expected", {
  hhmodel <- new_housepi_model("blank")

  # define Box-Muller transformation
  hhmodel$load_augmented_list(
    norm = list(model = "normal"),
    u1 = list(model = "uniform"),
    u2 = list(model = "uniform"),
    z0 = sqrt(- 2 * log(u1)) * cos(6.283185 * u2)
  )

  results <- hhmodel$run_inference(steps = 1e5, progress = 0)
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
  hhmodel <- new_housepi_model("blank")
  
  external_list <- list(z0 = ~ date@D, z1 = ~ sqrt(log(z0)))
  # raise error for on injection
  expect_error(hhmodel$load_augmented_list(external_list),
               regexp = "injection operator")
  
  tmp <- hhmodel$load_augmented_list(!!!external_list)$augmented
  hhmodel$load_augmented_list(z0 = date@D, z1 = sqrt(log(z0)))
  expect_equal(tmp, hhmodel$augmented)
})
