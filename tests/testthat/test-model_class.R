library(testthat)

test_that("Model class inconsistent calls", {
  # Default blank
  expect_warning(new_housepi_model(), regexp = "Falling back to")

  # blank household model object
  hhmodel <- new_housepi_model("blank")

  # inference without definitions
  expect_error(hhmodel$infer(), regexp = "Nothing to sample")
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

  # zero loglikelihood
  expect_warning(x <- hhmodel$loglikelihood(), regexp = "blank")
  expect_equal(x, 0)

  # blank simulation
  hhmodel <- new_housepi_model("blank")
  hhmodel$add_augmented(name = "dummy", value = 3)
  expect_warning(out <- hhmodel$simulate(), regexp = "blank")
  expect_true(
    all(sapply(out$dummy[out$status == 0], function(x) identical(x, NA_real_)))
  )
})

test_that("Discretization is enabled as expected", {
  # blank household model object
  hhmodel <- new_housepi_model("blank")
  hhmodel$load_parameter_list(
    sample = list(prior = "gamma", discr = TRUE),
    uniform = list(prior = "gamma", discr = TRUE, rule = "uniform"),
    normal = list(prior = "gamma", discr = TRUE, rule = "normal"),
    lognormal = list(prior = "gamma", discr = TRUE, rule = "lognormal")
  )

  y <- posterior::as_draws_matrix(get_param_draws(hhmodel$infer()))
  expect_true(all(y == floor(y)))
})

test_that("Nested definitions raise errors", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$data_wide <- data.frame(house_id = 1, id = seq(4), status = seq(4))

  expect_no_error(
    hhmodel$load_augmented_list(
      ..layer = list(model = "uniform", status = 1),
      ..layer = list(expr = ~..dummy, status = c(2, 3, 4)),
      ..dummy = list(model = "uniform"),
      proba = ~..layer
    )$loglikelihood()
  )

  # simple nest
  expect_error(
    hhmodel$load_augmented_list(..dummy = ~proba,
                                proba = ~..dummy)$loglikelihood(),
    regexp = "nested definitions"
  )

  # triangular nest
  expect_error(
    hhmodel$load_augmented_list(proba = ~..layer,
                                ..layer = ~..dummy,
                                ..dummy = ~proba)$loglikelihood(),
    regexp = "nested definitions"
  )

  # first generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel$load_augmented_list(
      ..layer = list(model = "normal", status = 1),
      ..layer = list(expr = ~..dummy, status = 2),
      #> NO ..layer provided for status 3 and 4
      ..dummy = list(model = "lognormal", status = 2),
      proba = ~..layer
    )$loglikelihood()
  ))
  expect_match(err_msg, "Could not find a definition")
  expect_match(err_msg, "status 3")

  # second generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel$load_augmented_list(
      ..layer = list(model = "normal", status = 1),
      ..layer = list(expr = ~..dummy, status = c(2, 3, 4)),
      #> NO ..dummy provided for status 3, 4
      ..dummy = list(model = "lognormal", status = c(2, 3)),
      proba = ~..layer
    )$loglikelihood()
  ))
  expect_match(err_msg, "Could not find a definition")
  expect_match(err_msg, "status 4")
})

test_that("Missing function arguments raise errors", {
  hhmodel <- new_housepi_model("augmented_example")
  # data with NA values
  hhmodel$data_wide <- data.frame(house_id = seq(3), id = seq(3),
                                  status = 1, test = c(5.5, NA, -1.5))

  hhmodel$add_augmented(name = "proba", func = "gamma")
  expect_error(hhmodel$loglikelihood(), regexp = "Impossible to return value")

  hhmodel$add_augmented(name = "gamma", func = "gamma")
  expect_warning(
    hhmodel$add_augmented(name = "proba", expr = ~exp(-gamma)),
    "Overwriting existing"
  )
  expect_error(hhmodel$loglikelihood(), regexp = "Impossible to associate")

  # check that NaN + something causes errors
  hhmodel$load_augmented_list(noise = list(model = "normal"),
                              proba = ~exp(test@D + noise))
  expect_error(hhmodel$loglikelihood(), regexp = "Impossible to associate")

  # check that log(NaN) = NaN
  hhmodel$load_augmented_list(proba = ~log(test@D))
  expect_error(hhmodel$infer(), regexp = "Nothing to sample")
  expect_error(hhmodel$loglikelihood(), regexp = "Impossible to return value")
})

test_that("Augmented loglikelihood seeding is fine", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$data_wide <- data.frame(house_id = 1, id = seq(4), status = seq(4))
  hhmodel$load_augmented_list(
    proba = list(model = "uniform", status = 1),
    proba = list(expr = ~..layer, status = c(2, 3, 4)),
    ..layer = list(model = "uniform")
  )

  ll1 <- hhmodel$loglikelihood(seed = 1)
  ll2 <- hhmodel$loglikelihood(seed = 2)
  expect_true(ll1 != ll2)
  expect_identical(ll2, hhmodel$loglikelihood(seed = 2))
  expect_identical(ll2, sum(log(hhmodel$export_augmented_data()$proba)))
})

test_that("Augmented aliasing works as expected", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$data_wide <- data.frame(
    house_id = "hh1", id = seq(1, 3), x = 2 + runif(3), status = 1
  )
  hhmodel$add_parameter(name = "scale", value = 1)
  hhmodel$load_augmented_list(
    ..dummy = ~x @D,
    ..gamma = list(func = "gamma",
                   func_args = c(x = ~..dummy, shape = 2, scale = ~ scale@P)),
    proba = ~..gamma
  )

  ll <- sum(log(dgamma(hhmodel$data_wide$x, shape = 2, scale = 1)))
  expect_equal(hhmodel$loglikelihood(), ll)

  # check that order does not matter
  hhmodel$load_augmented_list(
    ..gamma = list(func = "gamma",
                   func_args = c(x = ~..dummy, shape = 2, scale = ~ scale@P)),
    proba = ~..gamma,
    ..dummy = ~x @D
  )
  expect_equal(hhmodel$loglikelihood(), ll)
  
  # check normalization of binomial with augmented example
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$data_wide <- data.frame(
    id = seq(0, 100), house_id = "hh1", status = 1, X = seq(0, 100))
  hhmodel$load_augmented_list(
    pdf = list(func = "binomial", 
               func_args = c(x = ~ X@D, size = 10, prob = 0.5)),
    proba = ~exp(pdf)
  )
  expect_equal(hhmodel$loglikelihood(), 1)
  
  # check that model arguments are correctly integrated
  hhmodel <- new_housepi_model("blank") |> 
    add_data_wide(data.frame(id = seq(1, 1000), house_id = 1, status = 1)) |> 
    load_augmented_list(
      delay = list(model = "uniform", args = c(min = ~ day, max = ~ day + 1)),
      day = list(model = "gamma", args = c(shape = 2, scale = 1))
    )
  expect_warning(simulated <- simulate(hhmodel), regexp = "blank")
  expect_true(all(simulated$day < simulated$delay))
})

test_that("Augmentation example on incubation from delays is coherent", {
  # loglikelihood from simulated delays
  shape <- 2.0
  scale <- 1.1
  df <- data.frame(status = 1,
                   ip = rgamma(223, shape = shape, scale = scale))
  exact <- sum(log(stats::dgamma(df$ip, shape = shape, scale = scale)))

  # generate an instance of the household model class
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$load_augmented_list(
    proba = list(func = "gamma",
                 func_args = c(x = ~ip@D, shape = shape, scale = scale))
  )
  # warning for no structure
  expect_warning(hhmodel$data_wide <- df,
                 regexp = "Ignored household structure")
  expect_equal(hhmodel$loglikelihood(), exact)
})

test_that("Nested definittions works as expected", {
  hhmodel <- new_housepi_model("blank")

  # define Box-Muller transformation
  hhmodel$load_augmented_list(
    norm = list(model = "normal"),
    u1 = list(model = "uniform"),
    u2 = list(model = "uniform"),
    z0 = ~sqrt(- 2 * log(u1)) * cos(6.283185 * u2)
  )

  results <- hhmodel$infer(steps = 1e5, progress = 0)
  sample_norm <- results$augm_draws$hh_01_1_norm
  sample_boxmul <- results$augm_draws$hh_01_1_z0

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
