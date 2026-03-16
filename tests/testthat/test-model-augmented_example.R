library(testthat)

test_that("Nested definitions raise errors", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$add_data_wide(
    data.frame(house_id = 1, id = seq(4), status = seq(4))
  )

  expect_no_error(
    hhmodel$load_augmented_list(
      .layer = list(model = "uniform", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      .dummy = list(model = "uniform"),
      proba = .layer
    )$evaluate_log_lik()
  )

  # simple nest
  expect_error(
    hhmodel$load_augmented_list(.dummy = proba,
                                proba = .dummy)$evaluate_log_lik(),
    regexp = "cyclic definitions"
  )

  # triangular nest
  expect_error(
    hhmodel$load_augmented_list(proba = .layer,
                                .layer = .dummy,
                                .dummy = proba)$evaluate_log_lik(),
    regexp = "cyclic definitions"
  )

  # first generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel$load_augmented_list(
      .layer = list(model = "normal", status = 1),
      .layer = list(value = .dummy, status = 2),
      #> NO .layer provided for status 3 and 4
      .dummy = list(model = "lognormal", status = 2),
      proba = .layer
    )$evaluate_log_lik()
  ))
  expect_match(err_msg, "Could not find a definition")
  expect_match(err_msg, "status")

  # second generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel$load_augmented_list(
      .layer = list(model = "normal", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      #> NO .dummy provided for status 3, 4
      .dummy = list(model = "lognormal", status = c(2, 3)),
      proba = .layer
    )$evaluate_log_lik()
  ))
  expect_match(err_msg, "Could not find a definition")
  expect_match(err_msg, "status")
})

test_that("Missing function arguments raise errors", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$add_data_wide(data.frame(house_id = seq(3), id = seq(3),
                                   status = 1, test = c(5.5, NA, -1.5)))

  # non default arguments
  expect_error(
    hhmodel$add_augmented(name = "proba", func = "gamma"),
    regexp = "Missing non default"
  )

  # data with NA values
  hhmodel$drop_augmented()
  hhmodel$add_augmented(name = "proba", func = "gamma",
                        func_args = c(shape = 2))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to return")

  # concatenated data with NA values
  hhmodel$load_augmented_list(
    gamma = list(func = "gamma",  func_args = c(shape = 2)),
    proba =  exp(-gamma)
  )
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")

  # check that NaN + something causes errors
  hhmodel$load_augmented_list(noise = list(model = "normal"),
                              proba = exp(test@D + noise))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")

  # check that log(NaN) causes error
  hhmodel$load_augmented_list(proba = log(test@D))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")

  # check that log(0) causes error
  hhmodel$load_augmented_list(noise = 0., proba = log(noise))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")
})

test_that("Augmented log-likelihood seeding is fine", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$add_data_wide(
    data.frame(house_id = 1, id = seq(4), status = seq(4))
  )
  hhmodel$load_augmented_list(
    proba = list(model = "uniform", status = 1),
    proba = list(value = .layer, status = c(2, 3, 4)),
    .layer = list(model = "uniform")
  )

  ll1 <- hhmodel$evaluate_log_lik(seed = 1)
  ll2 <- hhmodel$evaluate_log_lik(seed = 2)
  expect_true(ll1 != ll2)
  expect_equal(ll2, hhmodel$evaluate_log_lik(seed = 2))
  expect_equal(ll2, sum(log(hhmodel$export_augmented_data()$proba)))
})

test_that("Augmented aliasing works as expected", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$add_data_wide(data.frame(
    house_id = "hh1", id = seq(1, 3), x = 2 + runif(3), status = 1
  ))
  hhmodel$add_parameter(name = "scale", value = 1)
  hhmodel$load_augmented_list(
    .dummy = x@D,
    .layer = list(func = "gamma",
                  func_args = c(x = .dummy, shape = 2, scale = scale@P)),
    proba = .layer
  )

  ll <- sum(log(dgamma(hhmodel$data_wide$x, shape = 2, scale = 1)))
  expect_equal(hhmodel$evaluate_log_lik(), ll)

  # check that order does not matter
  hhmodel$load_augmented_list(
    .layer = list(func = "gamma",
                  func_args = c(x = .dummy, shape = 2, scale = scale@P)),
    proba = .layer,
    .dummy = x@D
  )
  expect_equal(hhmodel$evaluate_log_lik(), ll)

  # check normalization of binomial with augmented example
  size <- 10
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$add_data_wide(data.frame(
    id = seq(0, size), house_id = "hh1", status = 1, X = seq(0, size)
  ))
  hhmodel$load_augmented_list(
    pdf = list(func = "binomial",
               func_args = c(x = X@D, size = !!size, prob = 0.5)),
    proba = exp(pdf)
  )
  expect_equal(hhmodel$evaluate_log_lik(), 1)
})

test_that("Augmentation example on incubation from delays is coherent", {
  shape <- 1.95
  scale <- 1.13
  exact_val <- function(x, shape, scale) {
    sum(log(stats::dgamma(x, shape = shape, scale = scale)))
  }

  # log-likelihood from simulated delays

  df <- data.frame(status = 1,
                   ip = rgamma(223, shape = shape, scale = scale))

  # generate an instance of the household model class
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$load_augmented_list(
    proba = list(func = "gamma",
                 func_args = c(x = ip@D, shape = !!shape, scale = !!scale))
  )
  # warning for no structure
  expect_warning(hhmodel$add_data_wide(df),
                 regexp = "Ignored household structure")
  expect_equal(hhmodel$evaluate_log_lik(), exact_val(df$ip, shape, scale))

  # repeat with degenerated ip values to trigger multiplicity
  df$ip <- floor(df$ip)
  hhmodel$add_data_wide(df)
  expect_equal(hhmodel$evaluate_log_lik(), exact_val(df$ip, shape, scale))
})

test_that("Known errors are unchanged", {
  # Known infinity (values below minimum, i. e. proba 0)
  hhmodel <- new_housepi_model("augmented_example") |>
    add_data_wide(housepi::survey) |>
    load_augmented_list(
      value = list(model = "uniform", min = 0, max = 1),
      proba = list(func = "lognormal",
                   func_args = c(x = value, meanlog = 1, sdlog = 1),
                   min = 1, max = 15, discr = TRUE)
    ) |>
    evaluate_log_lik() |>
    expect_equal(-Inf)

  # Known NaN: qlnorm(15, meanlog = 4, sdlog = 0.01) = 0, impossible truncation
  # UPDATE: this should raise a warning instead
  new_housepi_model("augmented_example") |>
    add_data_wide(housepi::survey) |>
    load_augmented_list(
      value = list(model = "uniform", min = 1, max = 10, discr = TRUE),
      proba = list(func = "lognormal",
                   func_args = c(x = value, meanlog = 4, sdlog = 0.01),
                   min = 1, max = 15, discr = TRUE)
    ) |>
    evaluate_log_lik() |>
    expect_error(regexp = "returned NaN")
})

test_that("Augmented model log-likelihoods are correctly computed", {
  # dummy dataset
  house <- data.frame(id = sprintf("hh01m%s", seq(1, 7)),
                      house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      x = c(1.1, 2.2, 3.3, NA, NA, NA, NA))
  shape <- 2
  scale <- 1.1
  id_draws <- house[house$status > 0, ]$id
  results <- new_housepi_model("augmented_example") |>
    add_data_wide(house) |>
    load_augmented_list(proba = list(model = "gamma", x = x@D,
                                     shape = !!shape, scale = !!scale)) |>
    run_inference(steps = 5, thin = 1, adapt = 0, chain = 1,
                  progress = FALSE, id_draws = id_draws)

  # coherency between pointwise and total log-likelihood
  expect_equal(get_total_log_lik_draws(results)$log_lik,
               get_log_lik_draws(results)$hh01)

  # coherency between augmented draws and log-likelihood
  df <- as.data.frame(get_augm_draws(results))
  df <- as.matrix(df[, !startsWith(names(df), "."), drop = FALSE])
  expect_equal(get_log_lik_draws(results)$hh01,
               rowSums(log(df * dgamma(df, shape = shape, scale = scale))))
})
