library(testthat)

test_that("Nested definitions raise errors", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("augmented_example")

  expect_no_error(
    hhmodel |>
      add_augmented_list(
        .layer = list(model = "uniform", status = 1),
        .layer = list(value = .dummy, status = c(2, 3, 4)),
        .dummy = list(model = "uniform"),
        proba = .layer
      ) |>
      evaluate_log_lik()
  )

  # simple nest
  expect_error(
    hhmodel |>
      add_augmented_list(.dummy = proba, proba = .dummy) |>
      evaluate_log_lik(),
    regexp = "cyclic definition"
  )

  # triangular nest
  expect_error(
    hhmodel |>
      add_augmented_list(proba = .layer, .layer = .dummy, .dummy = proba) |>
      evaluate_log_lik(),
    regexp = "cyclic definition"
  )

  # first generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel$add_augmented_list(
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
    hhmodel$add_augmented_list(
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
  hhmodel <- data.frame(house_id = seq(3), id = seq(3),
                        status = 1, test = c(5.5, NA, -1.5)) |>
    new_housepi_model("augmented_example")

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
  hhmodel$add_augmented_list(
    gamma = list(func = "gamma",  func_args = c(shape = 2)),
    proba =  exp(-gamma)
  )
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")

  # check that NaN + something causes errors
  hhmodel$add_augmented_list(noise = list(model = "normal"),
                             proba = exp(test@D + noise))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")

  # check that log(NaN) causes error
  hhmodel$add_augmented_list(proba = log(test@D))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")

  # check that log(0) causes error
  hhmodel$add_augmented_list(noise = 0., proba = log(noise))
  expect_error(hhmodel$evaluate_log_lik(), regexp = "Impossible to associate")
})

test_that("Augmented log-likelihood seeding is fine", {
  hhmodel <- data.frame(house_id = 1, id = seq(4),
                        status = seq(4)) |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      proba = list(model = "uniform", status = 1),
      proba = list(value = .layer, status = c(2, 3, 4)),
      .layer = list(model = "uniform")
    )

  withr::with_seed(1, ll1 <- hhmodel$evaluate_log_lik())
  # check if result is correct
  expect_equal(ll1, sum(log(hhmodel$export_augmented_data()$proba)))
  # verify it's different if different seed
  withr::with_seed(2, ll2 <- hhmodel$evaluate_log_lik())
  expect_true(ll1 != ll2)
  # opposite if same seed
  withr::with_seed(2, ll2encore <- hhmodel$evaluate_log_lik())
  expect_identical(ll2, ll2encore)
})

test_that("Status-specific augmented variables draw from per-status ranges", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      proba = list(status = 1, model = "uniform", min = 0, max = 1),
      proba = list(status = 2, model = "uniform", min = 1, max = 2),
      proba = list(status = 3, model = "uniform", min = 2, max = 3),
      proba = list(status = 4, model = "uniform", min = 3, max = 4)
    )

  # simulation also routes draws into per-status supports
  withr::with_seed(1, sim <- hhmodel$run_simulation(nrep = 50))
  sim <- do.call(rbind, sim)
  for (s in seq(4)) {
    vals <- sim$proba[sim$status == s]
    expect_true(all(vals >= s - 1 & vals <= s))
  }

  # run_inference with adapt > 0 but no parametric proposals: the adapt budget
  # is silently repurposed as burn-in; only the short-run diagnostic fires
  # (steps = 5 < 400).
  expect_warning(
    results <- hhmodel |>
      run_inference(steps = 5, thin = 1, adapt = 2,
                    chain = 1, threads = 1, progress = FALSE,
                    id_draws = hhmodel$data_wide$id),
    regexp = "below 400"
  )
  draws <- as.data.frame(get_augmented_draws(results))
  for (s in seq(4)) {
    col <- draws[[as.character(s)]]
    expect_true(all(col >= s - 1 & col <= s))
  }
})

test_that("Augmented aliasing works as expected", {
  hhmodel <- data.frame(house_id = "hh1",
                        id = seq(1, 3),
                        x = 2 + runif(3),
                        status = 1) |>
    new_housepi_model("augmented_example")
  hhmodel$add_parameter(name = "scale", value = 1)
  hhmodel$add_augmented_list(
    .dummy = x@D,
    .layer = list(func = "gamma",
                  func_args = c(x = .dummy, shape = 2, scale = scale@P)),
    proba = .layer
  )

  ll <- sum(log(dgamma(hhmodel$data_wide$x, shape = 2, scale = 1)))
  expect_equal(hhmodel$evaluate_log_lik(), ll)

  # check that order does not matter
  hhmodel$add_augmented_list(
    .layer = list(func = "gamma",
                  func_args = c(x = .dummy, shape = 2, scale = scale@P)),
    proba = .layer,
    .dummy = x@D
  )
  expect_equal(hhmodel$evaluate_log_lik(), ll)

  # check normalization of binomial with augmented example
  size <- 10
  hhmodel <- data.frame(id = seq(0, size),
                        house_id = "hh1", status = 1,
                        X = seq(0, size)) |>
    new_housepi_model("augmented_example")
  hhmodel$add_augmented_list(
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

  build_model <- function(data) {
    hh <- new_housepi_model(data, "augmented_example", quiet = TRUE)
    hh$add_augmented_list(
      proba = list(func = "gamma",
                   func_args = c(x = ip@D, shape = !!shape, scale = !!scale))
    )
    hh
  }

  # warning for no structure
  expect_warning(new_housepi_model(df, "augmented_example"),
                 regexp = "Ignored household structure")
  hhmodel <- build_model(df)
  expect_equal(hhmodel$evaluate_log_lik(), exact_val(df$ip, shape, scale))

  # repeat with degenerated ip values to trigger multiplicity
  df$ip <- floor(df$ip)
  hhmodel <- build_model(df)
  expect_equal(hhmodel$evaluate_log_lik(), exact_val(df$ip, shape, scale))
})

test_that("Known errors are unchanged", {
  # Known infinity (values below minimum, i. e. proba 0)
  hhmodel <- new_housepi_model(sarscov2_2021_icofs, "augmented_example")
  hhmodel |>
    add_augmented_list(
      value = list(model = "uniform", min = 0, max = 1),
      proba = list(func = "lognormal",
                   func_args = c(x = value, meanlog = 1, sdlog = 1),
                   min = 1, max = 15, discr = TRUE)
    ) |>
    evaluate_log_lik() |>
    expect_equal(-Inf)

  # Known NaN: qlnorm(15, meanlog = 4, sdlog = 0.01) = 0, impossible truncation
  # UPDATE: shouldn't this raise a warning instead?
  hhmodel |>
    add_augmented_list(
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

  # Looping over thread counts exercises the RcppParallel reduction in
  # commit_total_log_lik_; both invariants must hold for any thread count.
  for (t in c(1L, 2L)) {
    results <- suppressWarnings(
      new_housepi_model(house, "augmented_example") |>
        add_augmented_list(proba = list(model = "gamma", x = x@D,
                                        shape = !!shape, scale = !!scale)) |>
        run_inference(steps = 5, thin = 1, adapt = 0, chain = 1, threads = t,
                      progress = FALSE, id_draws = id_draws)
    )

    # coherency between pointwise and total log-likelihood
    expect_equal(get_total_log_lik_draws(results)$log_lik,
                 get_log_lik_draws(results)$hh01)

    # coherency between augmented draws and log-likelihood
    df <- as.data.frame(get_augmented_draws(results))
    df <- as.matrix(df[, !startsWith(names(df), "."), drop = FALSE])
    expect_equal(get_log_lik_draws(results)$hh01,
                 rowSums(log(df * dgamma(df, shape = shape, scale = scale))))
  }
})

test_that("export_augmented_data works after evaluate_log_lik", {
  hhmodel <- new_housepi_model(
    data.frame(house_id = 1, id = c("a", "b"), status = c(1, 0)),
    "augmented_example"
  )
  hhmodel$add_augmented(name = "proba", model = "uniform")
  hhmodel$evaluate_log_lik()
  result <- hhmodel$export_augmented_data()
  expect_true(is.data.frame(result))
  expect_true("proba" %in% colnames(result))
})

test_that("export_house_snapshot returns correct structure", {
  hhmodel <- new_housepi_model(flu_1918_frederick, "chain_binomial")
  hhmodel$add_parameter_template()
  hhmodel$evaluate_log_lik()

  hid <- hhmodel$data_wide$house_id[1]
  snapshot <- hhmodel$export_house_snapshot(hid)
  expect_true(is.list(snapshot))
  expect_true("data_wide" %in% names(snapshot))
  expect_true("param" %in% names(snapshot))
  expect_true(is.data.frame(snapshot$data_wide))
})
