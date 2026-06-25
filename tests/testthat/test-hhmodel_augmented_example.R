library(testthat)

test_that("Nested definitions raise errors", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("augmented_example")

  expect_no_error(
    hhmodel |>
      add_augmented_list(
        .layer = list(prior = "uniform", status = 1),
        .layer = list(value = .dummy, status = c(2, 3, 4)),
        .dummy = list(prior = "uniform"),
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
      .layer = list(prior = "normal", status = 1),
      .layer = list(value = .dummy, status = 2),
      #> NO .layer provided for status 3 and 4
      .dummy = list(prior = "lognormal", status = 2),
      proba = .layer
    )$evaluate_log_lik()
  ))
  expect_match(err_msg, "Could not find a definition")
  expect_match(err_msg, "status")

  # second generation of definitions
  err_msg <- as.character(expect_error(
    hhmodel$add_augmented_list(
      .layer = list(prior = "normal", status = 1),
      .layer = list(value = .dummy, status = c(2, 3, 4)),
      #> NO .dummy provided for status 3, 4
      .dummy = list(prior = "lognormal", status = c(2, 3)),
      proba = .layer,
      append = FALSE
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
  hhmodel |>
    add_augmented(name = "proba", func = "gamma") |>
    expect_error(regexp = "Missing non default")

  # data with NA values: the NA flows through to the log-probability NaN guard
  hhmodel |>
    drop_augmented() |>
    add_augmented(name = "proba", func = "gamma",
                  func_args = c(shape = 2)) |>
    evaluate_log_lik() |>
    expect_error(regexp = "returned NaN")

  # concatenated data with NA values
  hhmodel |>
    add_augmented_list(
      gamma = list(func = "gamma",  func_args = c(shape = 2)),
      proba =  exp(-gamma)
    ) |>
    evaluate_log_lik() |>
    expect_error(regexp = "Impossible to associate")

  # check that NaN + something causes errors
  hhmodel |>
    add_augmented_list(noise = list(prior = "normal"),
                       proba = exp(test@D + noise)) |>
    evaluate_log_lik() |>
    expect_error(regexp = "Impossible to associate")

  # check that log(NaN) causes error
  hhmodel |>
    add_augmented_list(proba = log(test@D)) |>
    evaluate_log_lik() |>
    expect_error(regexp = "Impossible to associate")

  # check that log(0) causes error
  hhmodel |>
    add_augmented_list(noise = 0., proba = log(noise)) |>
    evaluate_log_lik() |>
    expect_error(regexp = "Impossible to associate")
})

test_that("Augmented log-likelihood seeding is fine", {
  hhmodel <- data.frame(house_id = 1, id = seq(4),
                        status = seq(4)) |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      proba = list(prior = "uniform", status = 1),
      proba = list(value = .layer, status = c(2, 3, 4)),
      .layer = list(prior = "uniform")
    )

  withr::with_seed(1, ll1 <- eval_log_density(hhmodel))
  # check if result is correct
  expect_equal(ll1, sum(log(hhmodel$export_augmented_data()$proba)))
  # verify it's different if different seed
  withr::with_seed(2, ll2 <- eval_log_density(hhmodel))
  expect_true(ll1 != ll2)
  # opposite if same seed
  withr::with_seed(2, ll2encore <- eval_log_density(hhmodel))
  expect_identical(ll2, ll2encore)
})

test_that("Status-specific augmented variables draw from per-status ranges", {
  hhmodel <- data.frame(house_id = 1, id = seq(4), status = seq(4)) |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      proba = list(status = 1, prior = "uniform", min = 0, max = 1),
      proba = list(status = 2, prior = "uniform", min = 1, max = 2),
      proba = list(status = 3, prior = "uniform", min = 2, max = 3),
      proba = list(status = 4, prior = "uniform", min = 3, max = 4)
    )

  # simulation also routes draws into per-status supports
  withr::with_seed(1, sim <- hhmodel$run_simulation(nrep = 50))
  sim <- do.call(rbind, sim)
  for (s in seq(4)) {
    vals <- sim$proba[sim$status == s]
    expect_true(all(vals >= s - 1 & vals <= s))
  }

  # case adapt > 0 but no parametric proposals
  # adapt budget is silently repurposed as burn-in
  # short-run diagnostic fires (steps = 5 < 400)
  expect_warning(
    results <- hhmodel |>
      run_inference(steps = 5, thin = 1, adapt = 2,
                    chain = 1, threads = 1,
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
  hhmodel <- data.frame(house_id = "hh1", id = seq(1, 3),
                        x = 2 + runif(3), status = 1) |>
    new_housepi_model("augmented_example") |>
    add_parameter(name = "scale", value = 1) |>
    add_augmented_list(
      .dummy = x@D,
      .layer = list(func = "gamma",
                    func_args = c(x = .dummy, shape = 2, scale = scale@P)),
      proba = .layer
    )

  ll <- sum(log(dgamma(hhmodel$data_wide$x, shape = 2, scale = 1)))
  expect_equal(eval_log_density(hhmodel), ll)

  # check that order does not matter
  hhmodel |>
    drop_augmented() |> 
    add_augmented_list(
      .layer = list(func = "gamma",
                    func_args = c(x = .dummy, shape = 2, scale = scale@P)),
      proba = .layer,
      .dummy = x@D
    ) |>
    eval_log_density() |>
    expect_equal(ll)

  # check normalization of binomial with augmented example
  size <- 10
  data.frame(id = seq(0, size), house_id = "hh1",
             status = 1, X = seq(0, size)) |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      pdf = list(func = "binomial",
                 func_args = c(x = X@D, size = !!size, prob = 0.5)),
      proba = exp(pdf)
    ) |>
    eval_log_density("hh1") |>
    expect_equal(1.)
})

test_that("Augmentation example on incubation from delays is coherent", {
  shape <- 1.95
  scale <- 1.13
  get_exact_val <- function(x, shape, scale) {
    sum(log(stats::dgamma(x, shape = shape, scale = scale)))
  }

  # log-likelihood from simulated delays

  df <- data.frame(status = 1,
                   ip = rgamma(223, shape = shape, scale = scale))

  build_model <- function(data) {
    data |>
      new_housepi_model("augmented_example", quiet = TRUE) |>
      add_augmented(
        "proba", func = "gamma",
        func_args = c(x = ip@D, shape = !!shape, scale = !!scale)
      )
  }

  # warning for no structure
  new_housepi_model(df, "augmented_example") |>
    expect_warning(regexp = "Ignored household structure")
  hhmodel <- build_model(df)
  expect_equal(eval_log_density(hhmodel), get_exact_val(df$ip, shape, scale))

  # repeat with degenerated ip values to trigger multiplicity
  df$ip <- floor(df$ip)
  hhmodel <- build_model(df)
  expect_equal(eval_log_density(hhmodel), get_exact_val(df$ip, shape, scale))
})

test_that("Known errors are unchanged", {
  # known infinity: values below minimum, i. e. probability is 0, log -Inf
  sarscov2_2021_icofs |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      value = list(prior = "uniform", min = 0, max = 1),
      proba = list(func = "lognormal",
                   func_args = c(x = value, meanlog = 1, sdlog = 1),
                   min = 1, max = 15, discr = TRUE)
    ) |>
    eval_log_density() |>
    expect_equal(-Inf)

  # Known NaN: qlnorm(15, meanlog = 4, sdlog = 0.01) = 0, impossible truncation
  # UPDATE: shouldn't this raise a warning instead?
  sarscov2_2021_icofs |>
    new_housepi_model("augmented_example") |>
    add_augmented_list(
      value = list(prior = "uniform", min = 1, max = 10, discr = TRUE),
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
      house |>
        new_housepi_model("augmented_example") |>
        add_augmented("proba", prior = "gamma", x = x@D,
                      shape = !!shape, scale = !!scale) |>
        run_inference(steps = 5, thin = 1, adapt = 0, chain = 1,
                      threads = t, id_draws = id_draws)
    )

    # coherency between augmented draws and log-density: the pointwise
    # likelihood plus the augmented prior recompose log(proba) +
    # log(dgamma(proba)). With no parameter priors, lprior == augm_lprior.
    df <- as.data.frame(get_augmented_draws(results))
    df <- as.matrix(df[, !startsWith(names(df), "."), drop = FALSE])
    augm_prior <- get_lprior_draws(results)$lprior
    expect_equal(get_log_lik_draws(results)$hh01 + augm_prior,
                 rowSums(log(df * dgamma(df, shape = shape, scale = scale))))
  }
})

test_that("export_augmented_data works after evaluate_log_lik", {
  hhmodel <- data.frame(house_id = 1, id = c("a", "b"), status = c(1, 0)) |>
    new_housepi_model("augmented_example") |>
    add_augmented(name = "proba", prior = "uniform")
  # UPDATE: I can test here that before export_augmented_data fails...
  hhmodel$evaluate_log_lik()
  result <- hhmodel$export_augmented_data()
  expect_true(is.data.frame(result))
  expect_true("proba" %in% colnames(result))
})

test_that("export_house_snapshot returns correct structure", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  # UPDATE: I can test here that before export_house_snapshot fails...
  hhmodel$evaluate_log_lik()

  hid <- hhmodel$data_wide$house_id[1]
  snapshot <- hhmodel$export_house_snapshot(hid)
  expect_true(is.list(snapshot))
  expect_true("data_wide" %in% names(snapshot))
  expect_true("param" %in% names(snapshot))
  expect_true(is.data.frame(snapshot$data_wide))
})
