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

test_that("Deep clone works as expected", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$data_wide <- housepi::survey
  hhmodel$load_parameter_template()
  expect_no_error(out <- hhmodel$infer(steps = 2, progress = FALSE))

  # shallow clone
  hhclone <- hhmodel$clone(deep = FALSE)

  # check backward changes
  hhclone$add_group_rate("infectivity", by = "sex")
  expect_equal(hhmodel$group_rates, hhclone$group_rates)

  # deep clone
  hhclone <- hhmodel$clone(deep = TRUE)

  # data
  expect_equal(hhmodel$data_wide, hhclone$data_wide)
  expect_equal(hhmodel$data_long, hhclone$data_long)

  # definitions
  expect_equal(hhmodel$setup, hhclone$setup)
  expect_equal(hhmodel$parameters, hhclone$parameters)
  expect_equal(hhmodel$augmented, hhclone$augmented)
  expect_equal(hhmodel$group_rates, hhclone$group_rates)

  # check no backward changes
  hhmodel$add_group_rate("susceptibility", by = "sex")
  expect_failure(expect_equal(hhmodel$group_rates, hhclone$group_rates))
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

test_that("Loglikelihood traces and pointwise traces are coherent", {
  results <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    add_data_wide(housepi::survey) |> 
    infer(steps = 10, thin = 1, adapt = 0, burnin = 0)

  pointwise <- results$house_draws |>
    posterior::as_draws_matrix() |>
    rowSums() |>
    unname()
  expect_equal(results$loglike_draws$loglikelihood, pointwise)
})

test_that("Chain binomial models are coherent", {
  # loglikelihood from chain binomial model
  hhmodel1 <- new_housepi_model("chain_binomial")

  hhmodel1$parameters <- list(alpha = 0.15, beta = 0.22)
  hhmodel1$data_wide <- housepi::survey
  ll1 <- hhmodel1$loglikelihood()

  # loglikelihood from heterogeneous chain binomial model
  # with homogenous rates
  hhmodel2 <- new_housepi_model("chain_binomial_heterogeneous")

  hhmodel2$parameters <- list(alpha = 0.15, beta = 0.22)
  hhmodel2$data_wide <- housepi::survey
  ll2 <- hhmodel2$loglikelihood()

  expect_equal(ll1, ll2)
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

test_that("Loglikelihood of augmented infection time model is correct", {
  # instantaneous infectivity profile
  infec_profile_inst <- function(x, param) {
    dgamma(x, shape = param$shape, scale = param$scale)
  }
  # cumulative infectivity profile
  infec_profile_cumul <- function(x, param) {
    pgamma(x, shape = param$shape, scale = param$scale)
  }
  # transmission rate
  transmission_rate <- function(house, param) {
    param$beta * (param$size_ref / dim(house)[1]) ^ param$gamma
  }
  # loglikelihood contribution of susceptible members
  loglike_susceptible <- function(house, param) {
    i0 <- which(house$infec_time == min(house$infec_time, na.rm = TRUE))
    member_0 <- house[i0, ]

    alpha_cumul <- param$alpha * (param$study_period + member_0$incub_period)
    beta <- transmission_rate(house, param)

    output <- 0
    for (i in which(house$status == 0)) {
      roi_cumul <- 0
      for (k in which(house$status == 1)) {
        t_infec_k <- house[k, "infec_time"]
        exposure <- param$study_period + member_0$onset_time - t_infec_k
        roi_cumul <- roi_cumul + infec_profile_cumul(exposure, param)
      }
      output <- output - alpha_cumul - beta * roi_cumul
    }

    output
  }
  # loglikelihood contribution of infectious members
  loglike_infectious <- function(house, param) {
    output <- 0.
    i0 <- which(house$infec_time == min(house$infec_time, na.rm = TRUE))
    member_0 <- house[i0, ]
    t_intro <- member_0$infec_time

    alpha <- param$alpha
    beta <- transmission_rate(house, param)

    for (i in which(house$status == 1)) {
      if (i != i0) {
        t_infec_i <- house[i, "infec_time"]
        roi_inst_i <- 0.
        roi_cumul_i <- 0.
        for (k in which(house$status == 1)) {
          t_infec_k <- house[k, "infec_time"]
          gentime <- t_infec_i - t_infec_k
          if (gentime > 0) {
            roi_inst_i <- roi_inst_i + infec_profile_inst(gentime, param)
            roi_cumul_i <- roi_cumul_i + infec_profile_cumul(gentime, param)
          }
        }

        roi_inst_i <- roi_inst_i * beta + alpha
        roi_cumul_i <- roi_cumul_i * beta + alpha * (t_infec_i - t_intro)
        output <- output + log(roi_inst_i) - roi_cumul_i
      }
    }

    output
  }
  # dummy dataset
  house <- data.frame(id = sprintf("hh_01_%s", seq(1, 7)),
                      house_id = "hh_01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0.0, 2.7, 3.6, NA, NA, NA, NA),
                      incub_period = c(2.5, 3.1, 2.0, NA, NA, NA, NA))
  house$onset_time <- house$infec_time + house$incub_period
  # chioce for the parameters
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                size_ref = 4, study_period = 45, shape = 2, scale = 1.1)

  hhmodel <- new_housepi_model("augmented_infection_time")
  hhmodel$data_wide <- house
  hhmodel$parameters <- param
  hhmodel$augmented <- list(infec_time = ~infec_time@D,
                            incub_period = ~incub_period@D,
                            onset_time = ~onset_time@D)
  hhmodel$add_infec_profile(name = "gamma",
                            args = c(shape = ~shape@P, scale = ~scale@P))

  expect_equal(
    hhmodel$loglikelihood(),
    loglike_susceptible(house, param) + loglike_infectious(house, param)
  )
})

test_that("Simulations behave as expected", {
  hhsimul <- new_housepi_model("augmented_infection_time")
  hhsimul$data_wide <- dplyr::mutate(housepi::survey, status = NA)
  hhsimul$load_parameter_list(alpha = 0,
                              beta = 0,
                              study_period = 21,
                              gamma = 1,
                              size_ref = 4)
  # UPDATE: these dependencies should not be in template

  hhsimul$load_augmented_list(
    incub_period = list(model = "gamma", args = c(shape = 2, scale = 1)),
    infec_time = ~- incub_period,
    onset_time = ~infec_time + incub_period
  )

  # simulate just introduction cases
  simulated <- hhsimul$simulate()
  expect_identical(
    length(unique(simulated$house_id)),
    sum(simulated$onset_time == 0, na.rm = TRUE)
  )

  # simulate all infected
  expect_warning(hhsimul$add_parameter(name = "beta", value = 100))
  simulated <- hhsimul$simulate()
  expect_true(all(simulated$status > 0))
  expect_false(any(is.na(simulated$infec_time)))
  expect_equal(
    simulated$infec_time,
    simulated$onset_time - simulated$incub_period
  )
})
