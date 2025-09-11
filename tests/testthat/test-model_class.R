library(testthat)

test_that("Model class inconsistent calls", {
  # Default blank
  expect_warning(new_housepi_model(), regexp = "Falling back to")

  # blank household model object
  hhmodel <- new_housepi_model("blank")

  # inference without definitions
  expect_error(hhmodel$infer(),
               regexp = "Nothing to sample")
})

test_that("Model class blank calls", {
  # blank household model object
  hhmodel <- new_housepi_model("blank")

  # blank default data wide format
  data_wide <- housepi:::template_blank_data_wide_()
  expect_equal(hhmodel$data_wide, data_wide)

  # correct data long format
  data_long <- data.frame(id = data_wide$id[1],
                          time = seq(3),
                          measure = seq(3)**3)
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
  expect_error(hhmodel$simulate())
})

test_that("Deep clone works as expected", {
  # UPDATE: change to a model with data augmentation
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$load_parameter_template()
  out <- hhmodel$infer(housepi::survey, progress = FALSE)

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

test_that("Loglikelihood traces and pointwise traces are coherent", {
  results <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    infer(housepi::survey, sample_size = 10, thin = 5,
          n_adapt = 50, n_burnin = 50)

  pointwise <- results$house_draws |>
    posterior::as_draws_matrix() |>
    rowSums() |>
    unname()
  expect_equal(results$loglike_draws$loglikelihood, pointwise)
})

test_that("Chain binomial models are coherent", {
  # loglikelihood from chain binomial model
  hhmodel1 <- new_housepi_model("chain_binomial")

  hhmodel1$load_parameter_list(alpha = 0.15, beta = 0.22)
  ll1 <- hhmodel1$loglikelihood(housepi::survey)

  # loglikelihood from heterogeneous chain binomial model
  # with homogenous rates
  hhmodel2 <- new_housepi_model("chain_binomial_heterogeneous")

  hhmodel2$load_parameter_list(alpha = 0.15, beta = 0.22)
  ll2 <- hhmodel2$loglikelihood(housepi::survey)

  expect_equal(ll1, ll2)
})

test_that("Nested definitions are raise errors", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$data_wide <- data.frame(house_id = 1, id = seq(3), status = 1)

  # simple nest
  expect_error(
    hhmodel |>
      load_augmented_list(dummy = ~log_proba,
                          log_proba = ~dummy) |>
      loglikelihood(),
    regexp = "nested definitions"
  )

  # triangular nest
  expect_error(
    hhmodel |>
      load_augmented_list(log_proba = ~layer,
                          layer = ~dummy,
                          dummy = ~log_proba) |>
      loglikelihood(),
    regexp = "nested definitions"
  )
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
                 func_args = c(x = ~ip@D, shape = shape, scale = scale)),
    log_proba = ~log(proba)
  )

  # warning for no structure
  expect_warning(x <- hhmodel$loglikelihood(df))
  expect_equal(x, exact)
})

test_that("Arithmetic function with undefined x return NaN", {
  hhmodel <- new_housepi_model("augmented_example")

  # NA values
  hhmodel$data_wide <- data.frame(house_id = seq(3), id = seq(3),
                                  status = 1, test = c(5.5, NA, -1.5))
  # check that NaN + something = NaN
  hhmodel$load_augmented_list(noise = list(model = "norm"),
                              log_proba = ~test@D + noise)
  # check that NA loglikelihood
  expect_error(hhmodel$infer(), regexp = "returned NaN")

  # check that log(NaN) = NaN
  hhmodel$load_augmented_list(log_proba = ~log(test@D))
  expect_error(hhmodel$infer(), regexp = "Nothing to sample")
  expect_error(hhmodel$loglikelihood(), regexp = "returned NaN")
})

test_that("Augmented function example returns expected distribution", {
  hhmodel <- new_housepi_model("blank")

  # define Box-Muller transformation
  # UPDATE: there should be a developer effort to improve this
  hhmodel$load_augmented_list(
    norm = list(model = "norm"),
    u1 = list(model = "unif"),
    u2 = list(model = "unif"),
    log_u1 = ~log(u1),
    two_log_u1 = ~2 * log_u1,
    m_two_log_u1 = ~ - two_log_u1,
    sqrt_m_two_log_u1 = ~sqrt(m_two_log_u1),
    two_pi_u2 = ~6.283185 * u2,
    cos_two_pi_u2 = ~cos(two_pi_u2),
    z0 = ~sqrt_m_two_log_u1 * cos_two_pi_u2
  )

  results <- hhmodel$infer(n_steps = 5e5, progress = 0)
  sample_norm <- results$augm_draws$hh_01_1_norm
  sample_boxmul <- results$augm_draws$hh_01_1_z0

  expect_gt(stats::ks.test(sample_norm, sample_boxmul)$p.value, 0.05)
  expect_gt(stats::t.test(sample_norm, sample_boxmul)$p.value, 0.05)
  expect_gt(stats::var.test(sample_norm, sample_boxmul)$p.value, 0.05)
})

test_that("Augmented infection times model has correct loglikelihood", {
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

    return(output)
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

    return(output)
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
