library(testthat)

test_that("Log-likelihood of augmented infection time model is correct", {
  # cumulative infectivity profile
  infec_profile_cumul <- function(x, param, discr) {
    if (!discr) return(pgamma(x, shape = param$shape, scale = param$scale))
    z <- pgamma(1, shape = param$shape, scale = param$scale)
    (pgamma(x + 1, shape = param$shape, scale = param$scale) - z) / (1 - z)
  }

  # instantaneous infectivity profile
  infec_profile_inst <- function(x, param, discr) {
    if (!discr) return(dgamma(x, shape = param$shape, scale = param$scale))
    z <- pgamma(1, shape = param$shape, scale = param$scale)
    (pgamma(x + 1, shape = param$shape, scale = param$scale) -
       pgamma(x, shape = param$shape, scale = param$scale)) / (1 - z)
  }

  # extract individual introduction
  find_introduction <- function(house) {
    subset(subset(house, status > 0), infec_time == min(infec_time))
  }

  # transmission rate
  transmission_rate <- function(house, param) {
    param$beta * (param$size_ref / dim(house)[1]) ^ param$gamma
  }

  # log-likelihood contribution of susceptible members
  loglike_susceptible <- function(house, param, discr) {
    member_0 <- find_introduction(house)

    alpha_cumul <- param$alpha * (param$study_period + member_0$incub_period)
    beta <- transmission_rate(house, param)

    output <- 0
    for (i in which(house$status == 0)) {
      roi_cumul <- 0
      for (k in which(house$status == 1)) {
        t_infec_k <- house[k, "infec_time"]
        exposure <- param$study_period + member_0$onset_time - t_infec_k
        roi_cumul <- roi_cumul + infec_profile_cumul(exposure, param, discr)
      }
      output <- output - alpha_cumul - beta * roi_cumul
    }
    output
  }

  # log-likelihood contribution of infectious members
  loglike_infectious <- function(house, param, discr) {
    member_0 <- find_introduction(house)
    i0 <- rownames(member_0)
    t_intro <- member_0$infec_time

    alpha <- param$alpha
    beta <- transmission_rate(house, param)

    output <- 0.
    for (i in which(house$status == 1)) {
      if (i != i0) {
        t_infec_i <- house[i, "infec_time"]
        roi_inst_i <- 0.
        roi_cumul_i <- 0.
        for (k in which(house$status == 1)) {
          t_infec_k <- house[k, "infec_time"]
          gentime <- t_infec_i - t_infec_k
          if (gentime > 0) {
            roi_inst_i <- roi_inst_i +
              infec_profile_inst(gentime, param, discr)
            roi_cumul_i <- roi_cumul_i +
              infec_profile_cumul(gentime - discr, param, discr)
          }
        }
        roi_cumul_i <- roi_cumul_i * beta + alpha * (t_infec_i - t_intro)
        roi_inst_i <- roi_inst_i * beta + alpha
        if (discr) {
          output <- output + log(1. - exp(-roi_inst_i)) - roi_cumul_i + alpha
        } else {
          output <- output + log(roi_inst_i) - roi_cumul_i
        }
      }
    }

    output
  }

  house_obs_loglike <- function(house, param, discr) {
    loglike_susceptible(house, param, discr) +
      loglike_infectious(house, param, discr)
  }

  # dummy dataset
  house <- data.frame(id = sprintf("hh_01_%s", seq(1, 7)),
                      house_id = "hh_01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0, 3, 4, NA, NA, NA, NA),
                      incub_period = c(2, 3, 2, NA, NA, NA, NA))
  house$onset_time <- house$infec_time + house$incub_period

  # choice for the parameters
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                size_ref = 4, study_period = 45, shape = 2, scale = 1.1)

  for (discr in c(FALSE, TRUE)) {
    hhmodel <- new_housepi_model("augmented_infection_time", 
                                 discr = discr, quiet = TRUE) |>
      add_data_wide(house) |>
      load_parameter_list(from = param) |>
      load_augmented_list(infec_time = ~infec_time@D,
                          incub_period = ~incub_period@D,
                          onset_time = ~onset_time@D) |>
      add_infec_profile(name = "gamma", discr = discr,
                        shape = ~shape@P, scale = ~scale@P)

    expect_equal(evaluate_log_lik(hhmodel),
                 house_obs_loglike(house, param, discr))
  }
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
    incub_period = list(model = "gamma", shape = 2, scale = 1),
    infec_time = ~ -incub_period,
    onset_time = ~ infec_time + incub_period
  )

  # simulate just introduction cases
  simulated <- hhsimul$run_simulation()
  expect_identical(
    length(unique(simulated$house_id)),
    sum(simulated$onset_time == 0, na.rm = TRUE)
  )

  # simulate all infected
  expect_warning(hhsimul$add_parameter(name = "beta", value = 100))
  simulated <- hhsimul$run_simulation()
  expect_true(all(simulated$status > 0))
  expect_false(any(is.na(simulated$infec_time)))
  expect_equal(
    simulated$infec_time,
    simulated$onset_time - simulated$incub_period
  )
})
