library(testthat)

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

# household log-likelihood of augmented infection time model
augm_infec_obs_log_lik <- function(house, param, discr = FALSE) {
  loglike_susceptible(house, param, discr) +
    loglike_infectious(house, param, discr)
}

test_that("Simulations behave as expected", {
  hhsimul <- new_housepi_model("augmented_infection_time")
  hhsimul$data_wide <- dplyr::mutate(housepi::survey, status = NA)
  hhsimul$load_parameter_list(alpha = 0,
                              beta = 0,
                              study_period = 21,
                              gamma = 1,
                              size_ref = 4)
  # UPDATE: these dependencies should not be in template...

  hhsimul$load_augmented_list(
    incub_period = list(model = "gamma", shape = 2, scale = 1),
    infec_time = -incub_period,
    onset_time = infec_time + incub_period
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

test_that("Log-likelihood of augmented infection time model is correct", {
  # dummy dataset
  house <- data.frame(id = sprintf("hh01m%s", seq(1, 7)),
                      house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0, 3, 4, NA, NA, NA, NA),
                      onset_time = c(2, 6, 6, NA, NA, NA, NA))
  house$incub_period <- house$onset_time - house$infec_time

  # choice for the parameters
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                size_ref = 4, study_period = 45, shape = 2, scale = 1.1)
  # choice for the augmented variables
  augm <- list(infec_time = ~ infec_time@D,
               incub_period = ~ incub_period@D,
               onset_time = ~ onset_time@D)

  for (discr in c(FALSE, TRUE)) {
    hhmodel <- new_housepi_model("augmented_infection_time",
                                 discr = discr, quiet = TRUE) |>
      add_data_wide(house) |>
      load_parameter_list(!!!param) |>
      load_augmented_list(!!!augm) |>
      add_infec_profile(func = "gamma", discr = discr,
                        shape = shape@P, scale = scale@P)

    expect_equal(evaluate_log_lik(hhmodel),
                 augm_infec_obs_log_lik(house, param, discr))
  }
})

test_that("Augmented moves preserve log-likelihood", {
  # dummy dataset
  house <- data.frame(id = sprintf("hh01m%s", seq(1, 7)),
                      house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA))

  hhmodel <- new_housepi_model("augmented_infection_time") |>
    add_data_wide(house) |>
    load_parameter_template() |>
    load_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
    load_augmented_list(
      infec_time = onset_time@A - incub_period@A,
      onset_time = list(value = onset_date@D, 
                        model = "normal", mean = onset_date@D),
      incub_period = list(model = "weibull", shape = 2, scale = 1)
    ) |>
    add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)

  n_steps <- 5
  id_draws <- house[house$status > 0,]$id
  results <- hhmodel |>
    run_inference(steps = n_steps, thin = 1, adapt = 0, chain = 1,
                  progress = FALSE, id_draws = id_draws)

  tot_log_lik <- get_total_log_lik_draws(results)$log_lik
  expect_equal(tot_log_lik, get_log_lik_draws(results)$hh01)

  param_draws <- get_param_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()

  augm_draws <- get_augm_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()

  # constant parameters
  consts <- hhmodel$parameters[hhmodel$parameters$const, c("name", "value")] |>
    tidyr::pivot_wider(names_from = "name", values_from = "value") |>
    as.list()

  # state of augmented values per iteration
  get_augm_state <- function(df, index) {
    df[index, ] |>
      tidyr::pivot_longer(dplyr::everything(),
                          names_to = "name", values_to = "value") |>
      tidyr::separate(name, into = c("id", "augm"),
                      sep = "_", extra = "merge") |>
      tidyr::pivot_wider(names_from = augm, values_from = value) |>
      as.data.frame()
  }
  
  join_augm_data <- function(i) {
    dplyr::left_join(house, get_augm_state(augm_draws, i), by = "id")
  }

  # get total log-likelihood (observation + augmented)
  get_log_lik_i <- function(i) {
    df <- join_augm_data(i)
    out <- augm_infec_obs_log_lik(df, c(param_draws[i, ], consts))
    df <- df[df$status > 0,]
    out + 
      sum(dnorm(df$onset_time - df$onset_date, log = TRUE)) +
      sum(dweibull(df$incub_period, shape = 2, scale = 1, log = TRUE))
  }
  
  # check that snapshot is last state
  snapshot <- hhmodel$export_house_snapshot("hh01")
  expect_identical(unlist(c(param_draws[n_steps, ], consts)), 
                   snapshot$param)
  expect_identical(get_augm_state(augm_draws, n_steps), 
                   snapshot$augm[house$status > 0,])
  
  # check if log-likelihood is coherent
  hhmodel |>
    add_data_wide(join_augm_data(n_steps)) |>
    load_parameter_list(!!!snapshot$param) |>
    load_augmented_list(
      infec_time = infec_time@D,
      onset_time = list(value = onset_time@D, const = TRUE,
                        model = "normal", mean = onset_date@D),
      incub_period = list(value = incub_period@D, const = TRUE,
                          model = "weibull", shape = 2, scale = 1)
    ) |>
    evaluate_log_lik() |> 
    expect_equal(tot_log_lik[n_steps])

  # check that total log-likelihood is well computed
  expect_equal(tot_log_lik, sapply(seq_along(tot_log_lik), get_log_lik_i))
})
