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
  t_intro <- member_0$infec_time
  alpha <- param$alpha
  beta <- transmission_rate(house, param)

  output <- 0
  for (k in which(house$status == 0)) {
    t_end_k <- house[k, "end_followup"]
    roi_cumul <- 0
    alpha_cumul <- alpha * (t_end_k - t_intro)
    for (i in which(house$status == 1)) {
      t_infec_i <- house[i, "infec_time"]
      exposure <- t_end_k - t_infec_i
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
        output <- output + log1p(- exp(-roi_inst_i)) - roi_cumul_i + alpha
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
  hhsimul <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 200)) |>
    dplyr::mutate(status = NA) |>
    new_housepi_model("augmented_infection_time")
  hhsimul$add_parameter_list(alpha = 0, beta = 0,
                              gamma = 1, size_ref = 4)
  hhsimul$add_augmented_list(
    incub_period = list(model = "gamma", shape = 2, scale = 1),
    infec_time = -incub_period
  )

  obs_process <- function(df) {
    df$onset_time <- df$infec_time + df$incub_period
    df
  }
  # simulate just introduction cases
  simulated <- hhsimul$run_simulation(obs_process = obs_process)
  expect_identical(
    length(unique(simulated$house_id)),
    sum(simulated$onset_time == 0, na.rm = TRUE)
  )

  # simulate all infected
  expect_warning(hhsimul$add_parameter(name = "beta", value = 100),
                 regexp = "Overwriting")
  simulated <- hhsimul$run_simulation(obs_process = obs_process)
  expect_true(all(simulated$status > 0))
  expect_false(any(is.na(simulated$infec_time)))
  expect_equal(
    simulated$infec_time,
    simulated$onset_time - simulated$incub_period
  )
})

test_that("Log-likelihood of augmented infection time model is correct", {
  # dummy dataset
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0, 3, 4, NA, NA, NA, NA),
                      end_followup = 21)
  house$id <- sprintf("hh01m%s", seq(1, dim(house)[1]))

  # choice for the parameters
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                size_ref = 4, shape = 2, scale = 1.1)

  for (discr in c(FALSE, TRUE)) {
    hhmodel <- house |>
      new_housepi_model("augmented_infection_time",
                        discr = discr, quiet = TRUE) |>
      add_parameter_list(!!!param) |>
      add_augmented_list(infec_time = infec_time@D) |>
      add_infec_profile(func = "gamma", discr = discr,
                        shape = shape@P, scale = scale@P)

    expect_equal(evaluate_log_lik(hhmodel),
                 augm_infec_obs_log_lik(house, param, discr))
  }
})

# state of augmented values per iteration
get_augm_state <- function(augm_draws, index) {
  augm_draws[index, , drop = FALSE] |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "name", values_to = "value") |>
    tidyr::separate(name, into = c("id", "augm"),
                    sep = "_", extra = "merge") |>
    tidyr::pivot_wider(names_from = "augm", values_from = "value") |>
    as.data.frame()
}

join_augm_data <- function(house, augm_data) {
  dplyr::left_join(house, augm_data, by = "id")
}

# get total log-likelihood (observation + augmented)
get_log_lik_i <- function(df, par) {
  out <- augm_infec_obs_log_lik(df, par)
  df <- df[df$status > 0, ]
  out +
    sum(dnorm(df$onset_time - df$onset_date, log = TRUE)) +
    sum(dweibull(df$incub_period, shape = 2, scale = 1, log = TRUE))
}

test_that("Parallel augmented-move reduction preserves log-lik invariant", {
  # Exercises move_all_augm_'s RcppParallel reduction (only fired when
  # augmented variables are sampled). We reuse the cheap internal-consistency
  # check (pointwise sum == total) so this stays decoupled from the analytic
  # log-lik comparisons done in the dedicated test below.
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      end_followup = 21)
  house$id <- sprintf("hh01m%s", seq_len(nrow(house)))

  build_hh <- function() {
    new_housepi_model(house, "augmented_infection_time", quiet = TRUE) |>
      add_parameter_template() |>
      add_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
      add_augmented_list(
        infec_time = onset_time@A - incub_period@A,
        onset_time = list(model = "normal", mean = onset_date@D),
        incub_period = list(model = "weibull", shape = 2, scale = 1)
      ) |>
      add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)
  }
  id_draws <- house[house$status > 0, ]$id

  for (t in c(1L, 2L)) {
    results <- suppressWarnings(
      build_hh() |>
        run_inference(steps = 20, thin = 1, adapt = 0, chain = 1, threads = t,
                      progress = FALSE, id_draws = id_draws)
    )
    expect_equal(get_total_log_lik_draws(results)$log_lik,
                 get_log_lik_draws(results)$hh01)
  }
})

test_that("Augmented moves preserve log-likelihood", {
  # dummy dataset
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      end_followup = 21)
  house$id <- sprintf("hh01m%s", seq(1, dim(house)[1]))

  hhmodel <- house |>
    new_housepi_model("augmented_infection_time") |>
    add_parameter_template() |>
    add_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented_list(
      infec_time = onset_time@A - incub_period@A,
      onset_time = list(model = "normal", mean = onset_date@D),
      incub_period = list(model = "weibull", shape = 2, scale = 1)
    ) |>
    add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)

  n_steps <- 100
  id_draws <- house[house$status > 0, ]$id
  results <- suppressWarnings(
    hhmodel |>
      run_inference(steps = n_steps, thin = 1, adapt = 0, chain = 1,
                    progress = FALSE, id_draws = id_draws)
  )

  tot_log_lik <- get_total_log_lik_draws(results)$log_lik
  expect_equal(tot_log_lik, get_log_lik_draws(results)$hh01)

  param_draws <- get_parameter_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()

  augm_draws <- get_augmented_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()

  # constant parameters
  consts <- hhmodel$parameters[hhmodel$parameters$const, c("name", "value")] |>
    tidyr::pivot_wider(names_from = "name", values_from = "value") |>
    as.list()

  # check that snapshot is last state
  snapshot <- hhmodel$export_house_snapshot("hh01")
  expect_identical(unlist(c(param_draws[n_steps, ], consts)),
                   snapshot$param)
  expect_identical(get_augm_state(augm_draws, n_steps),
                   snapshot$augm[house$status > 0, ])

  check_eval <- function(i) {
    df <- join_augm_data(house, get_augm_state(augm_draws, i))
    par <- c(param_draws[i, , drop = FALSE], consts)
    get_log_lik_i(df, par)
  }
  ll_check <- sapply(seq_len(n_steps), check_eval)
  expect_equal(tot_log_lik, ll_check)

  # check if log-likelihood is coherent
  housepi_eval <- function(i) {
    params_i <- c(param_draws[i, , drop = FALSE], consts)
    join_augm_data(house, get_augm_state(augm_draws, i)) |>
      new_housepi_model("augmented_infection_time", quiet = TRUE) |>
      add_augmented_list(
        infec_time = onset_time@A - incub_period@A,
        onset_time = list(value = onset_time@D, const = TRUE,
                          model = "normal", mean = onset_date@D),
        incub_period = list(value = incub_period@D, const = TRUE,
                            model = "weibull", shape = 2, scale = 1)
      ) |>
      add_parameter_list(!!!params_i) |>
      add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P) |>
      evaluate_log_lik()
  }
  expect_equal(tot_log_lik[n_steps], housepi_eval(n_steps))
})
