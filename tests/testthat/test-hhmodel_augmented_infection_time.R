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
    t_end_k <- house[k, "start_date"] + house[k, "study_period"]
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
  obs_process <- function(df) {
    df$onset_time <- df$infec_time + df$incub_period
    df
  }
  # simulate with no transmission risk
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 200)) |>
    # pick introduction randomly
    dplyr::mutate(status = NA) |>
    new_housepi_model("augmented_infection_time") |>
    add_parameter_list(alpha = 0, beta = 0,
                       gamma = 1, size_ref = 4) |>
    add_augmented_list(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    run_simulation(obs_process = obs_process)
  # there are just introduction cases
  expect_identical(length(unique(simulated$house_id)),
                   sum(simulated$onset_time == 0, na.rm = TRUE))

  # simulate with saturated transmission risk
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 200)) |>
    # pick introduction informed by original data
    new_housepi_model("augmented_infection_time") |>
    add_parameter_list(alpha = 0, beta = 100,
                       gamma = 1, size_ref = 4) |>
    add_augmented_list(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    run_simulation(obs_process = obs_process)
  # all infected
  expect_true(all(simulated$status > 0))
  expect_false(any(is.na(simulated$infec_time)))
  expect_equal(simulated$infec_time,
               simulated$onset_time - simulated$incub_period)
})

test_that("Sentinel against no introduction throws error", {
  data.frame(id = seq(2), house_id = "h1", status = 0,
             start_date = 0, study_period = 21) |>
    new_housepi_model("augmented_infection_time") |>
    add_parameter_template() |>
    add_augmented_template() |>
    run_inference() |>
    expect_error(regexp = "at least one infected member")
})

test_that("Log-likelihood of augmented infection time model is correct", {
  # dummy dataset
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0, 3, 4, NA, NA, NA, NA),
                      start_date = 0, study_period = 21)
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

    expect_equal(eval_log_density(hhmodel),
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
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      start_date = 0, study_period = 21)
  house$id <- sprintf("hh01m%s", seq_len(nrow(house)))

  build_hh <- function() {
    house |>
      new_housepi_model("augmented_infection_time", quiet = TRUE) |>
      add_parameter_template() |>
      add_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
      add_augmented_list(
        infec_time = onset_time@A - incub_period@A,
        onset_time = list(prior = "normal", mean = onset_date@D),
        incub_period = list(prior = "weibull", shape = 2, scale = 1)
      ) |>
      add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)
  }
  id_draws <- house[house$status > 0, ]$id

  for (t in c(1L, 2L)) {
    hhmodel <- build_hh()
    suppressWarnings(
      results <- hhmodel |>
        run_inference(steps = 30, thin = 1, adapt = 20, chain = 1,
                      threads = t, id_draws = id_draws)
    )
    # pointwise log-lik (transmission only, no obs node) matches an independent
    # R recomputation for every retained draw, under any thread count
    param_draws <- get_parameter_draws(results) |>
      posterior::as_draws_matrix() |>
      as.data.frame()
    augm_draws <- get_augmented_draws(results) |>
      posterior::as_draws_matrix() |>
      as.data.frame()
    consts <- hhmodel$parameters[hhmodel$parameters$const, c("name", "value")] |>
      tidyr::pivot_wider(names_from = "name", values_from = "value") |>
      as.list()
    transm_ref <- sapply(seq_len(30), function(i) {
      df <- join_augm_data(house, get_augm_state(augm_draws, i))
      augm_infec_obs_log_lik(df, c(param_draws[i, , drop = FALSE], consts))
    })
    expect_equal(get_log_lik_draws(results)$hh01, transm_ref)
  }
})

test_that("Augmented moves preserve log-likelihood", {
  # dummy dataset
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      start_date = 0, study_period = 21)
  house$id <- sprintf("hh01m%s", seq(1, dim(house)[1]))

  hhmodel <- house |>
    new_housepi_model("augmented_infection_time") |>
    add_parameter_template() |>
    add_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented_list(
      infec_time = onset_time@A - incub_period@A,
      onset_time = list(prior = "normal", mean = onset_date@D),
      incub_period = list(prior = "weibull", shape = 2, scale = 1)
    ) |>
    add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)

  n_steps <- 100
  id_draws <- house[house$status > 0, ]$id
  results <- suppressWarnings(
    hhmodel |>
      run_inference(steps = n_steps, thin = 1, adapt = 0,
                    chain = 1, id_draws = id_draws)
  )

  # no obs node here, so the pointwise log-lik is the transmission likelihood
  tot_log_lik <- get_log_lik_draws(results)$hh01
  # full log-density recomposes the likelihood + augmented priors;
  # augmented prior = joint lprior minus the parameter log-prior. The
  # per-parameter prior is internal, so reach the private storage to sum it.
  param_prior <- results$.__enclos_env__$private$param_prior_draws_
  param_lprior <- if (is.null(param_prior)) 0 else
    rowSums(posterior::as_draws_matrix(param_prior))
  augm_prior <- get_lprior_draws(results)$lprior - param_lprior
  tot_log_dens <- unname(tot_log_lik + augm_prior)

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
  expect_equal(tot_log_dens, ll_check)

  # check if log-likelihood is coherent
  housepi_eval <- function(i) {
    params_i <- c(param_draws[i, , drop = FALSE], consts)
    join_augm_data(house, get_augm_state(augm_draws, i)) |>
      new_housepi_model("augmented_infection_time", quiet = TRUE) |>
      add_augmented_list(
        infec_time = onset_time@A - incub_period@A,
        onset_time = list(value = onset_time@D, const = TRUE,
                          prior = "normal", mean = onset_date@D),
        incub_period = list(value = incub_period@D, const = TRUE,
                            prior = "weibull", shape = 2, scale = 1)
      ) |>
      add_parameter_list(!!!params_i) |>
      add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P) |>
      eval_log_density()
  }
  # the fresh-model likelihood matches the tracked pointwise likelihood
  expect_equal(tot_log_lik[n_steps], housepi_eval(n_steps))
})

test_that("Observation node enters the pointwise log-likelihood", {
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      start_date = 0, study_period = 21)
  house$id <- sprintf("hh01m%s", seq_len(nrow(house)))
  id_draws <- house[house$status > 0, ]$id

  hhmodel <- house |>
    new_housepi_model("augmented_infection_time", quiet = TRUE) |>
    add_parameter_template() |>
    add_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented_list(
      incub_period = list(prior = "gamma", shape = 2, scale = 1, max = 14),
      infec_time   = list(prior = "uniform", min = -14, max = study_period@D,
                          proposal = "normal")
    ) |>
    add_observation("onset_date", obs = "normal",
                    mean = infec_time@A + incub_period@A, sd = 1) |>
    add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)

  n_steps <- 100
  results <- suppressWarnings(
    hhmodel |> run_inference(steps = n_steps, thin = 1, adapt = 0,
                             chain = 1, id_draws = id_draws)
  )

  param_draws <- get_parameter_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()
  augm_draws <- get_augmented_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()
  consts <- hhmodel$parameters[hhmodel$parameters$const, c("name", "value")] |>
    tidyr::pivot_wider(names_from = "name", values_from = "value") |>
    as.list()

  # the pointwise log_lik combines transmission and observation: it matches an
  # independent recomputation of both pieces per retained draw
  ref <- sapply(seq_len(n_steps), function(i) {
    df <- join_augm_data(house, get_augm_state(augm_draws, i))
    transm <- augm_infec_obs_log_lik(df, c(param_draws[i, , drop = FALSE],
                                           consts))
    dd <- df[df$status > 0, ]
    obs <- sum(dnorm(dd$onset_date, mean = dd$infec_time + dd$incub_period,
                     sd = 1, log = TRUE))
    transm + obs
  })
  expect_equal(get_log_lik_draws(results)$hh01, ref)
})

test_that("Discretized observation node initializes finite", {
  # regression: a discretized normal obs whose latent onset starts (from the
  # wide infec_time prior) far below the observed onset_date used to yield a
  # spurious -Inf at chain init (upper-tail CDF-difference cancellation).
  set.seed(101) # 'Liquid Snake'
  house <- data.frame(house_id = "shadowmoses",
                      status = c(1, 1, 1, 1, 0, 0, 0),
                      onset_date = c(2, 5, 6, 7, NA, NA, NA),
                      start_date = 0, study_period = 21)
  house$id <- sprintf("shadowmoses_m%s", seq_len(nrow(house)))
  id_draws <- house[house$status > 0, ]$id

  hhmodel <- house |>
    new_housepi_model("augmented_infection_time", quiet = TRUE) |>
    add_parameter_template() |>
    add_parameter_list(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented_list(
      incub_period = list(prior = "gamma", shape = 2, scale = 1, max = 14),
      infec_time   = list(prior = "uniform", min = -14, max = study_period@D,
                          proposal = "normal")
    ) |>
    # discretized observation: the path that previously broke initialization
    add_observation("onset_date", obs = "normal", discr = TRUE,
                    loc = infec_time@A + incub_period@A) |>
    add_infec_profile(func = "gamma", shape = shape@P, scale = scale@P)

  # initialization no longer returns infinity; the chain runs and the
  # complete-data log-likelihood stays finite throughout
  results <- expect_no_error(suppressWarnings(
    hhmodel |> run_inference(steps = 50, thin = 1, adapt = 0,
                             chain = 1, id_draws = id_draws)
  ))
  expect_true(all(is.finite(get_log_lik_draws(results)$shadowmoses)))
})

test_that("Observation node generates the data column in simulation", {
  set.seed(1)
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 100)) |>
    new_housepi_model("augmented_infection_time") |>
    add_parameter_list(alpha = 0, beta = 100, gamma = 1, size_ref = 4) |>
    add_augmented_list(
      incub_period = list(prior = "gamma", shape = 2, scale = 1),
      infec_time   = -incub_period
    ) |>
    # one definition set drives simulation too: no hand-written obs_process
    add_observation("onset_date", obs = "normal",
                    mean = infec_time@A + incub_period@A, sd = 1) |>
    run_simulation()

  # onset_date is generated for infected members and only those
  expect_true("onset_date" %in% names(simulated))
  expect_false(any(is.na(simulated$onset_date[simulated$status > 0])))
  expect_true(all(is.na(simulated$onset_date[simulated$status == 0])))

  # generated onset_date is a draw from normal(infec_time + incub_period, 1)
  inf <- simulated[simulated$status > 0, ]
  resid <- inf$onset_date - (inf$infec_time + inf$incub_period)
  expect_lt(abs(mean(resid)), 0.2)
  expect_lt(abs(sd(resid) - 1), 0.2)
})

test_that("Observation are reserved to add_observation", {
  hhmodel <- data.frame(id = seq(2), house_id = "h1", status = 0,
                        start_date = 0, study_period = 21) |>
    new_housepi_model("augmented_infection_time")

  # ensure augmented decorators reject 'obs ='
  hhmodel |>
    add_augmented("onset_date", obs = "normal") |>
    expect_error("add_observation")
  hhmodel |>
    add_augmented_list(onset_date = list(obs = "normal")) |>
    expect_error("add_observation")
})

test_that("Deterministic observation node simulates, inference rejects it", {
  set.seed(1)
  build <- function(data) {
    data |>
      new_housepi_model("augmented_infection_time") |>
      add_parameter_list(alpha = 0, beta = 100, gamma = 1, size_ref = 4) |>
      add_augmented_list(
        incub_period = list(prior = "gamma", shape = 2, scale = 1),
        infec_time   = -incub_period
      ) |>
      # bare expression (no 'obs =' distribution): a deterministic observation
      add_observation_list(onset_date = floor(infec_time + incub_period))
  }

  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 100)) |>
    build() |>
    run_simulation()

  # generated deterministically for infected members, as a single column
  expect_equal(sum(names(simulated) == "onset_date"), 1)
  inf <- simulated[simulated$status > 0, ]
  expect_equal(inf$onset_date, floor(inf$infec_time + inf$incub_period))
  expect_true(all(is.na(simulated$onset_date[simulated$status == 0])))

  # inference rejects the degenerate indicator likelihood with a clear message
  expect_error(evaluate_log_lik(build(simulated)),
               regexp = "Deterministic observation node")
})
