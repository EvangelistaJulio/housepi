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
  param$beta * (param$hhsize_ref / dim(house)[1]) ^ param$gamma
}

# log-likelihood contribution of escaped members
loglike_escaped <- function(house, param, discr) {
  member_0 <- find_introduction(house)
  t_intro <- member_0$infec_time
  alpha <- param$alpha
  beta <- transmission_rate(house, param)

  output <- 0
  for (k in which(house$status == 0)) {
    # continuous time observes the whole last day: window closes at end_date + 1
    t_end_k <- house[k, "end_date"] + !discr
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

# household log-likelihood of case ascertained model
augm_infec_obs_log_lik <- function(house, param, discr = FALSE) {
  loglike_escaped(house, param, discr) +
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
    # pick introduction randomly with onset_time at 0
    dplyr::mutate(status = NA, start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 0,
                  gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    run_simulation(obs_process = obs_process)
  # there are just introduction cases
  expect_identical(length(unique(simulated$house_id)),
                   sum(simulated$onset_time == 0, na.rm = TRUE))

  # simulate with saturated transmission risk
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 200)) |>
    # pick introduction informed by original data; relative-frame anchor at 0
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100,
                  gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
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
             start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template") |>
    run_inference() |>
    expect_error("at least one infected member")
})

test_that("Sentinel against a reversed follow-up window throws error", {
  # end_date < start_date has no valid exposure period (former study_period < 0)
  house <- data.frame(house_id = "codetalker",
                      status = c(1, 0),
                      infec_time = c(0, NA),
                      start_date = 5, end_date = 2)
  house$id <- sprintf("codetalker_m%s", seq_len(nrow(house)))

  house |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                  hhsize_ref = 4, shape = 2, scale = 1.1) |>
    add_augmented(infec_time = infec_time@D) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P) |>
    evaluate_log_lik() |>
    expect_error("negative follow-up period")
})

test_that("Log-likelihood of case ascertained model is correct", {
  # dummy dataset
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0, 3, 4, NA, NA, NA, NA),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("hh01m%s", seq(1, dim(house)[1]))

  # choice for the parameters
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                hhsize_ref = 4, shape = 2, scale = 1.1)

  for (discr in c(FALSE, TRUE)) {
    hhmodel <- house |>
      new_housepi_model("case_ascertained",
                        discr = discr, quiet = TRUE) |>
      add_parameter(import = param) |>
      add_augmented(infec_time = infec_time@D) |>
      add_infectivity_profile(func = "gamma", discr = discr,
                              shape = shape@P, scale = scale@P)

    evaluate_log_lik(hhmodel) |>
      expect_equal(augm_infec_obs_log_lik(house, param, discr))

    # test complex print, discrete and continuous time
    expect_snapshot(print(hhmodel))
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
                      start_date = 0, end_date = 21)
  house$id <- sprintf("hh01m%s", seq_len(nrow(house)))

  build_hh <- function() {
    house |>
      new_housepi_model("case_ascertained", quiet = TRUE) |>
      add_parameter(import = "template") |>
      add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
      add_augmented(
        infec_time = onset_time - incub_period,
        onset_time = list(model = "normal", mean = onset_date@D),
        incub_period = list(model = "weibull", shape = 2, scale = 1)
      ) |>
      add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)
  }
  id_draws <- house[house$status > 0, ]$id

  for (t in c(1L, 2L)) {
    hhmodel <- build_hh()
    suppressWarnings(
      results <- hhmodel |>
        run_inference(steps = 30, thin = 1, adapt = 20, chain = 1,
                      threads = t, id_draws = id_draws)
    )
    # pointwise log-lik (transmission plus the augmented model, no obs node)
    # matches an independent R recomputation for every retained draw, under any
    # thread count
    param_draws <- get_parameter_draws(results) |>
      posterior::as_draws_matrix() |>
      as.data.frame()
    augm_draws <- get_augmented_draws(results) |>
      posterior::as_draws_matrix() |>
      as.data.frame()
    consts <- hhmodel$parameters |>
      dplyr::filter(const) |>
      dplyr::select(name, value) |>
      tidyr::pivot_wider(names_from = "name", values_from = "value") |>
      as.list()
    ll_ref <- sapply(seq_len(30), function(i) {
      df <- join_augm_data(house, get_augm_state(augm_draws, i))
      get_log_lik_i(df, c(param_draws[i, , drop = FALSE], consts))
    })
    expect_equal(get_log_lik_draws(results)$hh01, ll_ref)
  }
})

test_that("Augmented moves preserve log-likelihood", {
  # dummy dataset
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("hh01m%s", seq(1, dim(house)[1]))

  hhmodel <- house |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented(
      infec_time = onset_time - incub_period,
      onset_time = list(model = "normal", mean = onset_date@D),
      incub_period = list(model = "weibull", shape = 2, scale = 1)
    ) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)

  n_steps <- 100
  id_draws <- house[house$status > 0, ]$id
  results <- suppressWarnings(
    hhmodel |>
      run_inference(steps = n_steps, thin = 1, adapt = 0,
                    chain = 1, id_draws = id_draws)
  )

  # the pointwise log-lik folds transmission and the augmented model; with no
  # obs node it is the full complete-data density (parameter prior aside)
  tot_log_dens <- unname(get_log_lik_draws(results)$hh01)

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
      new_housepi_model("case_ascertained", quiet = TRUE) |>
      add_augmented(
        infec_time = onset_time - incub_period,
        onset_time = list(value = onset_time@D, const = TRUE,
                          model = "normal", mean = onset_date@D),
        incub_period = list(value = incub_period@D, const = TRUE,
                            model = "weibull", shape = 2, scale = 1)
      ) |>
      add_parameter(import = params_i) |>
      add_infectivity_profile("gamma", shape = shape@P, scale = scale@P) |>
      evaluate_log_lik()
  }
  # the fresh-model likelihood matches the tracked pointwise likelihood
  expect_equal(tot_log_dens[n_steps], housepi_eval(n_steps))
})

test_that("Observation node enters the pointwise log-likelihood", {
  house <- data.frame(house_id = "hh01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("hh01m%s", seq_len(nrow(house)))
  id_draws <- house[house$status > 0, ]$id

  hhmodel <- house |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1, max = 14),
      infec_time = list(model = "uniform", min = -14, max = end_date@D)
    ) |>
    add_observation("onset_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)

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

  # the pointwise log_lik combines transmission, observation and the augmented
  # model: it matches a fresh-model evaluation with the augmented values pinned
  ref <- sapply(seq_len(n_steps), function(i) {
    df <- join_augm_data(house, get_augm_state(augm_draws, i))
    par <- c(param_draws[i, , drop = FALSE], consts)
    df |>
      new_housepi_model("case_ascertained", quiet = TRUE) |>
      add_augmented(
        incub_period = list(value = incub_period@D, const = TRUE,
                            model = "gamma", shape = 2, scale = 1, max = 14),
        infec_time = list(value = infec_time@D, const = TRUE,
                          model = "uniform", min = -14, max = end_date@D)
      ) |>
      add_observation("onset_date", model = "normal",
                      mean = infec_time + incub_period, sd = 1) |>
      add_parameter(import = par) |>
      add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P) |>
      evaluate_log_lik()
  })
  expect_equal(get_log_lik_draws(results)$hh01, ref)
})

test_that("Observation through a deterministic node constrains its ancestors", {
  # regression: an observation whose location is a deterministic node
  # (onset_time = infec_time + incub_period) must propagate the constraint to the
  # sampled ancestors. Otherwise moving infec_time / incub_period skips the
  # observation ratio and the latent onset drifts away from the datum.
  set.seed(42) # 'Solid Snake'
  house <- data.frame(house_id = "outerheaven",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      onset_date = c(2, 6, 6, NA, NA, NA, NA),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("outerheaven_m%s", seq_len(nrow(house)))
  id_draws <- house[house$status > 0, ]$id

  results <- suppressWarnings(
    house |>
      new_housepi_model("case_ascertained", quiet = TRUE) |>
      add_parameter(import = "template") |>
      add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
      add_augmented(
        incub_period = list(model = "gamma", shape = 2, scale = 1, max = 14),
        infec_time = list(model = "uniform", min = -14, max = end_date@D),
        onset_time = infec_time + incub_period
      ) |>
      # tightly informative observation reached through the deterministic node
      add_observation(
        "onset_date", model = "normal", mean = onset_time, sd = 0.2
      ) |>
      add_infectivity_profile("gamma", shape = shape@P, scale = scale@P) |>
      run_inference(steps = 400, burn = 100, thin = 1, adapt = 100,
                    chain = 1, id_draws = id_draws)
  )

  onset <- get_augmented_draws(results) |>
    posterior::as_draws_matrix() |>
    as.data.frame()
  # the sampled onset (infec_time + incub_period) tracks the observed onset_date
  recon <- unname(colMeans(onset[, grep("onset_time", names(onset))]))
  expect_lt(max(abs(recon - house$onset_date[house$status > 0])), 0.5)
})

test_that("Discretized observation node initializes finite", {
  # regression: a discretized normal obs whose latent onset starts (from the
  # wide infec_time prior) far below the observed onset_date used to yield a
  # spurious -Inf at chain init (upper-tail CDF-difference cancellation).
  set.seed(101) # 'Liquid Snake'
  house <- data.frame(house_id = "shadowmoses",
                      status = c(1, 1, 1, 1, 0, 0, 0),
                      onset_date = c(2, 5, 6, 7, NA, NA, NA),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("shadowmoses_m%s", seq_len(nrow(house)))
  id_draws <- house[house$status > 0, ]$id

  hhmodel <- house |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1, max = 14),
      infec_time = list(model = "uniform", min = -14, max = end_date@D)
    ) |>
    # discretized observation: the path that previously broke initialization
    add_observation("onset_date", model = "normal", discr = TRUE,
                    loc = infec_time + incub_period) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)

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
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time   = -incub_period
    ) |>
    # one definition set drives simulation too: no hand-written obs_process
    add_observation("onset_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1) |>
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

test_that("Household observation derives from simulated observations", {
  set.seed(1)
  # both the observation (test_date) and the derived column (start_date) are
  # produced internally: neither is provided in the input data
  reps <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 80)) |>
    dplyr::select(-test_date, -start_date) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time   = -incub_period
    ) |>
    add_observation("test_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1) |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    run_simulation(nrep = 2)

  # each repetition is self-consistent: start_date is the per-household na.rm min
  # of *that* repetition's generated test_date (a stale leak would break this),
  # one value per household replicated across rows, never NA (every house has a case)
  for (sim in reps) {
    expect_true(all(c("start_date", "test_date") %in% names(sim)))
    expect_false(any(is.na(sim$start_date)))
    expect_true(all(
      tapply(sim$start_date, sim$house_id, function(x) length(unique(x))) == 1
    ))
    inf <- sim[sim$status > 0, ]
    ref <- tapply(inf$test_date, inf$house_id, min)
    got <- tapply(sim$start_date, sim$house_id, function(x) x[[1]])
    expect_equal(as.numeric(got[order(names(got))]),
                 as.numeric(ref[order(names(ref))]))
  }
})

test_that("Simulation leaves overwritten data columns pristine", {
  set.seed(1)
  m <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% head(unique(house_id), 30)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # test_date is real provided data AND a generated observation: simulation
    # overwrites it with generated values, then must restore the input
    add_observation("test_date", model = "normal",
                    mean = infec_time + incub_period, sd = 1)

  before <- m$export_augmented_data()
  m |> run_simulation(nrep = 2)
  m |> run_simulation(nrep = 2)
  after <- m$export_augmented_data()

  before <- before[order(before$id), ]
  after <- after[order(after$id), ]
  # the overwritten column is restored, and no simulation-added column lingers
  expect_equal(as.numeric(after$test_date), as.numeric(before$test_date))
  expect_setequal(names(after), names(before))
})

test_that("Augmented variable reads a static data-bound household observation", {
  set.seed(1)
  # start_date derived from PROVIDED (static) test_date, read by infec_time's prior
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% head(unique(house_id), 30)) |>
    dplyr::mutate(test_date = 5, end_date = 26) |>
    dplyr::select(-start_date) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      infec_time = list(model = "uniform", min = start_date@D - 5,
                        max = start_date@D + 5)
    ) |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    run_simulation()

  # start_date is materialized statically (= 5) before the index is drawn
  expect_false(any(is.na(simulated$start_date)))
  expect_true(all(simulated$start_date == 5))
})

test_that("Simulation rejects an augmented variable reading a live household observation", {
  # infec_time depends on start_date, which is derived from the GENERATED test_date:
  # an unsatisfiable cycle, rejected with a clear message
  sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% head(unique(house_id), 20)) |>
    dplyr::select(-test_date, -start_date) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      infec_time = list(model = "uniform", min = start_date@D - 5,
                        max = start_date@D + 5)
    ) |>
    add_observation("test_date", model = "normal", mean = infec_time, sd = 1) |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    run_simulation() |>
    expect_error("[Cc]yclic")
})

test_that("Deterministic observation node simulates, inference rejects it", {
  set.seed(1)
  build <- function(data) {
    data |>
      new_housepi_model("case_ascertained") |>
      add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
      add_augmented(
        incub_period = list(model = "gamma", shape = 2, scale = 1),
        infec_time   = -incub_period
      ) |>
      # bare expression (no 'model =' distribution): a deterministic observation
      add_observation(onset_date = floor(infec_time + incub_period))
  }

  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 100)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    build() |>
    run_simulation()

  # generated deterministically for infected members, as a single column
  expect_equal(sum(names(simulated) == "onset_date"), 1)
  inf <- simulated[simulated$status > 0, ]
  expect_equal(inf$onset_date, floor(inf$infec_time + inf$incub_period))
  expect_true(all(is.na(simulated$onset_date[simulated$status == 0])))

  # inference rejects the degenerate indicator likelihood with a clear message
  evaluate_log_lik(build(simulated)) |>
    expect_error("Deterministic observation node")
})

test_that("Household observations aggregate a data column across members", {

  houses <- head(unique(sarscov2_2021_icofs$house_id), 60)
  # force one household to have all-NA test_date: the na.rm reduction is empty.
  # the reference below reads 'dat', so mutate it here to keep both in sync.
  na_house <- houses[[1]]
  dat <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% houses) |>
    dplyr::mutate(
      start_date = 0,
      test_date = dplyr::if_else(house_id == na_house, as.Date(NA), test_date)
    )

  simulated <- dat |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    add_observation(
      first_test = min(test_date@D),
      last_test = max(test_date@D),
      level = "house"
    ) |>
    run_simulation()

  expect_true(all(c("first_test", "last_test") %in% names(simulated)))

  # one value per household, replicated across every member row
  expect_true(all(
    tapply(simulated$first_test, simulated$house_id,
           function(x) length(unique(x))) == 1
  ))

  # values equal the na.rm min/max of test_date over each household's members,
  # NA when every member is missing (the all-NA household)
  ref <- function(f) {
    v <- tapply(as.numeric(dat$test_date), dat$house_id,
                function(x) if (all(is.na(x))) NA_real_ else f(x, na.rm = TRUE))
    v[order(names(v))]
  }
  got <- function(col) {
    v <- tapply(simulated[[col]], simulated$house_id, function(x) x[[1]])
    v[order(names(v))]
  }
  expect_equal(got("first_test"), ref(min))
  expect_equal(got("last_test"), ref(max))

  # singleton household (one member): the aggregate is that member's value
  sizes <- table(dat$house_id)
  singleton <- names(sizes)[sizes == 1]
  if (length(singleton) > 0) {
    one <- singleton[[1]]
    expect_equal(unique(simulated$first_test[simulated$house_id == one]),
                 as.numeric(dat$test_date[dat$house_id == one]))
  }

  # the all-NA household yields NA across its rows
  expect_true(all(is.na(simulated$first_test[simulated$house_id == na_house])))
})

test_that("Data-bound household observation supplies a required data column", {
  base <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% head(unique(house_id), 40)) |>
    dplyr::select(-start_date)

  # derive the required 'start_date' from data instead of providing the column
  build_derived <- function() {
    base |>
      new_housepi_model("case_ascertained") |>
      add_parameter(import = "template") |>
      add_augmented(import = "template") |>
      add_observation(start_date = min(test_date@D), level = "house")
  }
  # reference: the same quantity provided as an ordinary data column
  build_provided <- function() {
    base |>
      dplyr::group_by(house_id) |>
      dplyr::mutate(start_date = min(test_date, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      new_housepi_model("case_ascertained") |>
      add_parameter(import = "template") |>
      add_augmented(import = "template")
  }

  # deriving the required column lets inference run without it in the data
  set.seed(7)
  derived <- build_derived()
  expect_silent(derived |> evaluate_log_lik())

  # the materialized column equals the per-household na.rm min of test_date
  # (start_date is a derived numeric column; test_date round-trips as a Date)
  aug <- derived$export_augmented_data()
  ok <- tapply(seq_len(nrow(aug)), aug$house_id, function(ix) {
    isTRUE(all.equal(unique(aug$start_date[ix]),
                     min(as.numeric(aug$test_date[ix]), na.rm = TRUE)))
  })
  expect_true(all(unlist(ok)))

  # identical log-likelihood to providing start_date as a data column
  set.seed(13)
  ll_derived <- build_derived() |> evaluate_log_lik()
  set.seed(13)
  ll_provided <- build_provided() |> evaluate_log_lik()
  expect_equal(ll_derived, ll_provided)
})

test_that("Data-bound household observation rejects conflict and latent dependence", {
  base <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% head(unique(house_id), 20)) |>
    dplyr::select(-start_date)

  # conflict data and derived household observation
  base |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template") |>
    add_observation(start_date = min(test_date@D), level = "house") |>
    evaluate_log_lik() |>
    expect_error("already provided")

  # forbidden household observation depending on parameters
  base |>
    new_housepi_model("case_ascertained") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template") |>
    add_observation(start_date = min(test_date@D) + alpha@P, level = "house") |>
    evaluate_log_lik() |>
    expect_error("pure function of data")
})

test_that("add_status rejects an empty probability list", {
  new_housepi_model(model = "test") |>
    add_status() |>
    expect_error("[Nn]o status probability")
})

test_that("Status process generates a labelled status factor in simulation", {
  set.seed(42) # 'Solid Snake'
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 150)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # categorical drawn first; labels (declaration order) name the statuses
    add_status(symptomatic = 0.2, asymptomatic = 0.8) |>
    run_simulation()

  # status is a factor: code 0 -> escaped, 1 -> symptomatic, 2 -> asymptomatic
  expect_s3_class(simulated$status, "factor")
  expect_identical(levels(simulated$status),
                   c("escaped", "symptomatic", "asymptomatic"))
  inf <- as.character(simulated$status[simulated$status != "escaped"])
  # outcome split follows the probabilities (0.2 symptomatic)
  expect_lt(abs(mean(inf == "symptomatic") - 0.2), 0.08)
})

test_that("Status process drives status-specific augmented blocks", {
  set.seed(101) # 'Liquid Snake'
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 200)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    # status-specific block: symptomatic incubation is longer (character status=)
    add_augmented(status = "symptomatic",
                  incub_period = list(model = "gamma", shape = 8, scale = 1)) |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  by_status <- tapply(simulated$incub_period,
                      as.character(simulated$status), mean, na.rm = TRUE)
  # the resolved status selects the block: symptomatic ~ shape 8, asymptomatic ~ 2
  expect_gt(by_status[["symptomatic"]], 5)
  expect_lt(by_status[["asymptomatic"]], 4)
})

test_that("Status is a constant in inference", {
  # refined statuses already in the data; infec_time is read from data so the
  # density is deterministic and the only difference is the status process.
  house <- data.frame(house_id = "outerheaven",
                      status = c(1, 2, 1, 0, 0, 0, 0),
                      infec_time = c(0, 3, 4, NA, NA, NA, NA),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("outerheaven_m%s", seq_len(nrow(house)))
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                hhsize_ref = 4, shape = 2, scale = 1.1)

  build <- function(with_status) {
    m <- house |>
      new_housepi_model("case_ascertained", quiet = TRUE) |>
      add_parameter(import = param) |>
      add_augmented(infec_time = infec_time@D) |>
      add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)
    if (with_status) {
      m <- m |> add_status(symptomatic = 0.7, asymptomatic = 0.3)
    }
    m
  }

  # the status process drops out of the density (status is fixed data)
  expect_equal(evaluate_log_lik(build(TRUE)), evaluate_log_lik(build(FALSE)))
})

test_that("Status process rejects latent/parameter probabilities", {
  house <- data.frame(house_id = "zanzibar",
                      status = c(1, 2, 0, 0),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("zanzibar_m%s", seq_len(nrow(house)))

  house |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_augmented(import = "template") |>
    # probabilities must be a pure function of data/numerics
    add_status(symptomatic = alpha@P, asymptomatic = 3) |>
    evaluate_log_lik() |>
    expect_error("only on data or numeric")
})

test_that("Character status round-trips with stable codes", {
  set.seed(202) # 'Gray Fox'
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 120)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained") |>
    add_parameter(alpha = 0, beta = 100, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  # feed the factor-status output back for inference: factorization maps the
  # labels to the same codes (escaped -> 0), so the density is finite
  reinfer <- simulated |>
    dplyr::select(id, house_id, status, infec_time, start_date, end_date) |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                  hhsize_ref = 4, shape = 2, scale = 1.1) |>
    add_augmented(infec_time = infec_time@D) |>
    add_infectivity_profile(func = "gamma", shape = shape@P, scale = scale@P)
  expect_true(is.finite(evaluate_log_lik(reinfer)))

  # the data-derived dictionary matches the simulated infected labels
  gen <- reinfer$.__enclos_env__$private$generator_
  expect_identical(gen$status_data_labels, c("symptomatic", "asymptomatic"))
})

test_that("escaped label is configurable", {
  set.seed(303) # 'Big Boss'
  simulated <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% sample(unique(house_id), 60)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case_ascertained", escaped = "healthy") |>
    add_parameter(alpha = 0, beta = 0.001, gamma = 1, hhsize_ref = 4) |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    ) |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    run_simulation()

  # weak transmission leaves escapeds, labelled per the construction argument
  expect_true("healthy" %in% levels(simulated$status))
  expect_true(any(as.character(simulated$status) == "healthy"))
})

test_that("Escaped label colliding with an infected status is rejected", {
  house <- data.frame(house_id = "bigshell",
                      status = c(1, 2, 0, 0),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("bigshell_m%s", seq_len(nrow(house)))

  house |>
    new_housepi_model("case_ascertained", quiet = TRUE,
                      escaped = "symptomatic") |>
    add_parameter(import = "template") |>
    add_augmented(import = "template") |>
    add_status(symptomatic = 0.5, asymptomatic = 0.5) |>
    expect_error("collides with the escaped-status label")
})

test_that("Mixing integer and character status declarations is rejected", {
  house <- data.frame(house_id = "shadowmoses",
                      status = c(1, 2, 0, 0),
                      start_date = 0, end_date = 21)
  house$id <- sprintf("shadowmoses_m%s", seq_len(nrow(house)))

  # one call carries both an integer and a character status: caught at the wall
  house |>
    new_housepi_model("case_ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_augmented(import = "template") |>
    add_augmented(
      incub_period = list(model = "gamma", shape = 2, scale = 1,
                          status = "symptomatic"),
      onset = list(value = 1, status = 1)
    ) |>
    evaluate_log_lik() |>
    expect_error("[Mm]ix")
})

test_that("Model-required data can be postponed", {
  # drop a model-required data column: it may be supplied later
  # construction must tolerate the gap silently
  dat <- sarscov2_2021_icofs |>
    dplyr::filter(house_id %in% head(unique(house_id), 20))
  dat$start_date <- NULL

  hhmodel <- dat |>
    new_housepi_model("case_ascertained")

  # hitting a wall (evaluate_*/run_*) raises capital error
  hhmodel |>
    evaluate_log_lik() |>
    expect_error("start_date")
  hhmodel |>
    run_simulation() |>
    expect_error("start_date")
})
