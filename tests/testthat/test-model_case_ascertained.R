library(testthat)

# Over budget: the likelihood test builds the discrete-time model, a DLL
# compile the infectivity-profile tests then find cached.

# household of three cases and four escapees, the introduction at day 0
mk_house <- helper_onset_house

# simulation over real households, the transmission risk left to each test
mk_simulation <- local({
  hhmodel <- sarscov2_2021_icofs |>
    dplyr::filter(hhid %in% utils::head(unique(hhid), 100)) |>
    dplyr::mutate(start_date = 0, end_date = 21) |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_latent(
      incub_period = list(model = "gamma", shape = 2, scale = 1),
      infec_time = -incub_period
    )
  list(model = hhmodel)
})

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
  infected <- house[house$status > 0, ]
  infected[infected$infec_time == min(infected$infec_time), ]
}

# transmission rate
transmission_rate <- function(house, param) {
  param$beta * (param$hhsize_ref / dim(house)[1]) ^ param$gamma
}

# log-likelihood contribution of escaped members
loglike_escaped <- function(house, param, discr) {
  t_intro <- find_introduction(house)$infec_time
  alpha <- param$alpha
  beta <- transmission_rate(house, param)

  output <- 0
  for (k in which(house$status == 0)) {
    # continuous time observes the whole last day: window closes at end_date + 1
    t_end_k <- house[k, "end_date"] + !discr
    roi_cumul <- 0
    alpha_cumul <- alpha * (t_end_k - t_intro)
    for (i in which(house$status == 1)) {
      exposure <- t_end_k - house[i, "infec_time"]
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
        gentime <- t_infec_i - house[k, "infec_time"]
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
        output <- output + log(-expm1(-roi_inst_i)) - roi_cumul_i + alpha
      } else {
        output <- output + log(roi_inst_i) - roi_cumul_i
      }
    }
  }

  output
}

# household log-likelihood of case ascertained model
latent_infec_obs_log_lik <- function(house, param, discr = FALSE) {
  loglike_escaped(house, param, discr) +
    loglike_infectious(house, param, discr)
}

# model reading the infection times of `hhdata`, every parameter pinned
build_onset_model <- function(hhdata, param, discr) {
  # nolint start: object_usage_linter. DSL symbols, read by the encoder
  transform(hhdata, infec_date = infec_time) |>
    new_housepi_model("case-ascertained", discr = discr, quiet = TRUE) |>
    add_parameter(import = param) |>
    add_latent(infec_time = infec_date) |>
    add_infectivity_profile(func = "gamma", discr = discr,
                            shape = shape, scale = scale)
  # nolint end
}

# state of latent values per iteration
get_latent_state <- function(latent_draws, index) {
  latent_draws[index, , drop = FALSE] |>
    tidyr::pivot_longer(dplyr::everything(),
                        names_to = "name", values_to = "value") |>
    tidyr::separate("name", into = c("id", "latent"),
                    sep = "_", extra = "merge") |>
    tidyr::pivot_wider(names_from = "latent", values_from = "value") |>
    as.data.frame()
}

##################################################
## CASE ASCERTAINED HH MODEL : SIMULATION       ##
##################################################

test_that("Without transmission only the introductions are infected", {
  report_onset <- function(hhdata) {
    hhdata$onset_time <- hhdata$infec_time + hhdata$incub_period
    hhdata
  }
  # the introduction is the only case, its onset_time at 0
  simulated <- mk_simulation$model |>
    add_parameter(alpha = 1e-12, beta = 0, gamma = 1, hhsize_ref = 4) |>
    add_reporting_process(report_onset) |>
    run_simulation()
  length(unique(simulated$hhid)) |>
    expect_identical(sum(simulated$onset_time == 0, na.rm = TRUE))
})

test_that("A saturated transmission infects every member", {
  simulated <- mk_simulation$model |>
    add_parameter(alpha = 1e-12, beta = 100, gamma = 1, hhsize_ref = 4) |>
    run_simulation()

  (simulated$status > 0) |> expect_all_true()
  simulated$infec_time |> expect_no_na()
  simulated$incub_period |> expect_no_na()
})

##################################################
## CASE ASCERTAINED HH MODEL : LOG-LIKELIHOOD   ##
##################################################

test_that("Sentinel against no introduction throws error", {
  data.frame(id = seq(2), hhid = "h1", status = 0,
             start_date = 0, end_date = 21) |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    run_inference() |>
    expect_error("at least one infected member")
})

test_that("Sentinel against a reversed follow-up window throws error", {
  # end_date < start_date has no valid exposure period
  data.frame(id = c("m1", "m2"), hhid = "h1",
             status = c(1, 0), infec_date = c(0, NA),
             start_date = 5, end_date = 2) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(import = mk_house$param) |>
    add_latent(infec_time = infec_date) |>
    add_infectivity_profile(func = "gamma", shape = shape, scale = scale) |>
    evaluate_log_lik() |>
    expect_error("negative follow-up period")
})

test_that("Sentinel against a member missing a follow-up date throws error", {
  hhdata <- transform(mk_house$data, infec_date = infec_time)
  hhmodel <- new_housepi_model(model = "case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(infec_time = infec_date) |>
    add_infectivity_profile(import = "template")

  # escaped members without end_date throw
  hhmodel |>
    add_data_wide(transform(hhdata, end_date = c(rep(21, 6), NA))) |>
    evaluate_log_lik() |>
    expect_error("[Dd]ata .end_date. has no value for individual .m7.")

  # end_date computed from member observations have no value for escaped members
  hhmodel |>
    add_data_wide(dplyr::select(hhdata, -end_date)) |>
    add_observation(onset_date = list(model = "normal", mean = infec_time),
                    end_date = onset_date + 21) |>
    evaluate_log_lik() |>
    expect_error("[Dd]ata.*end_date. needs a value for every member")
})

test_that("Log-likelihood of case ascertained model is correct", {
  house <- mk_house$data[c("id", "hhid", "status", "infec_time",
                           "start_date", "end_date")]
  # risks far below the machine epsilon
  tiny <- modifyList(mk_house$param, list(alpha = 1e-20, beta = 1e-20))

  for (discr in c(FALSE, TRUE)) {
    hhmodel <- build_onset_model(house, mk_house$param, discr)
    evaluate_log_lik(hhmodel) |>
      expect_equal(latent_infec_obs_log_lik(house, mk_house$param, discr))

    build_onset_model(house, tiny, discr) |>
      evaluate_log_lik() |>
      expect_equal(latent_infec_obs_log_lik(house, tiny, discr))

    # test complex print, discrete and continuous time
    print(hhmodel) |> expect_snapshot()
  }
})

test_that("A case with no source of infection has zero likelihood", {
  # `m2` is the only child, out of reach of the community,
  # and no household transmission
  hhdata <- transform(mk_house$data[c("id", "hhid", "status", "infec_time",
                                      "start_date", "end_date")],
                      age = c("adult", "child", rep("adult", 5)))
  strata <- data.frame(age = c(NA, "child"),
                       param_name = c(".reference", "com_child"),
                       value = c(NA, 0), const = c(NA, TRUE))
  param <- modifyList(mk_house$param, list(beta = 0))

  for (discr in c(FALSE, TRUE)) {
    build_onset_model(hhdata, param, discr) |>
      add_group_risk("community", strata = strata) |>
      evaluate_log_lik() |>
      expect_identical(-Inf)
  }
})

test_that("Latent moves preserve log-likelihood", {
  house <- mk_house$data[c("id", "hhid", "status", "onset_date",
                           "start_date", "end_date")]
  hhmodel <- house |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_parameter(import = "template") |>
    add_parameter(shape = 2, scale = 1.1, append = TRUE) |>
    add_latent(
      infec_time = list(min = -14, max = end_date),
      incub_period = list(model = "weibull", shape = 2, scale = 1)
    ) |>
    # read by the incubation only: a move of the infection time skips it
    add_observation("onset_date", model = "normal",
                    mean = incub_period, sd = 1) |>
    add_infectivity_profile(func = "gamma", shape = shape, scale = scale)

  # the latent moves are reduced over two threads
  n_steps <- 30L
  results <- hhmodel |>
    run_inference(steps = n_steps, thin = 1L, adapt = 0L, chains = 1L,
                  threads = 2L, id_draws = house$id[house$status > 0]) |>
    suppressWarnings()

  # the pointwise log-lik folds transmission, the latent model and the
  # observation node: the complete-data density, parameter prior aside
  # the flat infection time adds no density: transmission gives it
  log_lik <- unname(get_log_lik_draws(results)$h1)
  param_draws <- get_parameter_draws(results) |>
    draws_to_frame()
  latent_draws <- get_latent_draws(results) |>
    draws_to_frame()
  consts <- get_parameters(hhmodel) |>
    subset(const, c("name", "value")) |>
    tidyr::pivot_wider(names_from = "name", values_from = "value") |>
    as.list()
  state_at <- function(i) {
    dplyr::left_join(house, get_latent_state(latent_draws, i), by = "id")
  }

  # the snapshot holds the last state
  snapshot <- inspect_house_state(hhmodel, "h1")
  unlist(c(param_draws[n_steps, ], consts)) |>
    expect_identical(snapshot$parameters)
  get_latent_state(latent_draws, n_steps) |>
    expect_identical(snapshot$latent[house$status > 0, ])

  # every retained draw matches an independent R recomputation
  expected <- vapply(seq_len(n_steps), function(i) {
    hhdata <- state_at(i)
    par <- c(param_draws[i, , drop = FALSE], consts)
    infected <- hhdata[hhdata$status > 0, ]
    latent_infec_obs_log_lik(hhdata, par) +
      sum(dnorm(infected$onset_date, infected$incub_period, log = TRUE)) +
      sum(dweibull(infected$incub_period, shape = 2, scale = 1, log = TRUE))
  }, numeric(1))
  log_lik |> expect_equal(expected)

  # the fresh-model likelihood matches the tracked pointwise likelihood
  state_at(n_steps) |>
    transform(infec_draw = infec_time, incub_draw = incub_period) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_latent(
      infec_time = list(value = infec_draw, const = TRUE),
      incub_period = list(value = incub_draw, const = TRUE,
                          model = "weibull", shape = 2, scale = 1)
    ) |>
    add_observation("onset_date", model = "normal",
                    mean = incub_period, sd = 1) |>
    add_parameter(import = c(param_draws[n_steps, , drop = FALSE], consts)) |>
    add_infectivity_profile("gamma", shape = shape, scale = scale) |>
    evaluate_log_lik() |>
    expect_equal(log_lik[n_steps])
})
