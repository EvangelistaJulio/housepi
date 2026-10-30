library(testthat)

####  Helpers  ####

# Over budget: the likelihood test builds the discrete-time model, a DLL
# compile the infectivity-profile tests then find cached.

# household of three cases and four escapees, the introduction at day 0
mk_house <- helper_onset_house

# simulation over real households, the transmission rate left to each test
mk_simulation <- local({
  hh_model <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 100)) |>
    dplyr::mutate(end_date = 21) |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template")
  list(model = hh_model)
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
  param$beta * (param$hh_size_ref / dim(house)[1]) ^ param$gamma
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

# model reading the infection times of `hh_data`, every parameter pinned
build_onset_model <- function(hh_data, param, discr) {
  # nolint start: object_usage_linter. DSL symbols, read by the encoder
  transform(hh_data, infec_date = infec_time) |>
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

####  Simulation  ####

test_that("With no transmission only a random household member is infected", {
  # the introduction is the only case, infected at the time origin
  # transmission generates `infec_time`: simulation needs no definition
  hh_model <- mk_simulation$model |>
    add_parameter(alpha = 1e-12, beta = 0, gamma = 1, hh_size_ref = 4)
  simulated <- hh_model |> run_simulation()
  infected <- simulated$status > 0
  as.vector(table(simulated$hh_id[infected])) |> expect_all_equal(1L)
  simulated$infec_time[infected] |> expect_all_equal(0)

  # inference fails without definition
  get_latent(hh_model) |> expect_nrow(0L)
  hh_model |>
    evaluate_log_lik() |>
    expect_error("[Mm]issing latent definition.*infec_time")

  # definitions stay unchanged
  # over 100 households the introduction is not always the same member
  member <- ave(seq_along(simulated$hh_id), simulated$hh_id, FUN = seq_along)
  unique(member[infected]) |> length() |> expect_gt(1L)
})

test_that("A saturated transmission infects every member", {
  simulated <- mk_simulation$model |>
    add_parameter(alpha = 1e-12, beta = 100, gamma = 1, hh_size_ref = 4) |>
    run_simulation()

  (simulated$status > 0) |> expect_all_true()
  simulated$infec_time |> expect_no_na()
})

####  Log-likelihood  ####

test_that("Sentinel against no introduction throws error", {
  data.frame(id = seq(2), hh_id = "h1", status = 0,
             end_date = 21) |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template") |>
    add_parameter(import = "template") |>
    add_latent(import = "template") |>
    run_inference() |>
    expect_error("at least one infected member")
})

test_that("Sentinel against a member missing a follow-up date throws error", {
  hh_data <- transform(mk_house$data, infec_date = infec_time)
  hh_model <- new_housepi_model(model = "case-ascertained") |>
    add_parameter(import = "template") |>
    add_latent(infec_time = infec_date) |>
    add_infectivity_profile(import = "template")

  # escaped members without end_date throw
  hh_model |>
    add_dataset(transform(hh_data, end_date = c(rep(21, 6), NA))) |>
    evaluate_log_lik() |>
    expect_error("[Dd]ata .end_date. has no value for individual .m7.")

  # end_date computed from member observations have no value for escaped members
  hh_model |>
    add_dataset(dplyr::select(hh_data, -end_date)) |>
    add_observation(onset_date = list(model = "normal", mean = infec_time),
                    end_date = onset_date + 21) |>
    evaluate_log_lik() |>
    expect_error("[Dd]ata.*end_date. needs a value for every member")
})

test_that("Log-likelihood of case ascertained model is correct", {
  house <- mk_house$data[c("id", "hh_id", "status", "infec_time", "end_date")]
  # rates far below the machine epsilon
  tiny <- modifyList(mk_house$param, list(alpha = 1e-20, beta = 1e-20))

  for (discr in c(FALSE, TRUE)) {
    hh_model <- build_onset_model(house, mk_house$param, discr)
    evaluate_log_lik(hh_model) |>
      expect_equal(latent_infec_obs_log_lik(house, mk_house$param, discr))

    build_onset_model(house, tiny, discr) |>
      evaluate_log_lik() |>
      expect_equal(latent_infec_obs_log_lik(house, tiny, discr))

    # test complex print, discrete and continuous time:
    # the discrete one has the community rate defined, the others stay default
    if (discr) {
      strata <- data.frame(age = c(NA, "child"),
                           param_name = c(".reference", "com_child"),
                           const = c(NA, 0))
      hh_model <- house |>
        transform(age = c("adult", "child", rep("adult", 5))) |>
        build_onset_model(mk_house$param, discr) |>
        add_group_rate("community", strata = strata)
    }
    print(hh_model) |> expect_snapshot()
  }
})

test_that("A case with no source of infection has zero likelihood", {
  # `m2` is the only child, out of reach of the community,
  # and no household transmission
  hh_data <- transform(mk_house$data[c("id", "hh_id", "status", "infec_time",
                                       "end_date")],
                       age = c("adult", "child", rep("adult", 5)))
  strata <- data.frame(age = c(NA, "child"),
                       param_name = c(".reference", "com_child"),
                       const = c(NA, 0))
  param <- modifyList(mk_house$param, list(beta = 0))

  for (discr in c(FALSE, TRUE)) {
    build_onset_model(hh_data, param, discr) |>
      add_group_rate("community", strata = strata) |>
      evaluate_log_lik() |>
      expect_identical(-Inf)
  }
})

test_that("Model quantities are the rates of the force of infection", {
  # `m2` is the only child: its own community, susceptibility, infectivity
  # and contact with the adults
  child <- c(FALSE, TRUE, rep(FALSE, 5))
  house <- transform(mk_house$data[c("id", "hh_id", "status", "infec_time",
                                     "end_date")],
                     age = ifelse(child, "child", "adult"))
  strata <- function(param_name, value) {
    data.frame(age = c(NA, "child"), param_name = c(".reference", param_name),
               const = c(NA, value))
  }
  hh_model <- build_onset_model(house, mk_house$param, FALSE) |>
    add_group_rate(community = list(strata = strata("com_child", 2)),
                   susceptibility = list(strata = strata("sus_child", 3)),
                   infectivity = list(strata = strata("inf_child", 0.5))) |>
    add_contact_rate(strata = data.frame(
      from_age = c("child", NA), to_age = c("adult", NA),
      param_name = c("child_to_adult", ".reference"),
      const = c(1.5, NA)
    ))
  evaluate_log_lik(hh_model)

  # infector by row, infectee by column, no self-contact
  par <- mk_house$param
  beta_star <- par$beta * (par$hh_size_ref / nrow(house))^par$gamma
  ctr <- ifelse(outer(child, child, xor), 1.5, 1)
  diag(ctr) <- 0
  inspect_state(hh_model, "h1")$model_quantities |>
    expect_equal(list(
      community_rate = setNames(par$alpha * ifelse(child, 2, 1), house$id),
      p2p_rate = matrix(beta_star * ifelse(child, 0.5, 1) * ctr *
                          rep(ifelse(child, 3, 1), each = nrow(house)),
                        nrow(house), dimnames = list(house$id, house$id))
    ))
})

test_that("Latent moves preserve log-likelihood", {
  house <- mk_house$data[c("id", "hh_id", "status", "onset_date", "end_date")]
  hh_model <- house |>
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
  hh_fit <- hh_model |>
    run_inference(steps = n_steps, thin = 1L, adapt = 0L, chains = 1L,
                  threads = 2L, id_draws = house$id[house$status > 0]) |>
    suppressWarnings()

  # the pointwise log-lik folds transmission, the latent model and the
  # observation node: the complete-data density, parameter prior aside
  # the bare infection time adds no density: transmission gives it
  log_lik <- unname(get_log_lik_draws(hh_fit)$h1)
  param_draws <- get_parameter_draws(hh_fit) |>
    draws_to_frame()
  latent_draws <- get_latent_draws(hh_fit) |>
    draws_to_frame()
  consts <- get_parameters(hh_model) |>
    subset(!is.na(const), c("name", "const")) |>
    tidyr::pivot_wider(names_from = "name", values_from = "const") |>
    as.list()
  state_at <- function(i) {
    dplyr::left_join(house, get_latent_state(latent_draws, i), by = "id")
  }

  # the snapshot holds the last state
  snapshot <- inspect_state(hh_model, "h1")
  unlist(c(param_draws[n_steps, ], consts))[names(snapshot$parameters)] |>
    expect_identical(snapshot$parameters)
  get_latent_state(latent_draws, n_steps) |>
    expect_identical(snapshot$latent[house$status > 0, ])

  # every retained draw matches an independent R recomputation
  expected <- vapply(seq_len(n_steps), function(i) {
    hh_data <- state_at(i)
    par <- c(param_draws[i, , drop = FALSE], consts)
    infected <- hh_data[hh_data$status > 0, ]
    latent_infec_obs_log_lik(hh_data, par) +
      sum(dnorm(infected$onset_date, infected$incub_period, log = TRUE)) +
      sum(dweibull(infected$incub_period, shape = 2, scale = 1, log = TRUE))
  }, numeric(1))
  log_lik |> expect_equal(expected)

  # the fresh-model likelihood matches the tracked pointwise likelihood
  state_at(n_steps) |>
    transform(infec_draw = infec_time, incub_draw = incub_period) |>
    new_housepi_model("case-ascertained", quiet = TRUE) |>
    add_latent(
      infec_time = list(const = infec_draw),
      incub_period = list(const = incub_draw,
                          model = "weibull", shape = 2, scale = 1)
    ) |>
    add_observation("onset_date", model = "normal",
                    mean = incub_period, sd = 1) |>
    add_parameter(import = c(param_draws[n_steps, , drop = FALSE], consts)) |>
    add_infectivity_profile("gamma", shape = shape, scale = scale) |>
    evaluate_log_lik() |>
    expect_equal(log_lik[n_steps])
})
