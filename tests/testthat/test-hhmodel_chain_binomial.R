library(testthat)

# definitions: `_1` : single infected group
#              `_3` : 3 infected groups
# group rate definitions
groups <- list(community = list("class", reference = "a"),
               susceptibility = list("class", reference = "a"),
               infectivity = list("class", reference = "a"))
# fixed group rates
# `a` is the reference class
grates_1 <- list(a = list(com = 1., sus = 1., inf = 1.),
                 b = list(com = 2., sus = 3., inf = 0.2))
grates_3 <- c(grates_1,
              list(c = list(com = 0.5, sus = 0.8, inf = 1.5)))

# members sorted infected, susceptible
house_3 <- data.frame(
  id = sprintf("hh_01_%s", seq(1, 4)),
  house_id = "hh_01",
  status = c(rep(1, 3), 0),
  age = rep(c("adult", "child"), 2),
  active = c(TRUE, TRUE, FALSE, FALSE),
  class = c("a", "b", rep("c", 2)),
  community = c(1, grates_3$b$com, rep(grates_3$c$com, 2)),
  susceptibility = c(1, grates_3$b$sus, rep(grates_3$c$sus, 2)),
  infectivity = c(1, grates_3$b$inf, rep(grates_3$c$inf, 2))
)
house_1 <- house_3[1:2, ]
house_1$status <- c(1, 0)

# contact rate definitions
contacts <- list(by = "class", symmetric = FALSE,
                 reference = c(from_class = "a", to_class = "b"))
# fixed contact rates
crates_1 <- list(class_b_to_class_a = 2.3)
crates_3 <- list(class_a_to_class_c = 0.7,
                 class_c_to_class_a = 0.6,
                 class_b_to_class_c = 1.1,
                 class_c_to_class_b = 0.9,
                 class_c_to_class_c = 1.2)

# parameters matching automatic names
par_0 <- list(alpha = 0.15, beta = 0.22, gamma = 1, size_ref = 4)
par_1 <- c(par_0, crates_1,
           list(com_class_b = grates_1$b$com,
                sus_class_b = grates_1$b$sus,
                inf_class_b = grates_1$b$inf))
par_3 <- c(par_1, crates_3,
           list(com_class_c = grates_3$c$com,
                sus_class_c = grates_3$c$sus,
                inf_class_c = grates_3$c$inf))
crates_3 <- c(crates_1, crates_3)

# helper building the symmetric contact rate matrix for a given house
make_ctr_mtx <- function(house, crates, prefix = "") {
  n <- nrow(house)
  output <- matrix(0., nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      pair <- c(house$class[i], house$class[j])
      key <- paste0(prefix, "class_", pair[1], "_to_class_", pair[2])
      if (key %in% names(crates)) output[i, j] <- crates[[key]]
      # NOTE: it must reference otherwise, errors are dealt elsewhere
      else output[i, j] <- 1
    }
  }
  output
}

##############################################
## CHAIN BINOMIAL HH MODEL : STANDARD CALLS ##
##############################################

test_that("Chain binomial gaurds against missusage", {
  # not discretizable
  blank_house |>
    new_housepi_model("chain_binomial", discr = TRUE) |>
    expect_error(regexp = "No visible discrete time option")

  # no support for missing status
  data.frame(id = c("m1", "m2"), house_id = "h1",
             status = c(1, NA)) |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template() |>
    run_inference() |>
    expect_error("incompatible")
})

test_that("Deep clone works as expected", {
  house <- data.frame(
    house_id = c(1, 1, 2, 2, 2),
    status = c(1, 0, 1, 0, 0),
    member = c("student", "other", "student", "mother", "father")
  )
  hhmodel <- house |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  out <- suppressWarnings(
    hhmodel$run_inference(steps = 2, chains = 1, adapt = 0, thin = 1)
  )

  # deep clone
  hhclone <- hhmodel$copy()

  # data
  expect_equal(hhmodel$data_wide, hhclone$data_wide)
  expect_equal(hhmodel$data_long, hhclone$data_long)

  # clone setup starts from default configuration
  expect_equal(hhclone$setup, new_housepi_model(house, "chain_binomial")$setup)

  # definitions
  expect_equal(hhmodel$parameters, hhclone$parameters)
  expect_equal(hhmodel$augmented, hhclone$augmented)
  expect_equal(hhmodel$group_rates, hhclone$group_rates)

  # check no backward changes
  hhmodel$add_group_rate("susceptibility", by = "member")
  hhclone <- house |>
    dplyr::mutate(status = NA) |>
    new_housepi_model("chain_binomial")
  expect_failure(expect_equal(hhmodel$group_rates, hhclone$group_rates))
  expect_failure(expect_equal(hhmodel$data_wide, hhclone$data_wide))
})

test_that("Pointwise log-lik is carried into the complete-draws bundle", {
  for (t in c(1L, 2L)) {
    suppressWarnings(
      results <- flu_1918_frederick |>
        new_housepi_model("chain_binomial") |>
        add_parameter_template() |>
        run_inference(steps = 10, thin = 2, adapt = 2,
                      burn = 2, chains = 1, threads = t)
    )

    pointwise <- results |>
      get_log_lik_draws() |>
      posterior::as_draws_matrix()
    # complete_draws renames the pointwise columns to the loo / priorsense
    # convention; their values must be the same pointwise log-likelihood
    bundle <- results |>
      get_complete_draws() |>
      posterior::as_draws_matrix()
    ll_cols <- grep("^log_lik\\[", colnames(bundle))
    expect_equal(length(ll_cols), ncol(pointwise))
    expect_equal(unname(bundle[, ll_cols]), unname(pointwise))
  }
})

test_that("Simulation with observation process is correct", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  rand_observ <- function(x) dplyr::mutate(dplyr::rowwise(x), flag = runif(1))
  withr::with_seed(
    1, result1 <- hhmodel$run_simulation(obs_process = rand_observ)
  )
  # ensure simulation is correct
  expect_true("flag" %in% colnames(result1))
  withr::with_seed(
    1, result1encore <- hhmodel$run_simulation(obs_process = rand_observ)
  )
  expect_identical(result1, result1encore)
})

test_that("Run simulation with nrep > 1 returns list of dataframes", {
  nrep <- 3
  result <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template() |>
    run_simulation(nrep = nrep)
  expect_true(is.list(result))
  expect_length(result, nrep)
  expect_true(is.data.frame(result[[1]]))
})

test_that("run_inference rejects invalid setup arguments", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  hhmodel$run_inference(steps = 5, burn = 10) |>
    expect_error(regex = "[Bb]urn")
  hhmodel$run_inference(steps = 5, thin = 10) |>
    expect_error(regex = "[Tt]hin")
  # threads-clamp warning and short-run diagnostic both fire here
  expect_warning(
    expect_warning(hhmodel$run_inference(steps = 5, chains = 1, adapt = 0,
                                         thin = 1, threads = 1000),
                   regex = "[Tt]hreads"),
    regex = "sample size was below"
  )
})

test_that("run_inference progress bar shows as expected", {
  # FIXME: timing-sensitive
  hhmodel <- flu_1918_baltimora |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  out <- capture.output(
    res <- suppressWarnings(
      run_inference(hhmodel, steps = 2000, chains = 1,
                    adapt = 0, thin = 1, progress = TRUE)
    ),
    type = "output"
  )
  expect_s3_class(res, "HousepiPosteriorDraws")
  expect_true(any(grepl("%|sample", out)),
              info = paste(out, collapse = " | "))
})

test_that("predicted_cases_by_house_size exercises parallel for_each", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  results <- suppressWarnings(
    run_inference(hhmodel, steps = 50, chains = 1, adapt = 0, thin = 1)
  )
  outs <- list()
  for (t in c(1L, 2L)) {
    suppressWarnings(
      outs[[as.character(t)]] <- predicted_cases_by_house_size(
        hhmodel, results, max_nrep = 4, threads = t
      )
    )
    expect_s3_class(outs[[as.character(t)]], "PredictedCasesBySize")
  }
  expect_identical(outs[["1"]]$cases_table_by_size$observed,
                   outs[["2"]]$cases_table_by_size$observed)
})

test_that("Run simulation thins param_table down to nrep when smaller", {
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  pt <- as.data.frame(stats::setNames(
    as.list(hhmodel$parameters$value), hhmodel$parameters$name
  ))
  pt <- pt[rep(1L, 10), , drop = FALSE]
  rownames(pt) <- NULL
  result <- hhmodel$run_simulation(param_table = pt, nrep = 4)
  expect_length(result, 4)
})

####################################################
## CHAIN BINOMIAL HH MODEL : SIMPLE CASE LOG-LIKE ##
####################################################

test_that("Debugger cannot be called on large numbers of infections", {
  n_infec <- 10L
  B <- runif(n_infec)
  Q <- matrix(runif(n_infec * n_infec), ncol = n_infec)
  expect_error(housepi:::eval_chain_binomial_debugger_(B, Q),
               regexp = "Maximum number")
})

test_that("Single infected member log-likelihood is correct", {
  par <- par_1
  house <- house_1
  n_full <- nrow(house)
  n_infec <- sum(house$status > 0)
  # fill house escape rate matrices
  commu_infec_rate <- par$alpha * house$community
  p2p_infec_rate <- par$beta * (par$size_ref / n_full) ^ par$gamma *
    house$infectivity %o% house$susceptibility *
    make_ctr_mtx(house, crates_1)

  # escape probabilities restricted to seropositives
  B <- exp(-commu_infec_rate[seq_len(n_infec)])
  Q <- exp(-p2p_infec_rate[seq_len(n_infec), seq_len(n_infec), drop = FALSE])
  # negative (seronegative) contribution
  neg_log <- -sum(commu_infec_rate[(n_infec + 1):n_full]) -
    sum(p2p_infec_rate[seq_len(n_infec), (n_infec + 1):n_full])
  # positive (seropositive) contribution from the debugger
  pos <- housepi:::eval_chain_binomial_debugger_(B, Q)
  log_lik <- neg_log + log(pos)

  # define housepi model
  hhmodel <- house |>
    new_housepi_model("chain_binomial") |>
    add_group_rate_list(!!!groups)
  do.call(hhmodel$add_contact_rate, contacts)
  hhmodel$add_parameter_list(!!!par, append = FALSE)
  expect_equal(log_lik, eval_log_density(hhmodel))
})

#####################################################
## CHAIN BINOMIAL HH MODEL : COMPLEX CASE LOG-LIKE ##
#####################################################

test_that("Three infected members likelihood is correct", {
  par <- par_3
  house <- house_3
  n_full <- nrow(house)
  # fill escape rate matrices
  commu_infec_rate <- par$alpha * house$community
  p2p_infec_rate <- par$beta * (par$size_ref / n_full) ^ par$gamma *
    house$infectivity %o% house$susceptibility *
    make_ctr_mtx(house, crates_3)

  # CASE 1: no infected
  hhmodel <- house |>
    dplyr::mutate(status = 0) |>
    new_housepi_model("chain_binomial") |>
    add_group_rate_list(!!!groups)
  do.call(hhmodel$add_contact_rate, contacts)
  hhmodel$add_parameter_list(!!!par, append = FALSE)
  log_lik <- -sum(commu_infec_rate)
  expect_equal(eval_log_density(hhmodel), log_lik)

  # CASE 2: mixed infected and susceptible
  n_infec <- sum(house$status > 0)
  hhmodel <- house |>
    new_housepi_model("chain_binomial") |>
    add_group_rate_list(!!!groups)
  do.call(hhmodel$add_contact_rate, contacts)
  hhmodel$add_parameter_list(!!!par, append = FALSE)

  # escape probabilities restricted to seropositives
  B <- exp(-commu_infec_rate[seq_len(n_infec)])
  Q <- exp(-p2p_infec_rate[seq_len(n_infec), seq_len(n_infec), drop = FALSE])
  # negative (seronegative) contribution
  neg_log <- -sum(commu_infec_rate[(n_infec + 1):n_full]) -
    sum(p2p_infec_rate[seq_len(n_infec), (n_infec + 1):n_full])
  # positive (seropositive) contribution from the debugger
  pos <- housepi:::eval_chain_binomial_debugger_(B, Q)
  log_lik <- neg_log + log(pos)

  expect_equal(eval_log_density(hhmodel), log_lik)
})

##################################################
## CHAIN BINOMIAL HH MODEL : REAL CASE LOG-LIKE ##
##################################################

test_that("Chain binomial runs ok on real case with partitions", {
  hhmodel <- flu_2014_matsumoto |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template() |>
    add_group_rate("community", by = "member", reference = ".none") |>
    add_contact_rate(
      strata = data.frame(
        from_member = c("student", "student", "sibling",
                        "mother", "mother",
                        "father",
                        "other",
                        NA),
        to_member = c("student", "sibling", "sibling",
                      "student", "sibling",
                      "mother",
                      "other",
                      NA),
        param_name = c("ctr_child_child", "ctr_child_child", "ctr_child_child",
                       "ctr_mother_child", "ctr_mother_child",
                       ".reference",
                       "ctr_other_other",
                       "ctr_cross_generational")
      )
    )
  expect_silent(hhmodel$evaluate_log_lik())
})

test_that("Chain binomial returns infinite if extreme values are provided", {
  extreme_params <- list(alpha = 1e-20, beta = 1e-20,
                         gamma = 1, size_ref = 4)
  # heterogeneous_house
  hhmodel <- data.frame(
    id = c(1, 2), house_id = "hh1", status = c(1, 1),
    age = c("child", "adult")
  ) |>
    new_housepi_model("chain_binomial") |>
    add_parameter_list(!!!extreme_params)

  # homogeneous branch
  expect_false(is.finite(eval_log_density(hhmodel)))

  # heterogeneous branch
  expect_false(is.finite(
    hhmodel |>
      add_group_rate("susceptibility", by = "age") |>
      eval_log_density()
  ))
})
