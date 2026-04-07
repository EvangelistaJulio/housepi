library(testthat)

# definitions: `_1` : single infected group
#              `_3` : 3 infected groups
# group rate definitions
groups <- list(community = list("class", reference = "a"),
               susceptibility = list("class", reference = "a"),
               infectivity = list("class", reference = "a"))
# fixed group rates
# `a` is the reference class
rates_1 <- list(a = list(com = 1., sus = 1., inf = 1.),
                b = list(com = 1.41, sus = 1.03, inf = 1.3))
rates_3 <- c(rates_1,
             list(c = list(com = 0.73, sus = 0.89, inf = 1.51)))
# parameters matching automatic names
par_1 <- list(alpha = 0.15, beta = 0.22,
              com_class_b = rates_1$b$com,
              sus_class_b = rates_1$b$sus,
              inf_class_b = rates_1$b$inf)
par_3 <- c(par_1,
           list(com_class_c = rates_3$c$com,
                sus_class_c = rates_3$c$sus,
                inf_class_c = rates_3$c$inf))
# members sorted infected, susceptible
house_3 <- data.frame(
  id = sprintf("hh_01_%s", seq(1, 4)),
  house_id = "hh_01",
  status = c(rep(1, 3), 0),
  class = c("a", "b", rep("c", 2)),
  community = c(1, rates_3$b$com, rep(rates_3$c$com, 2)),
  susceptibility = c(1, rates_3$b$sus, rep(rates_3$c$sus, 2)),
  infectivity = c(1, rates_3$b$inf, rep(rates_3$c$inf, 2))
)
house_1 <- house_3[1:2, ]
house_1$status <- c(1, 0)

test_that("Chain binomial is not discretizable", {
  expect_error(new_housepi_model("chain_binomial", discr = TRUE),
               regexp = "No visible discrete time option")
  expect_error(new_housepi_model("chain_binomial_heterogeneous", discr = TRUE),
               regexp = "No visible discrete time option")
})

test_that("Deep clone works as expected", {
  hhmodel <- new_housepi_model("chain_binomial_grouped")
  hhmodel$add_data_wide(helper_test_get_data("flu_2014_matsumoto"))
  hhmodel$load_parameter_template()
  expect_no_error(out <- hhmodel$run_inference(steps = 2, progress = FALSE))

  # deep clone
  hhclone <- hhmodel$copy()

  # data
  expect_equal(hhmodel$data_wide, hhclone$data_wide)
  expect_equal(hhmodel$data_long, hhclone$data_long)

  # definitions
  expect_equal(hhmodel$setup, hhclone$setup)
  expect_equal(hhmodel$parameters, hhclone$parameters)
  expect_equal(hhmodel$augmented, hhclone$augmented)
  expect_equal(hhmodel$group_rates, hhclone$group_rates)

  # check no backward changes
  hhmodel$add_group_rate("susceptibility", by = "member")
  hhclone$add_data_wide(
    dplyr::mutate(helper_test_get_data("flu_2014_matsumoto"),
                  status = NA)
  )
  expect_failure(expect_equal(hhmodel$group_rates, hhclone$group_rates))
  expect_failure(expect_equal(hhmodel$data_wide, hhclone$data_wide))
})

test_that("Log-likelihood traces and pointwise traces are coherent", {
  results <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    add_data_wide(helper_test_get_data("flu_1918_frederick")) |>
    run_inference(steps = 10, thin = 2, adapt = 2, burn = 2)

  pointwise <- results$log_lik_draws |>
    posterior::as_draws_matrix() |>
    rowSums() |>
    unname()
  expect_equal(results$total_log_lik_draws$log_lik, pointwise)
})

test_that("Single infected member log-likelihood is correct", {
  par <- par_1
  house <- house_1
  # fill house escape rate matrices
  commu_esc_rate <- par$alpha * house$community
  p2p_esc_rate <- par$beta * house$infectivity %o% house$susceptibility

  # susceptible member
  log_lik <- - commu_esc_rate[2] - p2p_esc_rate[1, 2] +
    log1p(-exp(-commu_esc_rate[1]))

  # define housepi model
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$add_data_wide(house)
  hhmodel$load_group_rate_list(!!!groups)
  hhmodel$load_parameter_list(!!!par)
  expect_equal(log_lik, hhmodel$evaluate_log_lik())
})

test_that("Chain binomial models are coherent", {
  par <- par_1[c("alpha", "beta")]

  # log-likelihood from chain binomial model
  hhmodel_1 <- new_housepi_model("chain_binomial")
  hhmodel_1$load_parameter_list(!!!par)
  hhmodel_1$add_data_wide(helper_test_get_data("sarscov2_2021_icofs"))

  # log-likelihood from heterogeneous chain binomial model
  # with homogenous rates
  hhmodel_2 <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel_2$load_parameter_list(!!!par)
  hhmodel_2$add_data_wide(helper_test_get_data("sarscov2_2021_icofs"))

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Three infected members likelihood is correct", {
  par <- par_3
  house <- house_3
  # fill house escape rate matrices
  commu_esc_rate <- par$alpha * house$community
  p2p_esc_rate <- par$beta * house$infectivity %o% house$susceptibility

  # ... escape probabilities of infected
  commu_esc_proba <- exp(-commu_esc_rate[1:3])
  p2p_esc_proba <- exp(-p2p_esc_rate[1:3, 1:3])
  diag(p2p_esc_proba) <- 1
  # ... infection probabilities
  commu_inf_proba <- 1 - commu_esc_proba

  # susceptible member
  log_lik <- -commu_esc_rate[4] -
    p2p_esc_rate[1, 4] - p2p_esc_rate[2, 4] - p2p_esc_rate[3, 4]

  # Generation-based decomposition of seropositive contribution.
  # gen [3]: all 3 infected by community
  lik_0 <- prod(commu_inf_proba)

  # gen [2, 1]: 2 in gen 0 (a,b from community), 1 in gen 1 (c from p2p)
  lik_1 <- function(a, b, c) {
    commu_inf_proba[a] * commu_inf_proba[b] * commu_esc_proba[c] *
      (1 - p2p_esc_proba[a, c] * p2p_esc_proba[b, c])
  }

  # gen [1, 2]: 1 in gen 0 (a from community), 2 in gen 1 (b,c from a)
  lik_2 <- function(a, b, c) {
    commu_inf_proba[a] * commu_esc_proba[b] * commu_esc_proba[c] *
      (1 - p2p_esc_proba[a, b]) * (1 - p2p_esc_proba[a, c])
  }

  # gen [1, 1, 1]: a from community, b from a (gen 1), c from b (gen 2)
  #   two sub-chains: a->b->c and a->c->b
  lik_3 <- function(a, b, c) {
    commu_inf_proba[a] * commu_esc_proba[b] * commu_esc_proba[c] *
      ((1 - p2p_esc_proba[a, b]) * p2p_esc_proba[a, c] * (1 - p2p_esc_proba[b, c]) +
        (1 - p2p_esc_proba[a, c]) * p2p_esc_proba[a, b] * (1 - p2p_esc_proba[c, b]))
  }

  log_lik <- log_lik +
    log(lik_0 + lik_1(1, 2, 3) + lik_1(2, 3, 1) + lik_1(3, 1, 2) +
          lik_2(1, 2, 3) + lik_2(2, 3, 1) + lik_2(3, 1, 2) +
          lik_3(1, 2, 3) + lik_3(2, 3, 1) + lik_3(3, 1, 2))

  # define housepi heterogeneous model
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$add_data_wide(house)
  hhmodel$load_group_rate_list(!!!groups)
  hhmodel$load_parameter_list(!!!par)
  expect_equal(log_lik, hhmodel$evaluate_log_lik())
})

test_that("Grouped model with single type matches homogeneous", {
  par <- par_1[c("alpha", "beta")]

  # log-likelihood from homogeneous chain binomial model
  hhmodel_1 <- new_housepi_model("chain_binomial")
  hhmodel_1$load_parameter_list(!!!par)
  hhmodel_1$add_data_wide(helper_test_get_data("flu_1918_frederick"))

  # log-likelihood from grouped chain binomial model
  # with homogeneous rates (single type)
  hhmodel_2 <- new_housepi_model("chain_binomial_grouped")
  hhmodel_2$load_parameter_list(!!!par)
  hhmodel_2$add_data_wide(helper_test_get_data("flu_1918_frederick"))

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Grouped model single infected member matches heterogeneous", {
  par <- par_1
  house <- house_1

  # heterogeneous model
  hhmodel_1 <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel_1$add_data_wide(house)
  hhmodel_1$load_group_rate_list(!!!groups)
  hhmodel_1$load_parameter_list(!!!par)

  # grouped model
  hhmodel_2 <- new_housepi_model("chain_binomial_grouped")
  hhmodel_2$add_data_wide(house)
  hhmodel_2$load_group_rate_list(!!!groups)
  hhmodel_2$load_parameter_list(!!!par)

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Grouped model with all unique types matches heterogeneous", {
  par <- par_3[!stringr::str_ends(names(par_3), "_c")]
  # 3 infected of type b; 1 susceptible of type a
  house <- house_3
  house[4, ] <- house_3[1, ]
  house[1, ] <- house_3[2, ]
  house[3, ] <- house_3[2, ]
  house$id <- sprintf("hh_01_%s", seq(1, dim(house)[1]))
  house$status <- c(rep(1, 3), 0)

  # heterogeneous model
  hhmodel_1 <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel_1$add_data_wide(house)
  hhmodel_1$load_group_rate_list(!!!groups)
  hhmodel_1$load_parameter_list(!!!par)

  # grouped model
  hhmodel_2 <- new_housepi_model("chain_binomial_grouped")
  hhmodel_2$add_data_wide(house)
  hhmodel_2$load_group_rate_list(!!!groups)
  hhmodel_2$load_parameter_list(!!!par)

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Grouped model matches heterogeneous for mixed types", {
  par <- par_3
  house <- house_3

  # heterogeneous model
  hhmodel_1 <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel_1$add_data_wide(house)
  hhmodel_1$load_group_rate_list(!!!groups)
  hhmodel_1$load_parameter_list(!!!par)

  # grouped model
  hhmodel_2 <- new_housepi_model("chain_binomial_grouped")
  hhmodel_2$add_data_wide(house)
  hhmodel_2$load_group_rate_list(!!!groups)
  hhmodel_2$load_parameter_list(!!!par)

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Grouped model matches heterogeneous for mixed types (size 6)", {
  par <- par_3
  # 5 infected: 3 of type a, 2 of type b; 1 susceptible of type c
  house <- house_3
  house[1, ] <- house_3[1, ]
  house[2, ] <- house_3[1, ]
  house[3, ] <- house_3[1, ]
  house[4, ] <- house_3[2, ]
  house[5, ] <- house_3[2, ]
  house[6, ] <- house_3[3, ]

  house$id <- sprintf("hh_01_%s", seq(1, dim(house)[1]))
  house$status <- c(rep(1, 5), 0)

  # heterogeneous model
  hhmodel_1 <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel_1$add_data_wide(house)
  hhmodel_1$load_group_rate_list(!!!groups)
  hhmodel_1$load_parameter_list(!!!par)

  # grouped model
  hhmodel_2 <- new_housepi_model("chain_binomial_grouped")
  hhmodel_2$add_data_wide(house)
  hhmodel_2$load_group_rate_list(!!!groups)
  hhmodel_2$load_parameter_list(!!!par)

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Grouped model matches heterogeneous for complex dataset", {
  more_groups <- list(infectivity = "isolation",
                      susceptibility = "vacc_status",
                      community = "sex")

  hhmodel_1 <- new_housepi_model("chain binomial grouped") |>
    load_group_rate_list(!!!more_groups) |>
    add_data_wide(housepi::sarscov2_2021_icofs)

  # UPDATE: this is a workaround to have constant group rates
  param_names <- hhmodel_1$parameters$name
  values <- runif(length(param_names), min = 0.2, max = 1.8)
  params <- c(par_1[c("alpha", "beta")],
              setNames(as.list(values), param_names))
  hhmodel_1$load_parameter_list(!!!params)

  hhmodel_2 <- new_housepi_model("chain binomial heterogeneous") |>
    load_group_rate_list(!!!more_groups) |>
    add_data_wide(housepi::sarscov2_2021_icofs) |>
    load_parameter_list(!!!params)

  expect_equal(hhmodel_1$evaluate_log_lik(), hhmodel_2$evaluate_log_lik())
})

test_that("Simulation with observation process is correct", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$load_parameter_template()
  hhmodel$add_data_wide(helper_test_get_data("flu_1918_frederick"))
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

test_that("Run simulation with nsim > 1 returns list of dataframes", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(helper_test_get_data("flu_1918_frederick"))
  hhmodel$load_parameter_template()
  nsim <- 3
  result <- hhmodel$run_simulation(nsim = nsim)
  expect_true(is.list(result))
  expect_length(result, nsim)
  expect_true(is.data.frame(result[[1]]))
})
