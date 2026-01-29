library(testthat)

test_that("Chain binomial is not discretizable", {
  expect_warning(
    new_housepi_model("chain_binomial", discr = TRUE),
    regexp = "No visible discrete time option"
  )
  expect_warning(
    new_housepi_model("chain_binomial_heterogeneous", discr = TRUE),
    regexp = "No visible discrete time option"
  )
})

test_that("Deep clone works as expected", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$data_wide <- housepi::survey
  hhmodel$load_parameter_template()
  expect_no_error(out <- hhmodel$run_inference(steps = 2, progress = FALSE))

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
  hhclone$data_wide <- dplyr::mutate(housepi::survey, status = NA)
  expect_failure(expect_equal(hhmodel$group_rates, hhclone$group_rates))
  expect_failure(expect_equal(hhmodel$data_wide, hhclone$data_wide))
})

test_that("Log-likelihood traces and pointwise traces are coherent", {
  results <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    add_data_wide(housepi::survey) |>
    run_inference(steps = 10, thin = 1, adapt = 0, burnin = 0)

  pointwise <- results$log_lik_draws |>
    posterior::as_draws_matrix() |>
    rowSums() |>
    unname()
  expect_equal(results$total_log_lik_draws$log_lik, pointwise)
})

test_that("Single infected member log-likelihood is correct", {
  # fixed group rates
  com_class_B <- 1.41
  sus_class_B <- 1.03
  inf_class_B <- 1.3
  # members sorted infected, susceptible
  # 'A' is the reference class
  house <- data.frame(id = sprintf("hh_01_%s", seq(1, 2)),
                      house_id = "hh_01",
                      status = c(1, 0),
                      class = c("A", "B"),
                      community = c(1, com_class_B),
                      susceptibility = c(1, sus_class_B),
                      infectivity = c(1, 1.3))

  # fill house escape rate matrices
  par <- list(alpha = 0.15, beta = 0.22)
  commu_esc_rate <- par$alpha * house$community
  p2p_esc_rate <- par$beta * house$infectivity %o% house$susceptibility

  # susceptible member
  log_lik <- - commu_esc_rate[2] - p2p_esc_rate[1, 2] + 
    log1p(-exp(-commu_esc_rate[1]))

  # define housepi model
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$data_wide <- house
  hhmodel$load_group_rate_list(community = "class",
                               susceptibility = "class",
                               infectivity = "class")
  hhmodel$parameters <- list(alpha = par$alpha, beta = par$beta,
                             com_class_B = com_class_B,
                             sus_class_B = sus_class_B,
                             inf_class_B = inf_class_B)
  expect_equal(log_lik, hhmodel$evaluate_log_lik())
})

test_that("Chain binomial models are coherent", {
  # log-likelihood from chain binomial model
  hhmodel1 <- new_housepi_model("chain_binomial")

  hhmodel1$parameters <- list(alpha = 0.15, beta = 0.22)
  hhmodel1$data_wide <- housepi::survey
  ll1 <- hhmodel1$evaluate_log_lik()

  # log-likelihood from heterogeneous chain binomial model
  # with homogenous rates
  hhmodel2 <- new_housepi_model("chain_binomial_heterogeneous")

  hhmodel2$parameters <- list(alpha = 0.15, beta = 0.22)
  hhmodel2$data_wide <- housepi::survey
  ll2 <- hhmodel2$evaluate_log_lik()

  expect_equal(ll1, ll2)
})

# UPDATE: TO BE CHECKED...
# test_that("Three infected members likelihood is correct", {
#   # fixed group rates
#   com_class_B <- 1.41
#   sus_class_B <- 1.03
#   inf_class_B <- 1.3
#   com_class_C <- 0.73
#   sus_class_C <- 0.89
#   inf_class_C <- 1.51
#   # members sorted infected, susceptible
#   # 'A' is the reference class
#   house <- data.frame(
#     id = sprintf("hh_01_%s", seq(1, 4)),
#     house_id = "hh_01",
#     status = c(1, 1, 1, 0),
#     class = c("A", "B", "C", "C"),
#     community = c(1, com_class_B, com_class_C, com_class_C),
#     susceptibility = c(1, sus_class_B, sus_class_C, sus_class_C),
#     infectivity = c(1, inf_class_B, inf_class_C, inf_class_C)
#   )
# 
#   # fill house escape rate matrices
#   par <- list(alpha = 0.15, beta = 0.22)
#   commu_esc_rate <- par$alpha * house$community
#   p2p_esc_rate <- par$beta * house$infectivity %o% house$susceptibility
# 
#   # ... escape probabilities of infected
#   commu_esc_proba <- exp(-commu_esc_rate[1:3])
#   p2p_esc_proba <- exp(-p2p_esc_rate[1:3,1:3])
#   diag(p2p_esc_proba) <- 1
#   # ... infection probabilities
#   commu_inf_proba <- 1 - commu_esc_proba
#   p2p_inf_proba <- 1 - p2p_esc_proba
# 
#   # susceptible member
#   log_lik <- - (commu_esc_rate[4] +
#     p2p_esc_rate[1, 4] + p2p_esc_rate[2, 4] + p2p_esc_rate[3, 4])
# 
#   # order 0: 0->1, 0->2, 0->3
#   lik_0 <- prod(commu_inf_proba)
# 
#   # order 1: 0->1, 0->2, (1->3 + 2->3)
#   lik_1 <- function(a, b, c) {
#     commu_inf_proba[a] * commu_inf_proba[b] * commu_esc_proba[c] *
#     (p2p_inf_proba[a, c] * p2p_esc_proba[b, c] +
#        p2p_inf_proba[b, c] * p2p_esc_proba[a, c])
#   }
# 
#   # order 2: 0->1, 1->2, 1->3
#   lik_2 <- function(a, b, c) {
#     commu_inf_proba[a] * commu_esc_proba[b] * commu_esc_proba[c] *
#     p2p_inf_proba[a, b] * p2p_inf_proba[c, b]
#   }
# 
#   # order 3: 0->1, ((1->2, 2->3) + (1->3, 3->2))
#   lik_3 <- function(a, b, c) {
#     commu_inf_proba[a] * commu_esc_proba[b] * commu_esc_proba[c] *
#     (p2p_inf_proba[a, b] * p2p_esc_proba[a, c] * p2p_inf_proba[b, c] +
#        p2p_inf_proba[a, c] * p2p_esc_proba[a, b] * p2p_inf_proba[c, b])
#   }
# 
#   log_lik <- log_lik +
#     log(lik_0 + lik_1(1, 2, 3) + lik_1(2, 3, 1) + lik_1(3, 1, 2) +
#           lik_2(1, 2, 3) + lik_2(2, 3, 1) + lik_2(3, 1, 2) +
#           lik_3(1, 2, 3) + lik_3(2, 3, 1) + lik_3(3, 1, 2))
# 
#   # define housepi model
#   hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
#   hhmodel$data_wide <- house
#   hhmodel$load_group_rate_list(community = list("class", reference = "A"),
#                                susceptibility = list("class", reference = "A"),
#                                infectivity = list("class", reference = "A"))
#   hhmodel$parameters <- list(alpha = par$alpha, beta = par$beta,
#                              com_class_B = com_class_B,
#                              sus_class_B = sus_class_B,
#                              inf_class_B = inf_class_B,
#                              com_class_C = com_class_C,
#                              sus_class_C = sus_class_C,
#                              inf_class_C = inf_class_C)
#   expect_equal(log_lik, hhmodel$evaluate_log_lik())
# })
