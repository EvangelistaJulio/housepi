library(testthat)

test_that("Deep clone works as expected", {
  hhmodel <- new_housepi_model("chain_binomial_heterogeneous")
  hhmodel$data_wide <- housepi::survey
  hhmodel$load_parameter_template()
  expect_no_error(out <- hhmodel$infer(steps = 2, progress = FALSE))

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
    infer(steps = 10, thin = 1, adapt = 0, burnin = 0)

  pointwise <- results$log_lik_draws |>
    posterior::as_draws_matrix() |>
    rowSums() |>
    unname()
  expect_equal(results$total_log_lik_draws$log_lik, pointwise)
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
