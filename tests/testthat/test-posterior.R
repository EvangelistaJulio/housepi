library(testthat)

# Run inference once and reuse across tests
local_results <- local({
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(housepi::survey)
  hhmodel$load_parameter_template()
  hhmodel$run_inference(steps = 20, burnin = 0, thin = 1,
                        chains = 1, progress = FALSE)
})

############################################################################
## Active fields
############################################################################

test_that("param_draws returns posterior::draws_df", {
  pd <- local_results$param_draws
  expect_true(posterior::is_draws(pd))
  expect_true(posterior::niterations(pd) == 20)
})

test_that("prior_draws returns draws_df with lprior column", {
  pr <- local_results$prior_draws
  expect_true(posterior::is_draws(pr))
  expect_true("lprior" %in% posterior::variables(pr))
})

test_that("augm_draws returns empty draws for non-augmented model", {
  ad <- local_results$augm_draws
  expect_true(posterior::is_draws(ad))
  expect_equal(posterior::nvariables(ad), 0)
})

test_that("log_lik_draws returns posterior::draws_df", {
  ll <- local_results$log_lik_draws
  expect_true(posterior::is_draws(ll))
  expect_true(posterior::nvariables(ll) > 0)
})

test_that("total_log_lik_draws returns posterior::draws_df", {
  tll <- local_results$total_log_lik_draws
  expect_true(posterior::is_draws(tll))
  expect_true("log_lik" %in% posterior::variables(tll))
})

test_that("complete_draws binds param, log_lik, and prior draws", {
  cd <- local_results$complete_draws
  expect_true(posterior::is_draws(cd))
  vars <- posterior::variables(cd)
  # should contain param names, log_lik[...], and lprior
  expect_true(any(grepl("log_lik\\[", vars)))
  expect_true("lprior" %in% vars)
})

############################################################################
## Read-only enforcement
############################################################################

test_that("active fields are read-only", {
  regexp <- "read only"
  expect_error(local_results$param_draws <- "foo", regexp = regexp)
  expect_error(local_results$prior_draws <- "foo", regexp = regexp)
  expect_error(local_results$augm_draws <- "foo", regexp = regexp)
  expect_error(local_results$log_lik_draws <- "foo", regexp = regexp)
  expect_error(local_results$total_log_lik_draws <- "foo", regexp = regexp)
  expect_error(local_results$complete_draws <- "foo", regexp = regexp)
})

############################################################################
## Exported getters
############################################################################

test_that("exported getter functions work", {
  expect_true(posterior::is_draws(get_param_draws(local_results)))
  expect_true(posterior::is_draws(get_prior_draws(local_results)))
  expect_true(posterior::is_draws(get_augm_draws(local_results)))
  expect_true(posterior::is_draws(get_log_lik_draws(local_results)))
  expect_true(posterior::is_draws(get_total_log_lik_draws(local_results)))
  expect_true(posterior::is_draws(get_complete_draws(local_results)))
})

test_that("exported getters reject wrong types", {
  expect_error(get_param_draws("not_draws"))
  expect_error(get_augm_draws(42))
  expect_error(get_log_lik_draws(list()))
  expect_error(get_total_log_lik_draws(NULL))
  expect_error(get_prior_draws(data.frame()))
  expect_error(get_complete_draws("x"))
})

############################################################################
## print and summarize
############################################################################

test_that("print works without error", {
  expect_snapshot(print(local_results))
})

test_that("summarize_all_draws returns summary table", {
  summary_df <- local_results$summarize_all_draws()
  expect_true(is.data.frame(summary_df))
  expect_true("variable" %in% colnames(summary_df))
  expect_true("mean" %in% colnames(summary_df))
  expect_true(nrow(summary_df) > 0)
})

############################################################################
## Multi-chain inference
############################################################################

test_that("multi-chain inference produces correct chain structure", {
  hhmodel <- new_housepi_model("chain_binomial")
  hhmodel$add_data_wide(housepi::survey)
  hhmodel$load_parameter_template()
  results <- hhmodel$run_inference(steps = 10, chains = 2, progress = FALSE)

  pd <- results$param_draws
  expect_true(posterior::is_draws(pd))
  expect_equal(posterior::nchains(pd), 2)
  expect_equal(posterior::niterations(pd), 10)
})

############################################################################
## Augmented draws (blank model)
############################################################################

test_that("augm_draws non-empty when id_draws specified", {
  hhmodel <- new_housepi_model("augmented_example")
  hhmodel$add_data_wide(
    data.frame(house_id = 1, id = c("a", "b"), status = c(1, 0))
  )
  hhmodel$add_augmented(name = "proba", model = "uniform")
  results <- hhmodel$run_inference(steps = 5, id_draws = "a",
                                   progress = FALSE)
  ad <- results$augm_draws
  expect_true(posterior::is_draws(ad))
  expect_true(posterior::nvariables(ad) > 0)
})
