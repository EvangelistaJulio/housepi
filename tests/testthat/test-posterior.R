library(testthat)

# Run inference once and reuse across tests.
# Short run trips the diagnostic short-run warning — irrelevant here.
local_results <- local({
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter_template()
  suppressWarnings(run_inference(hhmodel,
    steps = 100, chains = 1, adapt = 0, thin = 1, progress = FALSE
  ))
})

test_that("Active fields are read-only", {
  local_results |>
    helper_test_read_only(
      "parameter_draws", "augmented_draws",
      "log_lik_draws", "lprior_draws", "complete_draws"
    )
})

test_that("Invalid Credibile Interval level is rejected", {
  set_global_cri_level(0.5) |>
    expect_error("between 0.5 and 1")
  set_global_cri_level(1) |>
    expect_error("between 0.5 and 1")
  set_global_cri_level(NA_real_) |>
    expect_error("non-NA")
  set_global_cri_level("0.9") |>
    expect_error("numeric")
  set_global_cri_level(c(0.9, 0.95)) |>
    expect_error("single")
})

test_that("Get-global-cri-level behaves correctly", {
  withr::with_options(list(housepi.cri_level = NULL), {
    expect_equal(get_global_cri_level(), 0.95)
  })
  withr::with_options(list(housepi.cri_level = 0.95), {
    old <- set_global_cri_level(0.90)
    expect_equal(old, 0.95)
    expect_equal(get_global_cri_level(), 0.90)
  })
})

test_that("Summarise-draws default columns track cri_level", {
  skip_if_not_installed("posterior")
  plain <- posterior::example_draws()
  d <- plain
  class(d) <- unique(c("housepi_draws", class(d)))

  # spelling alias dispatches to the same handler
  expect_identical(posterior::summarize_draws(d), posterior::summarise_draws(d))

  # housepi method matches posterior::summarise_draws on every column
  # except the credible-interval quantile pair, which tracks cri_level.
  ref_cols <- setdiff(colnames(posterior::summarise_draws(plain)),
                      c("q5", "q95"))
  withr::with_options(list(housepi.cri_level = 0.95), {
    s <- posterior::summarise_draws(d)
    expect_true(all(ref_cols %in% colnames(s)))
    expect_true(all(c("q2.5", "q97.5") %in% colnames(s)))
  })
  withr::with_options(list(housepi.cri_level = 0.80), {
    s <- posterior::summarise_draws(d)
    expect_true(all(ref_cols %in% colnames(s)))
    expect_true(all(c("q10", "q90") %in% colnames(s)))
  })

  # custom summary functions fall through to vanilla posterior behaviour
  s <- posterior::summarise_draws(d, mean, sd)
  expect_setequal(colnames(s), c("variable", "mean", "sd"))

  # plain draws remain untouched by housepi (vanilla q5 / q95 quantiles)
  plain_cols <- colnames(posterior::summarise_draws(plain))
  expect_true(all(c("q5", "q95") %in% plain_cols))
  expect_false(any(c("q2.5", "q97.5", "q10", "q90") %in% plain_cols))
})

test_that("Draws active fields are correct", {
  pd <- local_results$parameter_draws
  expect_true(posterior::is_draws(pd))
  expect_equal(posterior::niterations(pd), 100L)

  # empty draws
  ad <- local_results$augmented_draws
  expect_true(posterior::is_draws(ad))
  expect_equal(posterior::nvariables(ad), 0)

  ll <- local_results$log_lik_draws
  expect_true(posterior::is_draws(ll))
  expect_true(posterior::nvariables(ll) > 0)

  lp <- local_results$lprior_draws
  expect_true(posterior::is_draws(lp))
  expect_true("lprior" %in% posterior::variables(lp))

  cd <- local_results$complete_draws
  expect_true(posterior::is_draws(cd))
  vars <- posterior::variables(cd)
  # should contain param names, log_lik[...], and lprior (loo / priorsense names)
  expect_true(any(grepl("log_lik\\[", vars)))
  expect_true("lprior" %in% vars)
})

test_that("Exported getter functions work correctly", {
  expect_true(posterior::is_draws(get_parameter_draws(local_results)))
  expect_true(posterior::is_draws(get_log_lik_draws(local_results)))
  expect_true(posterior::is_draws(get_complete_draws(local_results)))

  expect_error(get_parameter_draws("snake.eater"))
  expect_error(get_augmented_draws(42))
  expect_error(get_log_lik_draws(list()))
  expect_error(get_lprior_draws(NULL))
  expect_error(get_complete_draws(lnorm))
})

test_that("Empty draws getters warn", {
  # local_results has parameters but no augmented draws
  get_parameter_draws(local_results) |>
    expect_silent() |>
    posterior::is_draws() |>
    expect_true()

  get_augmented_draws(local_results) |>
    expect_warning("id_draws")
})

test_that("Print behavior is stable", {
  expect_snapshot(print(local_results))
})

test_that("Multi-chain inference produces correct chain structure", {
  hhmodel <- new_housepi_model(flu_1918_frederick, "final_size")
  hhmodel <- add_parameter_template(hhmodel)
  results <- suppressWarnings(
    run_inference(hhmodel, steps = 10, chains = 2, adapt = 0, thin = 1,
                  progress = FALSE)
  )

  pd <- results$parameter_draws
  expect_true(posterior::is_draws(pd))
  expect_equal(posterior::nchains(pd), 2L)
  expect_equal(posterior::niterations(pd), 10L)
})

test_that("Augm draws are not empty when id_draws is specified", {
  hhmodel <- data.frame(house_id = 1,
                        id = c("a", "b"),
                        status = c(1, 0)) |>
    new_housepi_model("test")
  hhmodel <- add_augmented(hhmodel, name = "proba", prior = "uniform")
  results <- suppressWarnings(
    run_inference(hhmodel, steps = 5, chains = 1, adapt = 0, thin = 1,
                  id_draws = "a", progress = FALSE)
  )
  ad <- results$augmented_draws
  expect_true(posterior::is_draws(ad))
  expect_true(posterior::nvariables(ad) > 0)
})

test_that("Diagnostic warning fires on short runs", {
  flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter_template() |>
    run_inference(steps = 100, chains = 1, adapt = 0, thin = 1,
                  progress = FALSE) |>
    expect_warning("sample size was below")
})

test_that("Diagnostic warning fires on poor mixing", {
  # explicitly define micro-steps
  expect_warning(
    results <- flu_1918_frederick |>
      new_housepi_model("final_size") |>
      add_parameter(
        alpha = list(prior = "uniform", proposal = "rw_log",
                     proposal_scale = 1e-4),
        beta  = list(prior = "uniform", proposal = "rw_log",
                     proposal_scale = 1e-4),
        gamma = 1, hhsize_ref = 5
      ) |>
      run_inference(steps = 1000, thin = 1, adapt = 0, chains = 2,
                    progress = FALSE)
  )
})
