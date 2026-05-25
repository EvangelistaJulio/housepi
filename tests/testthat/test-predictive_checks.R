library(testthat)

# Shared fixture: a freshly templated flu_1918_frederick chain-binomial model
# and one inference run, reused across PPC tests that don't mutate either.
local_ppc <- local({
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  results <- suppressWarnings(run_inference(
    hhmodel, steps = 100, chains = 1, adapt = 0, thin = 1, progress = FALSE
  ))
  list(hhmodel = hhmodel, results = results)
})

test_that("measure_sar_by_house_size computes SAR correctly", {
  data_wide <- data.frame(
    id = paste0("p", 1:7),
    house_id = c(rep("hh1", 2), rep("hh2", 3), rep("hh3", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0)
  )

  result <- measure_sar_by_house_size(data_wide)

  expect_true(is.data.frame(result))
  expect_true(all(c("house_size", "sar") %in% colnames(result)))
  # size 2: mean of SAR 1.0 (hh1) and 0.0 (hh3) = 0.5
  expect_equal(result[result$house_size == 2, "sar"], 0.5)
  # size 3: SAR = (1-1)/(3-1) = 0
  expect_equal(result[result$house_size == 3, "sar"], 0)

})

test_that("measure_cases_by_house_size produces correct contingency table", {
  data_wide <- data.frame(id = paste0("p", 1:5),
                          house_id = c(rep("hh1", 2), rep("hh2", 3)),
                          status = c(1, 1, 1, 0, 0))

  result <- measure_cases_by_house_size(data_wide)

  expect_true(is.data.frame(result))
  expect_true("house_size" %in% colnames(result))
  # grid should be complete
  expect_true(nrow(result) > 0)

})

test_that("validate_cri_level validates credible interval", {
  expect_equal(validate_cri_level(0.95), 0.95)
  expect_equal(validate_cri_level(0.8), 0.8)
  # boundary and invalid values
  expect_error(validate_cri_level(0.5), regexp = "between 0.5 and 1")
  expect_error(validate_cri_level(1.0), regexp = "between 0.5 and 1")
  expect_error(validate_cri_level(0.3), regexp = "between 0.5 and 1")
  expect_error(validate_cri_level(NA_real_), regexp = "non-NA")
})

test_that("ppc_frame_apply computes element-wise operations", {
  # list of data.frames
  df1 <- data.frame(a = c(1, 2), b = c(3, 4))
  df2 <- data.frame(a = c(5, 6), b = c(7, 8))
  result <- ppc_frame_apply(list(df1, df2), mean)
  expect_equal(result$a, c(3, 4))
  expect_equal(result$b, c(5, 6))

  # named vectors
  v1 <- c(x = 1, y = 2)
  v2 <- c(x = 3, y = 4)
  result <- ppc_frame_apply(list(v1, v2), mean)
  expect_equal(result, c(x = 2, y = 3))
})

test_that("ppc_frame_apply_summary handles mean/median, NA and non-numeric", {
  # multi-row data.frames required by internal apply() on matrix path
  df1 <- data.frame(a = c(1, 2), b = c(NA_real_, NA_real_))
  df2 <- data.frame(a = c(3, 4), b = c(NA_real_, NA_real_))

  result <- ppc_frame_apply_summary(list(df1, df2), "mean")
  expect_equal(result$a, c(2, 3))
  expect_true(all(is.na(result$b)))

  # median path
  df3 <- data.frame(a = c(10, 20), b = c(NA_real_, NA_real_))
  result_med <- ppc_frame_apply_summary(list(df1, df2, df3), "median")
  expect_equal(result_med$a, c(3, 4))

  # named vector path: non-numeric uniform values
  v1 <- c(x = 1, y = NA_real_)
  v2 <- c(x = 3, y = NA_real_)
  result <- ppc_frame_apply_summary(list(v1, v2), "mean")
  expect_equal(result[["x"]], 2)
  expect_true(is.na(result[["y"]]))

  # invalid summary stat
  expect_error(ppc_frame_apply_summary(list(df1, df2), "mode"))
})

test_that("ppc_frame_apply_cri computes quantiles", {
  # multi-row data.frames
  df1 <- data.frame(a = c(1, 10))
  df2 <- data.frame(a = c(2, 20))
  df3 <- data.frame(a = c(3, 30))

  upper <- ppc_frame_apply_cri(list(df1, df2, df3),
                               cri_level = 0.95, upper = TRUE)
  lower <- ppc_frame_apply_cri(list(df1, df2, df3),
                               cri_level = 0.95, upper = FALSE)

  expect_true(all(upper$a >= lower$a))

  # named vector path: all-NA returns NA
  v1 <- c(x = NA_real_)
  v2 <- c(x = NA_real_)
  result <- ppc_frame_apply_cri(list(v1, v2))
  expect_true(is.na(result[["x"]]))
})

test_that("Format-credible-interval works correctly on vectors", {
  # basic formatting
  result <- format_credible_interval(1.23, 0.9, 1.5)
  expect_true(is.character(result))
  expect_true(grepl("\\[", result))

  # as_integer
  expect_equal(
    format_credible_interval(2.1, 1.7, 2.3,
                             as_integer = TRUE,
                             fmt = "%s-%s-%s"),
    "2-2-2"
  )

  # NA mean returns NA_character_
  result <- format_credible_interval(NA, NA, NA)
  expect_true(is.na(result))

  # equal bounds
  result <- format_credible_interval(2, 2, 2)
  expect_true(is.character(result))

  # inconsistent interval warns
  expect_warning(format_credible_interval(5, 10, 20),
                 regexp = "inconsistent")

  # mismatched lengths error
  expect_error(format_credible_interval(c(1, 2), c(1), c(2)),
               regexp = "same length")

  # wrong str_format error
  expect_error(format_credible_interval(1, 0, 2, fmt = "%s only one"),
               regexp = "three.*placeholder")

  # data.frame input
  mn <- data.frame(a = c(1, 2), b = c(3, 4))
  lo <- data.frame(a = c(0, 1), b = c(2, 3))
  up <- data.frame(a = c(2, 3), b = c(4, 5))
  result <- format_credible_interval(mn, lo, up)
  expect_true(is.data.frame(result))
  expect_equal(dim(result), dim(mn))

  # mismatched data.frame dimensions
  lo_bad <- data.frame(a = c(0))
  expect_error(format_credible_interval(mn, lo_bad, up),
               regexp = "same dimensions")

  # mean == 0 branch for digit calculation
  result <- format_credible_interval(0, -0.1, 0.1)
  expect_true(is.character(result))
})

test_that("PPC extract-predicted-results returns correct structure", {
  df1 <- data.frame(a = c(1, 2))
  df2 <- data.frame(a = c(3, 4))
  df3 <- data.frame(a = c(5, 6))
  result <- ppc_extract_predicted_results(list(df1, df2, df3), cri = 0.95)
  expect_true(is.list(result))
  expect_true(all(c("predicted", "cri_lower", "cri_upper") %in% names(result)))
  expect_true(is.data.frame(result$predicted))
})

test_that("Helper to count non-zero-cases behaves correctly", {
  df <- data.frame(
    house_id = c(rep("hh1", 2), rep("hh2", 3), rep("hh3", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0)
  )
  df$id <- paste(df$house_id, seq_len(dim(df)[2]), sep = "_")
  result <- helper_count_non_zero_cases(df)
  expect_true(is.data.frame(result))
  expect_true(all(c("house_size", "count") %in% colnames(result)))
  # hh1 (size 2, 2 cases) and hh2 (size 3, 1 case) and hh3 (size 2, 1 case)
  expect_equal(result$count[result$house_size == 2], 2)
  expect_equal(result$count[result$house_size == 3], 1)
})

test_that("Helper to filter tables by max_house_size behaves correctly", {
  obs <- data.frame(house_size = 2:4,
                    `0` = c(1, 2, 3), `1` = c(4, 5, 0), `2` = c(0, 1, 0),
                    check.names = FALSE)
  pred <- obs * 1.1
  pred$house_size <- obs$house_size
  cri_lo <- obs * 0.9
  cri_lo$house_size <- obs$house_size
  cri_up <- obs * 1.2
  cri_up$house_size <- obs$house_size

  tbl <- list(observed = obs, predicted = pred,
              cri_lower = cri_lo, cri_upper = cri_up)

  # no filter
  result <- helper_filter_tables(tbl, NA)
  expect_equal(nrow(result$observed), 3)

  # filter to max size 3
  result2 <- helper_filter_tables(tbl, 3)
  expect_true(all(result2$observed$house_size <= 3))
})

test_that("Helper to combine observations and predictions interleaves rows", {
  obsdf <- data.frame(house_size = c(2, 3), `0` = c(1, 2), `1` = c(3, 4),
                      check.names = FALSE)
  preddf <- data.frame(`0` = c(1.1, 2.1), `1` = c(3.1, 4.1),
                       check.names = FALSE)
  result <- helper_combine_obs_pred(obsdf, preddf)
  expect_true(is.data.frame(result))
  expect_true("distr" %in% colnames(result))
  expect_true("house_size" %in% colnames(result))
  # alternating rows
  expect_equal(result$distr, rep(c("observed", "predicted"), 2))
  expect_equal(nrow(result), 4)
})

test_that("HousepiPrediction active fields are read-only", {
  hhmodel <- local_ppc$hhmodel
  results <- local_ppc$results
  suppressWarnings(
    ppc <- predicted_measures(hhmodel, results,
                              sar = measure_sar_by_house_size,
                              progress = FALSE)
  )
  expect_error(ppc$point_estimate <- "La", regexp = "read only")
  expect_error(ppc$cri_level <- "Li", regexp = "read only")
  expect_error(ppc$nrep <- "Lu", regexp = "read only")
  expect_error(ppc$thin_ess <- "Le", regexp = "read only")
  expect_error(ppc$is_posterior <- "Lo", regexp = "read only")
  expect_error(ppc$model_name <- "Patriots...", regexp = "read only")
  # read access works
  expect_silent(ppc$point_estimate)
  expect_silent(ppc$cri_level)
  expect_silent(ppc$nrep)
  expect_silent(ppc$thin_ess)
  expect_silent(ppc$is_posterior)
  expect_silent(ppc$model_name)
})

test_that("Predicted-measures raises error when no measure is provided", {
  expect_error(predicted_measures(local_ppc$hhmodel, local_ppc$results,
                                  progress = FALSE),
               regexp = "No measure was provided")
})

test_that("Predicted cases by house size works as expected", {
  hhmodel <- local_ppc$hhmodel
  results <- local_ppc$results
  suppressWarnings(
    ppc <- predicted_cases_by_house_size(hhmodel, results, progress = FALSE)
  )
  expect_s3_class(ppc, "PredictedCasesBySize")

  # formatted_sec_attack_rate
  sar <- ppc$formatted_sec_attack_rate
  expect_true(is.data.frame(sar))
  expect_true("house_size" %in% colnames(sar))
  expect_true("observed" %in% colnames(sar))

  # formatted_cases_table
  ct <- ppc$formatted_cases_table
  expect_true(is.list(ct))
  expect_true(length(ct) == 3)

  # subclass active fields and inherited dynamic measure fields are read-only
  expect_error(ppc$formatted_sec_attack_rate <- "La", regexp = "read only")
  expect_error(ppc$formatted_cases_table <- "Li", regexp = "read only")
  expect_error(ppc$sec_attack_rate <- "Lu", regexp = "read only")
  expect_error(ppc$cases_table_by_size <- "Le", regexp = "read only")

  # set_max_house_size and drop
  ppc$set_max_house_size(3)
  sar2 <- ppc$formatted_sec_attack_rate
  expect_true(all(sar2$house_size <= 3))
  ppc$drop_max_house_size()

  # chisq_test
  chi <- ppc$chisq_test()
  expect_true(is.list(chi))
  expect_true(all(c("chisq_value", "ddof", "p_value") %in% names(chi)))
})

test_that("PPC does not pollute the engine's data_wide between calls", {
  # repeated PPC calls must not mutate the model's data_wide
  hhmodel <- local_ppc$hhmodel
  results <- local_ppc$results

  dw_before <- hhmodel$data_wide
  suppressWarnings({
    ppc1 <- hhmodel |>
      predicted_cases_by_house_size(results, progress = FALSE,
                                    thin_ess = FALSE, max_nrep = 5)
  })
  expect_identical(hhmodel$data_wide, dw_before)

  suppressWarnings({
    ppc2 <- hhmodel |>
      predicted_cases_by_house_size(results, progress = FALSE,
                                    thin_ess = FALSE, max_nrep = 5)
  })
  expect_identical(ppc1$cases_table_by_size$observed,
                   ppc2$cases_table_by_size$observed)
  expect_identical(ppc1$sec_attack_rate$observed,
                   ppc2$sec_attack_rate$observed)
})

test_that("ppc_observed_outside_cri flags out-of-range numeric entries", {
  # numeric vectors: middle is below, last is above, NA propagates as FALSE
  obs <- c(1, 5, 10, NA)
  low <- c(0, 6, 0, 0)
  up  <- c(2, 8, 9, 10)
  expect_equal(ppc_observed_outside_cri(obs, low, up),
               c(FALSE, TRUE, TRUE, FALSE))

  # non-numeric input is treated as not-a-mismatch (skipped)
  expect_equal(ppc_observed_outside_cri(c("a", "b"), c("a", "b"), c("a", "b")),
               c(FALSE, FALSE))

  # data.frame: numeric column flagged, non-numeric column ignored
  obs_df <- data.frame(x = c(1, 5), label = c("a", "b"))
  low_df <- data.frame(x = c(0, 6), label = c("a", "b"))
  up_df  <- data.frame(x = c(2, 8), label = c("a", "b"))
  res <- ppc_observed_outside_cri(obs_df, low_df, up_df)
  expect_true(is.data.frame(res))
  expect_equal(res$x, c(FALSE, TRUE))
  expect_equal(res$label, c(FALSE, FALSE))
})

test_that("HousepiPrediction warns on observed values outside CrI", {
  # parameter set deliberately so that SAR~1
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("chain_binomial") |>
    add_parameter_template()
  param_table <- data.frame(alpha = 50, beta = 50, gamma = 1, size_ref = 4)
  expect_warning(
    predicted_measures(hhmodel, param_table,
                       sar = measure_sar_by_house_size,
                       progress = FALSE, thin_ess = FALSE, max_nrep = 20),
    regexp = ".*outside.*credible interval"
  )
})

test_that("predicted_measures works with custom measure", {
  hhmodel <- local_ppc$hhmodel
  results <- local_ppc$results
  suppressWarnings({
    ppc <- predicted_measures(hhmodel, results,
                              sar = measure_sar_by_house_size,
                              progress = FALSE,
                              thin_ess = FALSE,
                              max_nrep = 10)
  })
  expect_s3_class(ppc, "HousepiPrediction")
  expect_true(!is.null(ppc$sar))
  expect_true("observed" %in% names(ppc$sar))
  expect_true("predicted" %in% names(ppc$sar))
  expect_snapshot(print(ppc))
  # max_nrep caps the number of replicates used for predictions
  expect_identical(ppc$nrep, 10L)
})
