library(testthat)

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

  # empty input errors
  expect_error(measure_sar_by_house_size(data.frame()),
               regexp = "No valid data_wide")
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

  # empty input errors
  expect_error(measure_cases_by_house_size(data.frame()),
               regexp = "No valid data_wide")
})

test_that("ppc_set_cri_level validates credible interval", {
  expect_equal(ppc_set_cri_level(0.95), 0.95)
  expect_equal(ppc_set_cri_level(0.8), 0.8)
  # boundary and invalid values
  expect_error(ppc_set_cri_level(0.5), regexp = "between 0.5 and 1")
  expect_error(ppc_set_cri_level(1.0), regexp = "between 0.5 and 1")
  expect_error(ppc_set_cri_level(0.3), regexp = "between 0.5 and 1")
  expect_error(ppc_set_cri_level(NA), regexp = "cannot be NA")
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

test_that("ppc_frame_apply_mean handles NA and non-numeric", {
  # multi-row data.frames required by internal apply() on matrix path
  df1 <- data.frame(a = c(1, 2), b = c(NA_real_, NA_real_))
  df2 <- data.frame(a = c(3, 4), b = c(NA_real_, NA_real_))

  result <- ppc_frame_apply_mean(list(df1, df2))
  expect_equal(result$a, c(2, 3))
  expect_true(all(is.na(result$b)))

  # named vector path: non-numeric uniform values
  v1 <- c(x = 1, y = NA_real_)
  v2 <- c(x = 3, y = NA_real_)
  result <- ppc_frame_apply_mean(list(v1, v2))
  expect_equal(result[["x"]], 2)
  expect_true(is.na(result[["y"]]))
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

test_that("format_credible_intervals formats vectors correctly", {
  # basic formatting
  result <- format_credible_intervals(c(1.23), c(0.9), c(1.5))
  expect_true(is.character(result))
  expect_true(grepl("\\[", result))

  # NA mean returns NA_character_
  result <- format_credible_intervals(NA, NA, NA)
  expect_true(is.na(result))

  # equal bounds

  result <- format_credible_intervals(2, 2, 2)
  expect_true(is.character(result))

  # inconsistent interval warns
  expect_warning(format_credible_intervals(5, 10, 20),
                 regexp = "inconsistent")

  # as_integer
  result <- format_credible_intervals(1.7, 1.2, 2.3, as_integer = TRUE)
  expect_true(is.character(result))

  # mismatched lengths error
  expect_error(format_credible_intervals(c(1, 2), c(1), c(2)),
               regexp = "same length")

  # wrong str_format error
  expect_error(format_credible_intervals(1, 0, 2, str_format = "%s only one"),
               regexp = "three.*placeholder")

  # data.frame input
  mn <- data.frame(a = c(1, 2), b = c(3, 4))
  lo <- data.frame(a = c(0, 1), b = c(2, 3))
  up <- data.frame(a = c(2, 3), b = c(4, 5))
  result <- format_credible_intervals(mn, lo, up)
  expect_true(is.data.frame(result))
  expect_equal(dim(result), dim(mn))

  # mismatched data.frame dimensions
  lo_bad <- data.frame(a = c(0))
  expect_error(format_credible_intervals(mn, lo_bad, up),
               regexp = "same dimensions")

  # mean == 0 branch for digit calculation
  result <- format_credible_intervals(0, -0.1, 0.1)
  expect_true(is.character(result))
})

test_that("PPC print helper formats messages correctly", {
  sms <- ppc_helper_print("test_model", TRUE, c("sar", "cases"), 0.95)
  expect_true(is.character(sms))
  expect_true(any(grepl("Posterior predictions", sms)))
  expect_true(any(grepl("test_model", sms)))
  expect_true(any(grepl("sar", sms)))
  expect_true(any(grepl("cases", sms)))
  expect_true(any(grepl("95%", sms)))

  # non-posterior
  sms2 <- ppc_helper_print("m", FALSE, "sar", 0.8)
  expect_true(any(grepl("Predictions", sms2)))
  expect_true(any(grepl("80%", sms2)))
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

test_that("CRI set name convention renames correctly", {
  x <- c("predicted", "cri_lower", "cri_upper")
  result <- cri_set_name_convention(x, cri = 0.95)
  expect_equal(result, c("mean", "q5", "q95"))

  result2 <- cri_set_name_convention(x, cri = 0.9)
  expect_equal(result2, c("mean", "q10", "q90"))
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
  expect_equal(result$distr, rep(c("Observed", "Predicted"), 2))
  expect_equal(nrow(result), 4)
})

test_that("Predicted-measures raises error when no measure is provided", {
  hhmodel <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    add_data_wide(helper_test_get_data("flu_1918_frederick"))
  results <- run_inference(hhmodel, steps = 5, progress = FALSE)
  expect_error(predicted_measures(hhmodel, results, progress = FALSE),
               regexp = "No measure was provided")
})

test_that("Predicted cases by house size works as expected", {
  hhmodel <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    add_data_wide(helper_test_get_data("flu_1918_frederick"))

  results <- run_inference(hhmodel, steps = 10, progress = FALSE)

  ppc <- predicted_cases_by_house_size(hhmodel, results, progress = FALSE)
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

  # set_max_house_size and drop
  ppc$set_max_house_size(3)
  sar2 <- ppc$formatted_sec_attack_rate
  expect_true(all(sar2$house_size <= 3))
  ppc$drop_max_house_size()

  # chisq_test
  chi <- ppc$chisq_test()
  expect_true(is.list(chi))
  expect_true(all(c("chisq_value", "ddof", "p_value") %in% names(chi)))

  # print
  expect_snapshot(print(ppc))
})

test_that("predicted_measures works with custom measure", {
  hhmodel <- new_housepi_model("chain_binomial") |>
    load_parameter_template() |>
    add_data_wide(helper_test_get_data("flu_1918_frederick"))

  results <- run_inference(hhmodel, steps = 10, progress = FALSE)

  ppc <- predicted_measures(hhmodel, results,
                            sar = measure_sar_by_house_size,
                            progress = FALSE,
                            ess = FALSE)
  expect_s3_class(ppc, "HousepiPrediction")
  expect_true(!is.null(ppc$sar))
  expect_true("observed" %in% names(ppc$sar))
  expect_true("predicted" %in% names(ppc$sar))
  expect_snapshot(print(ppc))
})
