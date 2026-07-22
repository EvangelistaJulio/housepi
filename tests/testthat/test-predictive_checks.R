library(testthat)

# Shared fixture: a freshly templated flu_1918_frederick chain-binomial model
# and one inference run, reused across PPC tests that don't mutate either.
mk_ppc <- local({
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter_template()
  results <- suppressWarnings(run_inference(
    hhmodel, steps = 100, chains = 1, adapt = 0, thin = 1, progress = FALSE
  ))
  list(hhmodel = hhmodel, results = results)
})

test_that("measure-sar-by-hhsize computes SAR correctly", {
  # empty call: columns c(hhsize, 0))
  measure_sar_by_hhsize(NULL) |>
    dim() |>
    expect_equal(c(0L, 2L))

  # normal case
  data_wide <- data.frame(
    id = paste0("p", 1:7),
    house_id = c(rep("hh1", 2), rep("hh2", 3), rep("hh3", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0)
  )

  result <- measure_sar_by_hhsize(data_wide)
  expect_true(is.data.frame(result))
  expect_true(all(c("hhsize", "sar") %in% colnames(result)))
  # size 2: mean of SAR 1.0 (hh1) and 0.0 (hh3) = 0.5
  expect_equal(result[result$hhsize == 2, "sar"], 0.5)
  # size 3: SAR = (1-1)/(3-1) = 0
  expect_equal(result[result$hhsize == 3, "sar"], 0)

  # FIXME: what should be the behavior here?
  #data_wide$house_id <- NA
  #result <- measure_sar_by_hhsize(data_wide)
})

test_that("measure_ar_by_hhsize computes attack rate correctly", {
  # empty call: columns c(hhsize, attack_rate)
  measure_ar_by_hhsize(NULL) |>
    dim() |>
    expect_equal(c(0L, 2L))

  # normal case
  data_wide <- data.frame(
    id = paste0("p", 1:7),
    house_id = c(rep("hh1", 2), rep("hh2", 3), rep("hh3", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0)
  )

  result <- measure_ar_by_hhsize(data_wide)
  expect_true(is.data.frame(result))
  expect_true(all(c("hhsize", "ar") %in% colnames(result)))
  # size 2: mean of AR 2/2 (hh1) and 1/2 (hh3) = 0.75
  expect_equal(result[result$hhsize == 2, "ar"], 0.75)
  # size 3: AR = 1/3
  expect_equal(result[result$hhsize == 3, "ar"], 1 / 3)
})

test_that("measure_cases_by_hhsize produces correct contingency table", {
  # empty call: columns c(hhsize, 0))
  measure_cases_by_hhsize(NULL) |>
    dim() |>
    expect_equal(c(0L, 2L))

  data_wide <- data.frame(id = paste0("p", 1:5),
                          house_id = c(rep("hh1", 2), rep("hh2", 3)),
                          status = c(1, 1, 1, 0, 0))

  result <- measure_cases_by_hhsize(data_wide)

  expect_true(is.data.frame(result))
  expect_true("hhsize" %in% colnames(result))
  # grid should be complete
  expect_true(nrow(result) > 0)

})

test_that("validate_cri_level validates credible interval", {
  expect_equal(validate_cri_level(0.95), 0.95)
  expect_equal(validate_cri_level(0.8), 0.8)
  # boundary and invalid values
  validate_cri_level(0.5) |>
    expect_error("between 0.5 and 1")
  validate_cri_level(1.0) |>
    expect_error("between 0.5 and 1")
  validate_cri_level(0.3) |>
    expect_error("between 0.5 and 1")
  validate_cri_level(NA_real_) |>
    expect_error("non-NA")
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
  format_credible_interval(2.1, 1.7, 2.3,
                           as_integer = TRUE,
                           fmt = "%s-%s-%s") |>
    expect_equal("2-2-2")

  # NA mean returns NA_character_
  result <- format_credible_interval(NA, NA, NA)
  expect_true(is.na(result))

  # equal bounds
  result <- format_credible_interval(2, 2, 2)
  expect_true(is.character(result))

  # inconsistent interval warns
  format_credible_interval(5, 10, 20) |>
    expect_warning("inconsistent")

  # mismatched lengths error
  format_credible_interval(c(1, 2), c(1), c(2)) |>
    expect_error("same length")

  # wrong str_format error
  format_credible_interval(1, 0, 2, fmt = "%s only one") |>
    expect_error("three.*placeholder")

  # data.frame input
  mn <- data.frame(a = c(1, 2), b = c(3, 4))
  lo <- data.frame(a = c(0, 1), b = c(2, 3))
  up <- data.frame(a = c(2, 3), b = c(4, 5))
  result <- format_credible_interval(mn, lo, up)
  expect_true(is.data.frame(result))
  expect_equal(dim(result), dim(mn))

  # mismatched data.frame dimensions
  lo_bad <- data.frame(a = c(0))
  format_credible_interval(mn, lo_bad, up) |>
    expect_error("same dimensions")

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

test_that("Helper to count houses by household size behaves correctly", {
  df <- data.frame(
    house_id = c(rep("hh1", 2), rep("hh2", 3), rep("hh3", 2), rep("hh4", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0, 0, 0)
  )
  df$id <- paste(df$house_id, seq_len(dim(df)[2]), sep = "_")
  result <- helper_house_counts(df)
  expect_true(is.data.frame(result))
  expect_true(all(c("hhsize", "total", "at_least_1_case") %in% colnames(result)))
  # size 2: hh1 (2 cases), hh3 (1 case), hh4 (0 cases) -> 3 total, 2 affected
  expect_equal(result$total[result$hhsize == 2], 3)
  expect_equal(result$at_least_1_case[result$hhsize == 2], 2)
  # size 3: hh2 (1 case) -> 1 total, 1 affected
  expect_equal(result$total[result$hhsize == 3], 1)
  expect_equal(result$at_least_1_case[result$hhsize == 3], 1)
})

test_that("Helper to filter tables by max_hhsize behaves correctly", {
  obs <- data.frame(hhsize = 2:4,
                    `0` = c(1, 2, 3), `1` = c(4, 5, 0), `2` = c(0, 1, 0),
                    check.names = FALSE)
  pred <- obs * 1.1
  pred$hhsize <- obs$hhsize
  cri_lo <- obs * 0.9
  cri_lo$hhsize <- obs$hhsize
  cri_up <- obs * 1.2
  cri_up$hhsize <- obs$hhsize

  tbl <- list(observed = obs, predicted = pred,
              cri_lower = cri_lo, cri_upper = cri_up)

  # no filter
  result <- helper_filter_tables(tbl)
  expect_equal(nrow(result$observed), 3)

  # filter to max size 3
  result2 <- helper_filter_tables(tbl, max_hhsize = 3)
  expect_true(all(result2$observed$hhsize <= 3))

  # filter to min size 3
  result3 <- helper_filter_tables(tbl, min_hhsize = 3)
  expect_true(all(result3$observed$hhsize >= 3))

  # filter to a size band [3, 3]
  result4 <- helper_filter_tables(tbl, min_hhsize = 3, max_hhsize = 3)
  expect_true(all(result4$observed$hhsize == 3))
})

test_that("Helper to combine observations and predictions interleaves rows", {
  obsdf <- data.frame(hhsize = c(2, 3), `0` = c(1, 2), `1` = c(3, 4),
                      check.names = FALSE)
  preddf <- data.frame(`0` = c(1.1, 2.1), `1` = c(3.1, 4.1),
                       check.names = FALSE)
  result <- helper_combine_obs_pred(obsdf, preddf)
  expect_true(is.data.frame(result))
  expect_true("distr" %in% colnames(result))
  expect_true("hhsize" %in% colnames(result))
  # alternating rows
  expect_equal(result$distr, rep(c("observed", "predicted"), 2))
  expect_equal(nrow(result), 4)
})

test_that("HousepiCustomPPC active fields are read-only", {
  # predicted-measures raises error when no measure is provided
  mk_ppc$hhmodel |>
    run_custom_predictive_check(mk_ppc$results) |>
    expect_error("No measure was provided")

  hhmodel <- mk_ppc$hhmodel
  results <- mk_ppc$results

  suppressWarnings(
    ppc <- hhmodel |>
      run_custom_predictive_check(
        results, sar = measure_sar_by_hhsize,
        progress = FALSE, thin_ess = FALSE, max_nrep = 10
      )
  )

  expect_s3_class(ppc, "HousepiCustomPPC")
  expect_true(!is.null(ppc$sar))
  expect_true("observed" %in% names(ppc$sar))
  expect_true("predicted" %in% names(ppc$sar))
  expect_snapshot(print(ppc))
  expect_identical(ppc$nrep, 10L)

  # active fields are read-only
  ppc |>
    helper_test_read_only(
      "point_estimate", "cri_level", "nrep",
      "thin_ess", "is_posterior", "model_name"
    )

  # read access works
  expect_silent(ppc$point_estimate)
  expect_silent(ppc$cri_level)
  expect_silent(ppc$nrep)
  expect_silent(ppc$thin_ess)
  expect_silent(ppc$is_posterior)
  expect_silent(ppc$model_name)
})

test_that("Predicted cases by house size works as expected", {
  hhmodel <- mk_ppc$hhmodel
  results <- mk_ppc$results
  suppressWarnings(
    ppc <- run_epi_predictive_check(hhmodel, results)
  )
  expect_s3_class(ppc, "HousepiPPC")

  # formatted_attack_rate: counts are the raw house totals on the AR grid
  ar <- ppc$formatted_attack_rate
  expect_true(is.data.frame(ar))
  expect_true(all(c("hhsize", "observed", "counts") %in% colnames(ar)))
  hc <- helper_house_counts(hhmodel$data_wide)
  expect_equal(ar$counts, hc$total[match(ar$hhsize, hc$hhsize)])

  # formatted_sec_attack_rate: counts are houses with at least one case
  sar <- ppc$formatted_sec_attack_rate
  expect_true(is.data.frame(sar))
  expect_true("hhsize" %in% colnames(sar))
  expect_true("observed" %in% colnames(sar))
  expect_equal(sar$counts, hc$at_least_1_case[match(sar$hhsize, hc$hhsize)])

  # formatted_cases_count
  ct <- ppc$formatted_cases_count
  expect_true(is.list(ct))
  expect_true(length(ct) == 3)

  # subclass active fields and inherited dynamic measure fields are read-only
  ppc |>
    helper_test_read_only(
      "formatted_attack_rate", "attack_rate",
      "formatted_sec_attack_rate", "sec_attack_rate",
      "formatted_cases_count", "cases_count"
    )

  # set_max_hhsize and drop
  ppc$set_max_hhsize(3)
  sar2 <- ppc$formatted_sec_attack_rate
  expect_true(all(sar2$hhsize <= 3))
  ppc$drop_max_hhsize()

  # set_min_hhsize and drop
  ppc$set_min_hhsize(3)
  sar3 <- ppc$formatted_sec_attack_rate
  expect_true(all(sar3$hhsize >= 3))
  ar3 <- ppc$formatted_attack_rate
  expect_true(all(ar3$hhsize >= 3))
  ppc$drop_min_hhsize()

  # chisq_test
  chi <- ppc$chisq_test()
  expect_true(is.list(chi))
  expect_true(all(c("chisq_value", "ddof", "p_value") %in% names(chi)))
})

test_that("run_epi_predictive_check honours min/max_hhsize arguments", {
  hhmodel <- mk_ppc$hhmodel
  results <- mk_ppc$results
  suppressWarnings(
    ppc <- hhmodel |>
      run_epi_predictive_check(results, min_hhsize = 2L, max_hhsize = 3L)
  )
  expect_true(all(ppc$formatted_attack_rate$hhsize >= 2))
  expect_true(all(ppc$formatted_attack_rate$hhsize <= 3))
  expect_true(all(ppc$formatted_sec_attack_rate$hhsize >= 2))
  expect_true(all(ppc$formatted_sec_attack_rate$hhsize <= 3))
})

test_that("PPC does not pollute the engine's data_wide between calls", {
  # repeated PPC calls must not mutate the model's data_wide
  hhmodel <- mk_ppc$hhmodel
  results <- mk_ppc$results

  dw_before <- hhmodel$data_wide
  suppressWarnings({
    ppc1 <- hhmodel |>
      run_epi_predictive_check(results, thin_ess = FALSE, max_nrep = 5)
  })
  expect_identical(hhmodel$data_wide, dw_before)

  suppressWarnings({
    ppc2 <- hhmodel |>
      run_epi_predictive_check(results, thin_ess = FALSE, max_nrep = 5)
  })
  expect_identical(ppc1$cases_count$observed,
                   ppc2$cases_count$observed)
  expect_identical(ppc1$sec_attack_rate$observed,
                   ppc2$sec_attack_rate$observed)
})

test_that("ppc_observed_outside_cri flags out-of-range numeric entries", {
  # numeric vectors: middle is below, last is above, NA propagates as FALSE
  list(observed = c(1, 5, 10, NA),
       cri_lower = c(0, 6, 0, 0),
       cri_upper = c(2, 8, 9, 10)) |>
    ppc_observed_outside_cri() |>
    expect_equal(c(FALSE, TRUE, TRUE, FALSE))

  # non-numeric input is treated as not-a-mismatch (skipped)
  list(observed = c("a", "b"),
       cri_lower = c("a", "b"),
       cri_upper = c("a", "b")) |>
    ppc_observed_outside_cri() |>
    expect_equal(c(FALSE, FALSE))

  # data.frame: numeric column flagged, non-numeric column ignored
  res <- list(observed = data.frame(x = c(1, 5), label = c("a", "b")),
              cri_lower = data.frame(x = c(0, 6), label = c("a", "b")),
              cri_upper = data.frame(x = c(2, 8), label = c("a", "b"))) |>
    ppc_observed_outside_cri()
  expect_true(is.data.frame(res))
  expect_equal(res$x, c(FALSE, TRUE))
  expect_equal(res$label, c(FALSE, FALSE))
})

test_that("HousepiCustomPPC warns on observed values outside CrI", {
  # parameter set deliberately so that SAR~1
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final_size") |>
    add_parameter_template() |>
    run_custom_predictive_check(
      data.frame(alpha = 50., beta = 50., gamma = 1., hhsize_ref = 4L),
      sar = measure_sar_by_hhsize,
      thin_ess = FALSE, max_nrep = 20L
    ) |>
    expect_warning(".*outside.*CrI")
})
