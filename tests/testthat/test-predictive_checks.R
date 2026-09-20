library(testthat)

# templated chain-binomial model and one short inference run on it
mk_ppc <- local({
  hhmodel <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  results <- hhmodel |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings()
  list(model = hhmodel, results = results)
})

# three households of sizes two, three and two
mk_houses <- local({
  hhdata <- data.frame(
    id = paste0("m", 1:7),
    hhid = c(rep("h1", 2), rep("h2", 3), rep("h3", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0)
  )
  list(data = hhdata)
})

######################################
###  PREDICTIVE CHECKS : MEASURES  ###
######################################

test_that("Measures require the columns they aggregate over", {
  mk_houses$data[c("id", "hhid")] |>
    measure_sar_by_hhsize() |>
    expect_error("status")
  
  # without a household identifier every row would count as its own household
  mk_houses$data[c("id", "status")] |>
    measure_sar_by_hhsize() |>
    expect_error("hhid")
})

test_that("Attack rate is averaged per household size", {
  # empty call: columns c(hhsize, ar)
  measure_ar_by_hhsize(NULL) |>
    dim() |>
    expect_equal(c(0L, 2L))
  
  result <- measure_ar_by_hhsize(mk_houses$data)
  expect_named(result, c("hhsize", "ar"), ignore.order = TRUE)
  # size 2: mean of AR 2/2 (h1) and 1/2 (h3) = 0.75
  expect_equal(result[result$hhsize == 2, "ar"], 0.75)
  # size 3: AR = 1/3
  expect_equal(result[result$hhsize == 3, "ar"], 1 / 3)
})

test_that("Secondary attack rate is averaged per household size", {
  # empty call: columns c(hhsize, sar)
  measure_sar_by_hhsize(NULL) |>
    dim() |>
    expect_equal(c(0L, 2L))

  result <- measure_sar_by_hhsize(mk_houses$data)
  expect_named(result, c("hhsize", "sar"), ignore.order = TRUE)
  # size 2: mean of SAR 1.0 (h1) and 0.0 (h3) = 0.5
  expect_equal(result[result$hhsize == 2, "sar"], 0.5)
  # size 3: SAR = (1-1)/(3-1) = 0
  expect_equal(result[result$hhsize == 3, "sar"], 0)
})

test_that("Cases are tabulated on the complete household size grid", {
  # empty call: columns c(hhsize, 0)
  measure_cases_by_hhsize(NULL) |>
    dim() |>
    expect_equal(c(0L, 2L))

  result <- measure_cases_by_hhsize(mk_houses$data)
  expect_true("hhsize" %in% colnames(result))
  expect_setequal(result$hhsize, c(2, 3))
})

test_that("Measures keep households that all share one size", {
  hhdata <- data.frame(id = paste0("m", 1:3), hhid = "h1",
                       status = c(1, 1, 0))
  # AR = 2/3, SAR = (2-1)/(3-1)
  measure_ar_by_hhsize(hhdata) |>
    expect_equal(data.frame(hhsize = 3L, ar = 2 / 3))
  measure_sar_by_hhsize(hhdata) |>
    expect_equal(data.frame(hhsize = 3L, sar = 0.5))
  expect_equal(measure_cases_by_hhsize(hhdata)$hhsize, 3L)
})

test_that("Household size geometric mean is weighted by secondary cases", {
  # no household with a secondary case
  measure_geom_mean_hhsize(NULL) |>
    expect_equal(NA_real_)

  # only h1 (size 2, 2 cases) carries weight
  measure_geom_mean_hhsize(mk_houses$data) |>
    expect_equal(2)

  # h4 (size 4, 3 cases) weighs 2: exp((log 2 + 2 log 4) / 3) = 2^(5/3)
 mk_houses$data |> 
   rbind(data.frame(
     id = paste0("m", 8:11), hhid = "h4", status = c(1, 1, 1, 0)
     )) |>
   measure_geom_mean_hhsize() |>
    expect_equal(2^(5 / 3))
})

test_that("Measures decode a labelled status column", {
  ref_ar <- measure_ar_by_hhsize(mk_houses$data)
  ref_sar <- measure_sar_by_hhsize(mk_houses$data)

  # a factor status is decoded through its levels, not through level indices
  as_factor <- mk_houses$data
  as_factor$status = factor(ifelse(as_factor$status == 1,
                                   "symptomatic", "escaped"))
  expect_equal(measure_sar_by_hhsize(as_factor), ref_sar)
  expect_equal(measure_ar_by_hhsize(as_factor), ref_ar)

  # a character status is decoded, not coerced to NA
  as_char <- as_factor
  as_char$status <- as.character(as_factor$status)
  expect_equal(measure_sar_by_hhsize(as_char), ref_sar)
  expect_equal(measure_ar_by_hhsize(as_char), ref_ar)

  # an unnamed negative level warns, and still returns the (all-case) measure
  unnamed <- mk_houses$data
  unnamed$status <- ifelse(unnamed$status == 1, "positive", "negative")
  expect_warning(
    all_cases <- measure_ar_by_hhsize(unnamed),
    "escaped label"
  )
  expect_true(all(all_cases$ar == 1))

  # naming it silences the guard and recovers the coded result
  measure_ar_by_hhsize(unnamed, escaped = "negative") |>
    expect_equal(ref_ar)
})

##########################################
## PREDICTIVE CHECKS : SUMMARY HELPERS  ##
##########################################

test_that("Credible interval level is validated", {
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

test_that("Replicate frames are combined element-wise", {
  # list of data.frames
  result <- ppc_frame_apply(list(data.frame(a = c(1, 2), b = c(3, 4)),
                                 data.frame(a = c(5, 6), b = c(7, 8))), mean)
  expect_equal(result$a, c(3, 4))
  expect_equal(result$b, c(5, 6))

  # named vectors
  ppc_frame_apply(list(c(x = 1, y = 2), c(x = 3, y = 4)), mean) |>
    expect_equal(c(x = 2, y = 3))
})

test_that("Replicate summaries handle mean, median, NA and non-numeric", {
  # multi-row data.frames required by internal apply() on matrix path
  frame_1 <- data.frame(a = c(1, 2), b = c(NA_real_, NA_real_))
  frame_2 <- data.frame(a = c(3, 4), b = c(NA_real_, NA_real_))

  result <- ppc_frame_apply_summary(list(frame_1, frame_2), "mean")
  expect_equal(result$a, c(2, 3))
  expect_true(all(is.na(result$b)))

  # median path
  frame_3 <- data.frame(a = c(10, 20), b = c(NA_real_, NA_real_))
  ppc_frame_apply_summary(list(frame_1, frame_2, frame_3), "median")$a |>
    expect_equal(c(3, 4))

  # named vector path: non-numeric uniform values
  result <- ppc_frame_apply_summary(list(c(x = 1, y = NA_real_),
                                         c(x = 3, y = NA_real_)), "mean")
  expect_equal(result[["x"]], 2)
  expect_true(is.na(result[["y"]]))

  # invalid summary stat
  ppc_frame_apply_summary(list(frame_1, frame_2), "mode") |>
    expect_error("one between")
})

test_that("Replicate credible intervals bound each other", {
  frames <- list(data.frame(a = c(1, 10)), data.frame(a = c(2, 20)),
                 data.frame(a = c(3, 30)))
  upper <- ppc_frame_apply_cri(frames, cri_level = 0.95, upper = TRUE)
  lower <- ppc_frame_apply_cri(frames, cri_level = 0.95, upper = FALSE)
  expect_true(all(upper$a >= lower$a))

  # named vector path: all-NA returns NA
  ppc_frame_apply_cri(list(c(x = NA_real_), c(x = NA_real_)))[["x"]] |>
    is.na() |>
    expect_true()
})

test_that("Credible intervals are formatted on vectors and frames", {
  # basic formatting
  expect_match(format_credible_interval(1.23, 0.9, 1.5), "\\[")

  # as_integer
  format_credible_interval(2.1, 1.7, 2.3, as_integer = TRUE,
                           fmt = "%s-%s-%s") |>
    expect_equal("2-2-2")

  # NA mean returns NA_character_
  expect_equal(format_credible_interval(NA, NA, NA), NA_character_)

  # equal bounds, and a zero mean for the digit count
  expect_type(format_credible_interval(2, 2, 2), "character")
  expect_type(format_credible_interval(0, -0.1, 0.1), "character")

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
  mean_frame <- data.frame(a = c(1, 2), b = c(3, 4))
  upper_frame <- data.frame(a = c(2, 3), b = c(4, 5))
  result <- format_credible_interval(
    mean_frame, data.frame(a = c(0, 1), b = c(2, 3)), upper_frame
  )
  expect_s3_class(result, "data.frame")
  expect_equal(dim(result), dim(mean_frame))

  # mismatched data.frame dimensions
  format_credible_interval(mean_frame, data.frame(a = c(0)), upper_frame) |>
    expect_error("same dimensions")
})

test_that("Predicted results carry the point estimate and its interval", {
  result <- ppc_extract_predicted_results(
    list(data.frame(a = c(1, 2)), data.frame(a = c(3, 4)),
         data.frame(a = c(5, 6))),
    cri = 0.95
  )
  expect_true(all(c("predicted", "cri_lower", "cri_upper") %in% names(result)))
  expect_s3_class(result$predicted, "data.frame")
})

test_that("Houses are counted by household size", {
  hhdata <- data.frame(
    hhid = c(rep("h1", 2), rep("h2", 3), rep("h3", 2), rep("h4", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0, 0, 0)
  )
  hhdata$id <- paste0("m", seq_len(nrow(hhdata)))
  result <- helper_house_counts(hhdata)
  expect_true(all(c("hhsize", "total", "at_least_1_case") %in%
                    colnames(result)))
  # size 2: h1 (2 cases), h3 (1 case), h4 (0 cases) -> 3 total, 2 affected
  expect_equal(result$total[result$hhsize == 2], 3)
  expect_equal(result$at_least_1_case[result$hhsize == 2], 2)
  # size 3: h2 (1 case) -> 1 total, 1 affected
  expect_equal(result$total[result$hhsize == 3], 1)
  expect_equal(result$at_least_1_case[result$hhsize == 3], 1)
})

test_that("Tables are filtered by household size band", {
  # one case column per count a household of the largest size can hold
  observed <- data.frame(hhsize = 2:4,
                         `0` = c(1, 2, 3), `1` = c(4, 5, 0), `2` = c(0, 1, 0),
                         `3` = c(0, 1, 1), `4` = c(0, 0, 1),
                         check.names = FALSE)
  scaled <- function(k) transform(observed * k, hhsize = observed$hhsize)
  tables <- list(observed = observed, predicted = scaled(1.1),
                 cri_lower = scaled(0.9), cri_upper = scaled(1.2))

  # no filter
  expect_equal(nrow(helper_filter_tables(tables)$observed), 3)
  # upper, lower and both bounds
  helper_filter_tables(tables, max_hhsize = 3)$observed$hhsize |>
    expect_equal(2:3)
  helper_filter_tables(tables, min_hhsize = 3)$observed$hhsize |>
    expect_equal(3:4)
  helper_filter_tables(tables, min_hhsize = 3, max_hhsize = 3) |>
    (\(x) x$observed$hhsize)() |>
    expect_equal(3L)
})

test_that("Observations and predictions are interleaved row by row", {
  result <- helper_combine_obs_pred(
    data.frame(hhsize = c(2, 3), `0` = c(1, 2), `1` = c(3, 4),
               check.names = FALSE),
    data.frame(`0` = c(1.1, 2.1), `1` = c(3.1, 4.1), check.names = FALSE)
  )
  expect_true(all(c("distr", "hhsize") %in% colnames(result)))
  expect_equal(result$distr, rep(c("observed", "predicted"), 2))
})

test_that("Observed entries outside the interval are flagged", {
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
  result <- list(observed = data.frame(x = c(1, 5), label = c("a", "b")),
                 cri_lower = data.frame(x = c(0, 6), label = c("a", "b")),
                 cri_upper = data.frame(x = c(2, 8), label = c("a", "b"))) |>
    ppc_observed_outside_cri()
  expect_s3_class(result, "data.frame")
  expect_equal(result$x, c(FALSE, TRUE))
  expect_equal(result$label, c(FALSE, FALSE))
})

##########################################
## PREDICTIVE CHECKS : END-TO-END       ##
##########################################

test_that("HousepiCustomPPC active fields are read-only", {
  # predicted-measures raises error when no measure is provided
  mk_ppc$model |>
    run_custom_predictive_check(mk_ppc$results) |>
    expect_error("No measure was provided")

  ppc <- mk_ppc$model |>
    run_custom_predictive_check(
      mk_ppc$results, sar = measure_sar_by_hhsize,
      thin_ess = FALSE, max_nrep = 10
    ) |>
    suppressWarnings()

  expect_s3_class(ppc, "HousepiCustomPPC")
  expect_true(all(c("observed", "predicted") %in% names(ppc$sar)))
  expect_snapshot(print(ppc))

  # read access returns the run settings
  expect_identical(ppc$nrep, 10L)
  expect_identical(ppc$point_estimate, "median")
  expect_identical(ppc$cri_level, 0.95)
  # the summary settings are read off the global options at construction
  withr::with_options(list(housepi.point_estimate = "mean",
                           housepi.cri_level = 0.8), {
    averaged <- mk_ppc$model |>
      run_custom_predictive_check(
        mk_ppc$results, sar = measure_sar_by_hhsize,
        thin_ess = FALSE, max_nrep = 2
      ) |>
      suppressWarnings()
    expect_identical(averaged$point_estimate, "mean")
    expect_identical(averaged$cri_level, 0.8)
  })
  expect_false(ppc$thin_ess)
  expect_true(ppc$is_posterior)
  expect_identical(ppc$model_name, "final size")

  # active fields are read-only
  ppc |>
    helper_test_read_only(
      "point_estimate", "cri_level", "nrep",
      "thin_ess", "is_posterior", "model_name"
    )
})

test_that("Predicted cases by house size works as expected", {
  ppc <- run_epi_predictive_check(mk_ppc$model, mk_ppc$results) |>
    suppressWarnings()
  expect_s3_class(ppc, "HousepiPPC")
  counts <- helper_house_counts(mk_ppc$model$data_wide)

  # formatted_attack_rate: counts are the raw house totals on the AR grid
  ar <- ppc$formatted_attack_rate
  expect_true(all(c("hhsize", "observed", "counts") %in% colnames(ar)))
  expect_equal(ar$counts, counts$total[match(ar$hhsize, counts$hhsize)])

  # formatted_sec_attack_rate: counts are houses with at least one case
  sar <- ppc$formatted_sec_attack_rate
  expect_true(all(c("hhsize", "observed") %in% colnames(sar)))
  expect_equal(sar$counts,
               counts$at_least_1_case[match(sar$hhsize, counts$hhsize)])

  # formatted_cases_count
  expect_length(ppc$formatted_cases_count, 3)

  # subclass active fields and inherited dynamic measure fields are read-only
  ppc |>
    helper_test_read_only(
      "formatted_attack_rate", "attack_rate",
      "formatted_sec_attack_rate", "sec_attack_rate",
      "formatted_cases_count", "cases_count"
    )

  # set_max_hhsize and drop
  ppc$set_max_hhsize(3)
  expect_true(all(ppc$formatted_sec_attack_rate$hhsize <= 3))
  ppc$drop_max_hhsize()

  # set_min_hhsize and drop
  ppc$set_min_hhsize(3)
  expect_true(all(ppc$formatted_sec_attack_rate$hhsize >= 3))
  expect_true(all(ppc$formatted_attack_rate$hhsize >= 3))
  ppc$drop_min_hhsize()

  # chisq_test
  expect_true(all(c("chisq_value", "ddof", "p_value") %in%
                    names(ppc$chisq_test())))
})

test_that("Epi PPC honours min/max_hhsize arguments", {
  ppc <- mk_ppc$model |>
    run_epi_predictive_check(mk_ppc$results,
                             min_hhsize = 2L, max_hhsize = 3L) |>
    suppressWarnings()
  expect_true(all(ppc$formatted_attack_rate$hhsize %in% 2:3))
  expect_true(all(ppc$formatted_sec_attack_rate$hhsize %in% 2:3))
})

test_that("Epi PPC observed values are stable across calls and status codes", {
  data_wide <- mk_ppc$model$data_wide
  run_on <- function(...) {
    mk_ppc$model |>
      run_epi_predictive_check(mk_ppc$results, thin_ess = FALSE,
                               max_nrep = 5, ...) |>
      suppressWarnings()
  }
  ppc <- run_on()
  # repeated PPC calls must not mutate the model's data_wide
  expect_identical(mk_ppc$model$data_wide, data_wide)

  # without an observation process the replicates stay engine-coded
  observed <- measure_ar_by_hhsize(data_wide, mk_ppc$model$escaped)
  expect_equal(ppc$formatted_attack_rate$observed, observed$ar)
  # with one, they are relabelled before the measure sees them
  run_on(obs_process = function(x) x)$formatted_attack_rate$observed |>
    expect_equal(observed$ar)

  # observed values agree across thread counts
  ppc_2 <- run_on(threads = 2L)
  expect_identical(ppc_2$cases_count$observed, ppc$cases_count$observed)
  expect_identical(ppc_2$sec_attack_rate$observed,
                   ppc$sec_attack_rate$observed)
  expect_identical(ppc_2$cases_table_by_size$observed,
                   ppc$cases_table_by_size$observed)
})

test_that("Custom PPC warns on observed values outside the interval", {
  # parameter set deliberately so that SAR~1
  mk_ppc$model |>
    run_custom_predictive_check(
      data.frame(alpha = 50., beta = 50., gamma = 1., hhsize_ref = 4L),
      sar = measure_sar_by_hhsize,
      thin_ess = FALSE, max_nrep = 20L
    ) |>
    expect_warning(".*outside.*CrI")
})
