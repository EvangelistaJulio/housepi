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
    expect_identical(c(0L, 2L))

  result <- measure_ar_by_hhsize(mk_houses$data)
  result |> expect_named(c("hhsize", "ar"), ignore.order = TRUE)
  # size 2: mean of AR 2/2 (h1) and 1/2 (h3) = 0.75
  result[result$hhsize == 2, "ar"] |> expect_equal(0.75)
  # size 3: AR = 1/3
  result[result$hhsize == 3, "ar"] |> expect_equal(1 / 3)
})

test_that("Secondary attack rate is averaged per household size", {
  # empty call: columns c(hhsize, sar)
  measure_sar_by_hhsize(NULL) |>
    dim() |>
    expect_identical(c(0L, 2L))

  result <- measure_sar_by_hhsize(mk_houses$data)
  result |> expect_named(c("hhsize", "sar"), ignore.order = TRUE)
  # size 2: mean of SAR 1.0 (h1) and 0.0 (h3) = 0.5
  result[result$hhsize == 2, "sar"] |> expect_equal(0.5)
  # size 3: SAR = (1-1)/(3-1) = 0
  result[result$hhsize == 3, "sar"] |> expect_equal(0)
})

test_that("Cases are tabulated on the complete household size grid", {
  # empty call: columns c(hhsize, 0)
  measure_cases_by_hhsize(NULL) |>
    dim() |>
    expect_identical(c(0L, 2L))

  result <- measure_cases_by_hhsize(mk_houses$data)
  ("hhsize" %in% colnames(result)) |> expect_true()
  result$hhsize |> expect_setequal(c(2, 3))
})

test_that("Measures keep households that all share one size", {
  hhdata <- data.frame(id = paste0("m", 1:3), hhid = "h1",
                       status = c(1, 1, 0))
  # AR = 2/3, SAR = (2-1)/(3-1)
  measure_ar_by_hhsize(hhdata) |>
    expect_equal(data.frame(hhsize = 3L, ar = 2 / 3))
  measure_sar_by_hhsize(hhdata) |>
    expect_equal(data.frame(hhsize = 3L, sar = 0.5))
  measure_cases_by_hhsize(hhdata)$hhsize |> expect_identical(3L)
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
  measure_sar_by_hhsize(as_factor) |> expect_equal(ref_sar)
  measure_ar_by_hhsize(as_factor) |> expect_equal(ref_ar)

  # a character status is decoded, not coerced to NA
  as_char <- as_factor
  as_char$status <- as.character(as_factor$status)
  measure_sar_by_hhsize(as_char) |> expect_equal(ref_sar)
  measure_ar_by_hhsize(as_char) |> expect_equal(ref_ar)

  # an unnamed negative level warns, and still returns the (all-case) measure
  unnamed <- mk_houses$data
  unnamed$status <- ifelse(unnamed$status == 1, "positive", "negative")
  expect_warning(
    all_cases <- measure_ar_by_hhsize(unnamed),
    "escaped label"
  )
  all_cases$ar |> expect_all_equal(1)

  # naming it silences the guard and recovers the coded result
  measure_ar_by_hhsize(unnamed, escaped = "negative") |>
    expect_equal(ref_ar)
})

##########################################
## PREDICTIVE CHECKS : SUMMARY HELPERS  ##
##########################################

test_that("Credible interval level is validated", {
  validate_cri_level(0.95) |> expect_equal(0.95)
  validate_cri_level(0.8) |> expect_equal(0.8)
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
  result$a |> expect_equal(c(3, 4))
  result$b |> expect_equal(c(5, 6))

  # named vectors
  ppc_frame_apply(list(c(x = 1, y = 2), c(x = 3, y = 4)), mean) |>
    expect_equal(c(x = 2, y = 3))
})

test_that("Replicate summaries handle mean, median, NA and non-numeric", {
  # multi-row data.frames required by internal apply() on matrix path
  frame_1 <- data.frame(a = c(1, 2), b = c(NA_real_, NA_real_))
  frame_2 <- data.frame(a = c(3, 4), b = c(NA_real_, NA_real_))

  result <- ppc_frame_apply_summary(list(frame_1, frame_2), "mean")
  result$a |> expect_equal(c(2, 3))
  result$b |> expect_all_na()

  # median path
  frame_3 <- data.frame(a = c(10, 20), b = c(NA_real_, NA_real_))
  ppc_frame_apply_summary(list(frame_1, frame_2, frame_3), "median")$a |>
    expect_equal(c(3, 4))

  # named vector path: non-numeric uniform values
  result <- ppc_frame_apply_summary(list(c(x = 1, y = NA_real_),
                                         c(x = 3, y = NA_real_)), "mean")
  result[["x"]] |> expect_equal(2)
  result[["y"]] |> expect_all_na()

  # invalid summary stat
  ppc_frame_apply_summary(list(frame_1, frame_2), "mode") |>
    expect_error("one between")
})

test_that("Replicate credible intervals bound each other", {
  frames <- list(data.frame(a = c(1, 10)), data.frame(a = c(2, 20)),
                 data.frame(a = c(3, 30)))
  upper <- ppc_frame_apply_cri(frames, cri_level = 0.95, upper = TRUE)
  lower <- ppc_frame_apply_cri(frames, cri_level = 0.95, upper = FALSE)
  (upper$a >= lower$a) |> expect_all_true()

  # named vector path: all-NA returns NA
  ppc_frame_apply_cri(list(c(x = NA_real_), c(x = NA_real_)))[["x"]] |>
    is.na() |>
    expect_true()
})

test_that("Credible intervals are formatted on vectors and frames", {
  # basic formatting
  format_credible_interval(1.23, 0.9, 1.5) |> expect_match("\\[")

  # as_integer
  format_credible_interval(2.1, 1.7, 2.3, as_integer = TRUE,
                           fmt = "%s-%s-%s") |>
    expect_identical("2-2-2")

  # NA mean returns NA_character_
  format_credible_interval(NA, NA, NA) |> expect_identical(NA_character_)

  # equal bounds, and a zero mean for the digit count
  format_credible_interval(2, 2, 2) |> expect_type("character")
  format_credible_interval(0, -0.1, 0.1) |> expect_type("character")

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
  result |> expect_s3_class("data.frame")
  dim(result) |> expect_identical(dim(mean_frame))

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
  names(result) |> expect_contains(c("predicted", "cri_lower", "cri_upper"))
  result$predicted |> expect_s3_class("data.frame")
})

test_that("Houses are counted by household size", {
  hhdata <- data.frame(
    hhid = c(rep("h1", 2), rep("h2", 3), rep("h3", 2), rep("h4", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0, 0, 0)
  )
  hhdata$id <- paste0("m", seq_len(nrow(hhdata)))
  result <- helper_house_counts(hhdata)
  colnames(result) |> expect_contains(c("hhsize", "total", "at_least_1_case"))
  # size 2: h1 (2 cases), h3 (1 case), h4 (0 cases) -> 3 total, 2 affected
  result$total[result$hhsize == 2] |> expect_equal(3)
  result$at_least_1_case[result$hhsize == 2] |> expect_equal(2)
  # size 3: h2 (1 case) -> 1 total, 1 affected
  result$total[result$hhsize == 3] |> expect_equal(1)
  result$at_least_1_case[result$hhsize == 3] |> expect_equal(1)
})

test_that("Tables are filtered by household size band", {
  # one case column per count a household of the largest size can hold
  observed <- data.frame(hhsize = 2:4,
                         `0` = c(1, 2, 3), `1` = c(4, 5, 0), `2` = c(0, 1, 0),
                         `3` = c(0, 1, 1), `4` = c(0, 0, 1),
                         check.names = FALSE)
  scaled <- \(k) transform(observed * k, hhsize = observed$hhsize)
  tables <- list(observed = observed, predicted = scaled(1.1),
                 cri_lower = scaled(0.9), cri_upper = scaled(1.2))

  # no filter
  helper_filter_tables(tables)$observed |> expect_nrow(3L)
  # upper, lower and both bounds
  helper_filter_tables(tables, max_hhsize = 3)$observed$hhsize |>
    expect_identical(2:3)
  helper_filter_tables(tables, min_hhsize = 3)$observed$hhsize |>
    expect_identical(3:4)
  helper_filter_tables(tables, min_hhsize = 3, max_hhsize = 3) |>
    (\(x) x$observed$hhsize)() |>
    expect_identical(3L)
})

test_that("Observations and predictions are interleaved row by row", {
  result <- helper_combine_obs_pred(
    data.frame(hhsize = c(2, 3), `0` = c(1, 2), `1` = c(3, 4),
               check.names = FALSE),
    data.frame(`0` = c(1.1, 2.1), `1` = c(3.1, 4.1), check.names = FALSE)
  )
  colnames(result) |> expect_contains(c("distr", "hhsize"))
  result$distr |> expect_identical(rep(c("observed", "predicted"), 2))
})

test_that("Observed entries outside the interval are flagged", {
  # numeric vectors: middle is below, last is above, NA propagates as FALSE
  list(observed = c(1, 5, 10, NA),
       cri_lower = c(0, 6, 0, 0),
       cri_upper = c(2, 8, 9, 10)) |>
    ppc_observed_outside_cri() |>
    expect_identical(c(FALSE, TRUE, TRUE, FALSE))

  # non-numeric input is treated as not-a-mismatch (skipped)
  list(observed = c("a", "b"),
       cri_lower = c("a", "b"),
       cri_upper = c("a", "b")) |>
    ppc_observed_outside_cri() |>
    expect_identical(c(FALSE, FALSE))

  # data.frame: numeric column flagged, non-numeric column ignored
  result <- list(observed = data.frame(x = c(1, 5), label = c("a", "b")),
                 cri_lower = data.frame(x = c(0, 6), label = c("a", "b")),
                 cri_upper = data.frame(x = c(2, 8), label = c("a", "b"))) |>
    ppc_observed_outside_cri()
  result |> expect_s3_class("data.frame")
  result$x |> expect_identical(c(FALSE, TRUE))
  result$label |> expect_identical(c(FALSE, FALSE))
})

##########################################
## PREDICTIVE CHECKS : END-TO-END       ##
##########################################

test_that("Each kind of predictive check refuses the arguments of the other", {
  # household size bounds apply to the standard measures only
  mk_ppc$model |>
    check_predictive(mk_ppc$results, sar = measure_sar_by_hhsize,
                     max_hhsize = 3L) |>
    expect_error("[Mm]in_hhsize.*max_hhsize.*bound only the standard")
  # a long observation process needs a measure on data_long
  mk_ppc$model |>
    check_predictive(mk_ppc$results, obs_process_long = identity) |>
    expect_error("[Oo]bs_process_long.*has no measure")
})

test_that("Predictive check on custom measures has read-only fields", {
  ppc <- mk_ppc$model |>
    check_predictive(
      mk_ppc$results, sar = measure_sar_by_hhsize,
      thin_ess = FALSE, max_nrep = 10
    ) |>
    suppressWarnings()

  class(ppc)[1] |> expect_identical("HousepiPredictiveCheck")
  names(ppc$sar) |> expect_contains(c("observed", "predicted"))
  print(ppc) |> expect_snapshot()

  # read access returns the run settings
  ppc$nrep |> expect_identical(10L)
  ppc$point_estimate |> expect_identical("median")
  ppc$cri_level |> expect_identical(0.95)
  # the summary settings are read off the global options at construction
  withr::with_options(
    list(housepi.point_estimate = "mean", housepi.cri_level = 0.8), {
      averaged <- mk_ppc$model |>
        check_predictive(
          mk_ppc$results, sar = measure_sar_by_hhsize,
          thin_ess = FALSE, max_nrep = 2
        ) |>
        suppressWarnings()
      averaged$point_estimate |> expect_identical("mean")
      averaged$cri_level |> expect_identical(0.8)
    }
  )
  ppc$thin_ess |> expect_false()
  ppc$is_posterior |> expect_true()
  ppc$model_name |> expect_identical("final size")

  # active fields are read-only
  ppc |>
    expect_read_only_field(
      "point_estimate", "cri_level", "nrep",
      "thin_ess", "is_posterior", "model_name"
    )
})

test_that("Predicted cases by house size works as expected", {
  ppc <- check_predictive(mk_ppc$model, mk_ppc$results) |>
    suppressWarnings()
  class(ppc)[1] |> expect_identical("HousepiEpiPredictiveCheck")
  counts <- helper_house_counts(mk_ppc$model$data_wide)

  # formatted_attack_rate: counts are the raw house totals on the AR grid
  ar <- ppc$formatted_attack_rate
  colnames(ar) |> expect_contains(c("hhsize", "observed", "counts"))
  ar$counts |> expect_identical(counts$total[match(ar$hhsize, counts$hhsize)])

  # formatted_sec_attack_rate: counts are houses with at least one case
  sar <- ppc$formatted_sec_attack_rate
  colnames(sar) |> expect_contains(c("hhsize", "observed"))
  sar$counts |>
    expect_identical(counts$at_least_1_case[match(sar$hhsize, counts$hhsize)])

  # formatted_cases_count
  ppc$formatted_cases_count |> expect_length(3)

  # subclass active fields and inherited dynamic measure fields are read-only
  ppc |>
    expect_read_only_field(
      "formatted_attack_rate", "attack_rate",
      "formatted_sec_attack_rate", "sec_attack_rate",
      "formatted_cases_count", "cases_count"
    )

  # set_max_hhsize and drop
  ppc$set_max_hhsize(3)
  (ppc$formatted_sec_attack_rate$hhsize <= 3) |> expect_all_true()
  ppc$drop_max_hhsize()

  # set_min_hhsize and drop
  ppc$set_min_hhsize(3)
  (ppc$formatted_sec_attack_rate$hhsize >= 3) |> expect_all_true()
  (ppc$formatted_attack_rate$hhsize >= 3) |> expect_all_true()
  ppc$drop_min_hhsize()
})

test_that("Epi PPC honours min/max_hhsize arguments", {
  ppc <- mk_ppc$model |>
    check_predictive(mk_ppc$results, min_hhsize = 2L, max_hhsize = 3L) |>
    suppressWarnings()
  ppc$formatted_attack_rate$hhsize |> expect_in(2:3)
  ppc$formatted_sec_attack_rate$hhsize |> expect_in(2:3)
})

test_that("Epi PPC observed values are stable across calls and status codes", {
  data_wide <- mk_ppc$model$data_wide
  run_on <- function(...) {
    mk_ppc$model |>
      check_predictive(mk_ppc$results, thin_ess = FALSE, max_nrep = 5, ...) |>
      suppressWarnings()
  }
  ppc <- run_on()
  # repeated PPC calls must not mutate the model's data_wide
  mk_ppc$model$data_wide |> expect_identical(data_wide)

  # without an observation process the replicates stay engine-coded
  observed <- measure_ar_by_hhsize(data_wide, mk_ppc$model$escaped)
  ppc$formatted_attack_rate$observed |> expect_equal(observed$ar)
  # with one, they are relabelled before the measure sees them
  run_on(obs_process = \(x) x)$formatted_attack_rate$observed |>
    expect_equal(observed$ar)

  # observed values agree across thread counts
  ppc_2 <- run_on(threads = 2L)
  ppc_2$cases_count$observed |> expect_identical(ppc$cases_count$observed)
  ppc_2$sec_attack_rate$observed |>
    expect_identical(ppc$sec_attack_rate$observed)
  ppc_2$cases_table_by_size$observed |>
    expect_identical(ppc$cases_table_by_size$observed)
})

test_that("Custom PPC warns on observed values outside the interval", {
  # parameter set deliberately so that SAR~1
  param_table <- data.frame(alpha = 50, beta = 50, gamma = 1, hhsize_ref = 4)
  mk_ppc$model |>
    check_predictive(param_table, sar = measure_sar_by_hhsize,
                     thin_ess = FALSE, max_nrep = 20L) |>
    expect_warning(".*outside.*CrI")
})
