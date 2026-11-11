library(testthat)

####  Helpers  ####

# templated chain-binomial model and one short inference run on it
mk_ppc <- local({
  hh_model <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template")
  hh_fit <- hh_model |>
    run_inference(steps = 100L, chains = 1L, adapt = 0L, thin = 1L) |>
    suppressWarnings()
  list(model = hh_model, fit = hh_fit)
})

# three households of sizes two, three and two
mk_houses <- local({
  hh_data <- data.frame(
    id = paste0("m", 1:7),
    hh_id = c(rep("h1", 2), rep("h2", 3), rep("h3", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0)
  )
  list(data = hh_data)
})

####  Measures  ####

test_that("Measures require the columns they aggregate over", {
  mk_houses$data[c("id", "hh_id")] |>
    measure_sar_by_hh_size() |>
    expect_error("status")

  # without a household identifier every row would count as its own household
  mk_houses$data[c("id", "status")] |>
    measure_sar_by_hh_size() |>
    expect_error("hh_id")
})

test_that("Attack rate is averaged per household size", {
  # empty call: columns c(hh_size, ar)
  measure_ar_by_hh_size(NULL) |> expect_dim(c(0L, 2L))

  result <- measure_ar_by_hh_size(mk_houses$data)
  result |> expect_named(c("hh_size", "ar"), ignore.order = TRUE)
  # size 2: mean of AR 2/2 (h1) and 1/2 (h3) = 0.75
  result[result$hh_size == 2, "ar"] |> expect_equal(0.75)
  # size 3: AR = 1/3
  result[result$hh_size == 3, "ar"] |> expect_equal(1 / 3)
})

test_that("Secondary attack rate is averaged per household size", {
  # empty call: columns c(hh_size, sar)
  measure_sar_by_hh_size(NULL) |> expect_dim(c(0L, 2L))

  result <- measure_sar_by_hh_size(mk_houses$data)
  result |> expect_named(c("hh_size", "sar"), ignore.order = TRUE)
  # size 2: mean of SAR 1.0 (h1) and 0.0 (h3) = 0.5
  result[result$hh_size == 2, "sar"] |> expect_equal(0.5)
  # size 3: SAR = (1-1)/(3-1) = 0
  result[result$hh_size == 3, "sar"] |> expect_equal(0)
})

test_that("Cases are tabulated on the complete household size grid", {
  # empty call: columns c(hh_size, 0)
  measure_cases_by_hh_size(NULL) |> expect_dim(c(0L, 2L))

  result <- measure_cases_by_hh_size(mk_houses$data)
  colnames(result) |> expect_contains("hh_size")
  result$hh_size |> expect_setequal(c(2, 3))
})

test_that("Measures keep households that all share one size", {
  hh_data <- data.frame(id = paste0("m", 1:3), hh_id = "h1",
                        status = c(1, 1, 0))
  # AR = 2/3, SAR = (2-1)/(3-1)
  measure_ar_by_hh_size(hh_data) |>
    expect_equal(data.frame(hh_size = 3L, ar = 2 / 3))
  measure_sar_by_hh_size(hh_data) |>
    expect_equal(data.frame(hh_size = 3L, sar = 0.5))
  measure_cases_by_hh_size(hh_data)$hh_size |> expect_identical(3L)
})

test_that("Household size geometric mean is weighted by secondary cases", {
  # no household with a secondary case
  measure_geom_mean_hh_size(NULL) |>
    expect_equal(NA_real_)

  # only h1 (size 2, 2 cases) carries weight
  measure_geom_mean_hh_size(mk_houses$data) |>
    expect_equal(2)

  # h4 (size 4, 3 cases) weighs 2: exp((log 2 + 2 log 4) / 3) = 2^(5/3)
  mk_houses$data |>
    rbind(data.frame(
      id = paste0("m", 8:11), hh_id = "h4", status = c(1, 1, 1, 0)
    )) |>
    measure_geom_mean_hh_size() |>
    expect_equal(2^(5 / 3))
})

test_that("Measures decode a labelled status column", {
  ref_ar <- measure_ar_by_hh_size(mk_houses$data)
  ref_sar <- measure_sar_by_hh_size(mk_houses$data)

  # a factor status is decoded through its levels, not through level indices
  as_factor <- mk_houses$data
  as_factor$status = factor(ifelse(as_factor$status == 1,
                                   "symptomatic", "escaped"))
  measure_sar_by_hh_size(as_factor) |> expect_equal(ref_sar)
  measure_ar_by_hh_size(as_factor) |> expect_equal(ref_ar)

  # a character status is decoded, not coerced to NA
  as_char <- as_factor
  as_char$status <- as.character(as_factor$status)
  measure_sar_by_hh_size(as_char) |> expect_equal(ref_sar)
  measure_ar_by_hh_size(as_char) |> expect_equal(ref_ar)

  # an unnamed negative level warns, and still returns the (all-case) measure
  unnamed <- mk_houses$data
  unnamed$status <- ifelse(unnamed$status == 1, "positive", "negative")
  expect_warning(
    all_cases <- measure_ar_by_hh_size(unnamed),
    "escaped label"
  )
  all_cases$ar |> expect_all_equal(1)

  # naming it silences the guard and recovers the coded result
  measure_ar_by_hh_size(unnamed, escaped = "negative") |>
    expect_equal(ref_ar)
})

####  Summary helpers  ####

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
    expect_error("one of")
})

test_that("Replicate credible intervals bound each other", {
  frames <- list(data.frame(a = c(1, 10)), data.frame(a = c(2, 20)),
                 data.frame(a = c(3, 30)))
  upper <- ppc_frame_apply_cri(frames, cri_level = 0.95, upper = TRUE)
  lower <- ppc_frame_apply_cri(frames, cri_level = 0.95, upper = FALSE)
  upper$a |> expect_all_gte(lower$a)

  # named vector path: all-NA returns NA
  ppc_frame_apply_cri(list(c(x = NA_real_), c(x = NA_real_)))[["x"]] |>
    expect_all_na()
})

test_that("Predicted hh_fit carry the point estimate and its interval", {
  result <- ppc_extract_predicted_results(
    list(data.frame(a = c(1, 2)), data.frame(a = c(3, 4)),
         data.frame(a = c(5, 6))),
    cri = 0.95
  )
  names(result) |> expect_contains(c("predicted", "cri_lower", "cri_upper"))
  result$predicted |> expect_s3_class("data.frame")
})

test_that("Houses are counted by household size", {
  hh_data <- data.frame(
    hh_id = c(rep("h1", 2), rep("h2", 3), rep("h3", 2), rep("h4", 2)),
    status = c(1, 1, 1, 0, 0, 1, 0, 0, 0)
  )
  hh_data$id <- paste0("m", seq_len(nrow(hh_data)))
  result <- helper_house_counts(hh_data)
  colnames(result) |> expect_contains(c("hh_size", "total", "at_least_1_case"))
  # size 2: h1 (2 cases), h3 (1 case), h4 (0 cases) -> 3 total, 2 affected
  result$total[result$hh_size == 2] |> expect_equal(3)
  result$at_least_1_case[result$hh_size == 2] |> expect_equal(2)
  # size 3: h2 (1 case) -> 1 total, 1 affected
  result$total[result$hh_size == 3] |> expect_equal(1)
  result$at_least_1_case[result$hh_size == 3] |> expect_equal(1)
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

####  End-to-end  ####

test_that("Household size bounds need the standard measures", {
  mk_ppc$model |>
    check_predictive(mk_ppc$fit, sar = measure_sar_by_hh_size,
                     epi = FALSE, max_hh_size = 3L) |>
    expect_error("[Mm]in_hh_size.*max_hh_size.*bound.*only.*standard")

  # next to custom measures, the bounds restrict the standard outputs
  ppc <- mk_ppc$model |>
    check_predictive(mk_ppc$fit, overall = \(x) mean(x$status > 0),
                     max_hh_size = 3L, thin_ess = FALSE, nrep = 5) |>
    suppressWarnings()
  ppc$format_attack_rate()$hh_size |> expect_all_lte(3)
})

test_that("Custom measures are checked next to the standard ones", {
  ppc <- mk_ppc$model |>
    check_predictive(mk_ppc$fit, overall = \(x) mean(x$status > 0),
                     thin_ess = FALSE, nrep = 5) |>
    suppressWarnings()
  class(ppc)[1] |> expect_identical("HousepiEpiPredictiveCheck")
  names(ppc$overall) |> expect_contains(c("observed", "predicted"))
  names(ppc$attack_rate) |> expect_contains(c("observed", "predicted"))
  colnames(ppc$format_attack_rate()) |> expect_contains("observed")
  ppc |> expect_read_only_field("overall", "attack_rate")
  # the print links standard measures to their help page, custom ones stay plain
  withr::local_options(cli.hyperlink = TRUE, cli.hyperlink_run = TRUE)
  lines <- helper_print_lines(print(ppc))
  lines |>
    expect_match("x-r-run:utils::help\\(\"measure_ar_by_hh_size\"", all = FALSE)
  lines |> expect_match(", and overall$", all = FALSE)
})

test_that("Custom measures return values of one shape that summarise", {
  mk_ppc$model |>
    check_predictive(mk_ppc$fit, empty = \(x) numeric(), epi = FALSE,
                     thin_ess = FALSE, nrep = 5) |>
    expect_error("[Ee]mpty observation on data")
  # one value per case: replicates differ in length
  mk_ppc$model |>
    check_predictive(mk_ppc$fit, ragged = \(x) seq_len(sum(x$status > 0)),
                     epi = FALSE, thin_ess = FALSE, nrep = 5) |>
    expect_error("[Nn]ot all computed measures share identical dimensions") |>
    suppressWarnings()
  # a label that differs across replicates has no median
  mk_ppc$model |>
    check_predictive(mk_ppc$fit, label = \(x) as.character(sum(x$status > 0)),
                     epi = FALSE, thin_ess = FALSE, nrep = 5) |>
    expect_error("[Uu]ndefined .median.*non-numeric, non-unique") |>
    suppressWarnings()
})

test_that("Measures with an escaped argument receive the model's label", {
  hh_data <- flu_1918_frederick
  hh_data$status <- ifelse(hh_data$status == 0, "negative", "positive")
  hh_model <- hh_data |>
    new_housepi_model("final_size", escaped = "negative") |>
    add_parameter(import = "template")
  param_table <- data.frame(alpha = 0.1, beta = 0.3, gamma = 1,
                            hh_size_ref = 4L)
  run <- function(...) {
    hh_model |> check_predictive(param_table, ..., epi = FALSE, nrep = 5)
  }
  # the measure's own default yields to the model's label
  infected <- function(data, escaped = "escaped") mean(data$status != escaped)
  ppc <- run(sar = measure_sar_by_hh_size, overall = infected) |>
    expect_no_warning(message = "escaped label") |>
    suppressWarnings()
  ppc$sar$observed |>
    expect_equal(measure_sar_by_hh_size(hh_data, "negative"))
  ppc$overall$observed |> expect_equal(mean(hh_data$status == "positive"))
  ppc$overall$predicted |> expect_lt(1)
  # a measure without the argument keeps its own label
  run(sar = \(x) measure_sar_by_hh_size(x)) |>
    expect_warning("escaped label") |>
    suppressWarnings()
})

test_that("Standard measures can be dropped", {
  ppc <- mk_ppc$model |>
    check_predictive(mk_ppc$fit, overall = \(x) mean(x$status > 0),
                     epi = FALSE, thin_ess = FALSE, nrep = 5) |>
    suppressWarnings()
  class(ppc)[1] |> expect_identical("HousepiPredictiveCheck")
  ppc$format_attack_rate |> expect_null()
  ppc$attack_rate |> expect_null()
  ppc$overall$observed |> expect_type("double")

  mk_ppc$model |>
    check_predictive(mk_ppc$fit, epi = FALSE) |>
    expect_error("[Nn]o measure to check")
  mk_ppc$model |>
    check_predictive(mk_ppc$fit, epi = NA) |>
    expect_error("epi.*TRUE.*FALSE")
})

test_that("Measure names are unique, given and free of the members", {
  overall <- \(x) mean(x$status > 0)
  run <- \(...) mk_ppc$model |> check_predictive(mk_ppc$fit, ...)
  run(overall) |> expect_error("must be named")
  run(overall = overall, overall = overall) |>
    expect_error("[Rr]epeated measure.*overall")
  # standard names are free once the standard measures are dropped
  run(attack_rate = overall) |> expect_error("attack_rate.*used.*standard")
  ppc <- run(attack_rate = overall, epi = FALSE) |> suppressWarnings()
  ppc$attack_rate$observed |> expect_type("double")
  # a measure cannot take the name of a method or field of the check
  run(print = overall, epi = FALSE) |>  expect_error("reserved")
  run(format_attack_rate = overall) |> expect_error("reserved")
})

test_that("Long measures need a long dataset reporting process", {
  mk_ppc$model |>
    check_predictive(mk_ppc$fit, long_measures = list(n = nrow)) |>
    expect_error("[Ll]ong dataset reporting process")
})

test_that("Long measures are checked alongside those on the dataset", {
  cases_long <- \(x) data.frame(id = x$id, case = as.integer(x$status > 0))
  hh_model <- mk_ppc$model |>
    add_long_dataset(cases_long(get_dataset(mk_ppc$model))) |>
    add_reporting_process(long_dataset = cases_long)
  run <- function(...) {
    check_predictive(hh_model, mk_ppc$fit, ..., nrep = 3, thin_ess = FALSE,
                     progress = FALSE)
  }

  # one count read from either frame of the same replicates
  ppc <- run(long_measures = list(cases = \(x) sum(x$case)),
             wide_cases = \(x) sum(x$status > 0)) |>
    suppressWarnings()
  ppc$cases |> expect_equal(ppc$wide_cases)
  ppc$nrep |> expect_identical(3L)
  # measures on long_dataset alone, the standard ones dropped
  run(long_measures = list(cases = \(x) sum(x$case)), epi = FALSE) |>
    suppressWarnings() |>
    expect_s3_class("HousepiPredictiveCheck")
  # named as those in `...`, under names free across both
  run(long_measures = \(x) sum(x$case)) |>
    expect_error("long_measures.*must be named")
  run(long_measures = list(cases = nrow), cases = nrow) |>
    expect_error("[Rr]epeated measure.*cases")
  run(long_measures = list(attack_rate = nrow)) |>
    expect_error("attack_rate.*standard")
})

test_that("A predictive check takes time_step as a run setting", {
  run <- function(...) {
    check_predictive(mk_ppc$model, mk_ppc$fit, ..., nrep = 2, progress = FALSE)
  }
  # the step reaches the engine, as in run_simulation()
  run(time_step = -1) |> expect_error("time_step. must be strictly positive")
  ppc <- run(time_step = 0.5) |> suppressWarnings()
  ppc$time_step |> expect_null()
  # never a measure
  run(time_step = 0.5, epi = FALSE) |> expect_error("[Nn]o measure")
})

test_that("Predictive check on custom measures has read-only fields", {
  ppc <- mk_ppc$model |>
    check_predictive(
      mk_ppc$fit, sar = measure_sar_by_hh_size,
      epi = FALSE, thin_ess = FALSE, nrep = 10
    ) |>
    suppressWarnings()

  class(ppc)[1] |> expect_identical("HousepiPredictiveCheck")
  names(ppc$sar) |> expect_contains(c("observed", "predicted"))

  # read access returns the run settings
  ppc$nrep |> expect_identical(10L)
  ppc$point_estimate |> expect_identical("median")
  ppc$cri_level |> expect_identical(0.9)
  # the summary settings are read off the global options at construction
  withr::with_options(
    list(housepi.point_estimate = "mean", housepi.cri_level = 0.8), {
      averaged <- mk_ppc$model |>
        check_predictive(
          mk_ppc$fit, sar = measure_sar_by_hh_size,
          epi = FALSE, thin_ess = FALSE, nrep = 2
        ) |>
        suppressWarnings()
      averaged$point_estimate |> expect_identical("mean")
      averaged$cri_level |> expect_identical(0.8)
    }
  )
  ppc$thin_ess |> expect_false()
  ppc$source |> expect_identical("posterior")
  ppc$model_name |> expect_identical("final size")

  # active fields are read-only
  ppc |>
    expect_read_only_field(
      "point_estimate", "cri_level", "nrep",
      "thin_ess", "source", "model_name", "warnings"
    )
})

test_that("Observed values outside the interval are recorded as warnings", {
  hh_data <- get_dataset(mk_ppc$model)
  # only the observed data equal themselves, a constant always matches
  same <- \(x) as.numeric(identical(x, hh_data))
  run <- function(...) {
    mk_ppc$model |>
      check_predictive(mk_ppc$fit, ..., const = \(x) 1, epi = FALSE,
                       thin_ess = FALSE, nrep = 2, progress = FALSE)
  }
  expect_warning(ppc <- run(same = same), "outside.*CrI for measure same")
  ppc$warnings$message |> expect_match("measure same")
  ppc$warnings$count |> expect_identical(1L)
  helper_print_lines(print(ppc)) |>
    expect_match("^! Warning: 1 issue flagged, see .\\$warnings", all = FALSE)
  print(ppc) |> expect_snapshot(transform = helper_scrub_secs)

  # a check flagging nothing has no rows and no warning line
  ppc_2 <- run()
  ppc_2$warnings |> expect_nrow(0L)
  helper_print_lines(print(ppc_2)) |> expect_no_match("Warning")
})

test_that("The check names the source of the values it simulates at", {
  mk_ppc$model |> check_predictive() |> expect_error("[Pp]rior.*nrep")

  ppc <- mk_ppc$model |>
    check_predictive(nrep = 5, progress = FALSE) |>
    suppressWarnings()
  ppc$source |> expect_identical("prior")
  ppc$nrep |> expect_identical(5L)
  # prior draws are independent: no ESS thinning
  ppc$thin_ess |> expect_false()
  print(ppc) |> expect_message("Prior predictive check")

  # a plain table of values is neither prior nor posterior
  param_table <- data.frame(alpha = 1, beta = 1, gamma = 1, hh_size_ref = 4)
  ppc_2 <- mk_ppc$model |>
    check_predictive(param_table, progress = FALSE) |>
    suppressWarnings()
  ppc_2$source |> expect_identical("table")
  ppc_2$thin_ess |> expect_false()

  # draws taken out of a run stay a posterior
  ppc_3 <- mk_ppc$model |>
    check_predictive(get_parameter_draws(mk_ppc$fit), nrep = 2,
                     progress = FALSE) |>
    suppressWarnings()
  ppc_3$source |> expect_identical("posterior")
  ppc_3$thin_ess |> expect_true()
})

test_that("Draws of a run on other data are given values, not a posterior", {
  # households resampled: same definitions, other data
  hh_model <- get_dataset(mk_ppc$model) |>
    sample_house_structure(houses = 20) |>
    add_dataset(mk_ppc$model, data = _)
  ppc <- hh_model |>
    check_predictive(mk_ppc$fit, sar = measure_sar_by_hh_size,
                     epi = FALSE, nrep = 2, progress = FALSE) |>
    suppressWarnings()
  ppc$source |> expect_identical("table")
  # the draws of a run are thinned to their ESS all the same
  ppc$thin_ess |> expect_true()
  print(ppc) |> expect_message("Predictive check for model")
})

test_that("An epidemic predictive check formats counts per household size", {
  ppc <- check_predictive(mk_ppc$model, mk_ppc$fit) |>
    suppressWarnings()
  class(ppc)[1] |> expect_identical("HousepiEpiPredictiveCheck")
  counts <- helper_house_counts(get_dataset(mk_ppc$model))

  # formatted_attack_rate: counts are the raw house totals on the AR grid
  ar <- ppc$format_attack_rate()
  colnames(ar) |> expect_contains(c("hh_size", "observed", "counts"))
  ar$counts |> expect_identical(counts$total[match(ar$hh_size, counts$hh_size)])

  # formatted_sec_attack_rate: counts are houses with at least one case
  sar <- ppc$format_sec_attack_rate()
  colnames(sar) |> expect_contains(c("hh_size", "observed"))
  sar$counts |>
    expect_identical(counts$at_least_1_case[match(sar$hh_size, counts$hh_size)])

  # formatted_cases_distr: one row per household size and number of cases,
  # observed households of each size add up to its total
  cases <- ppc$format_cases_distr()
  colnames(cases) |> expect_identical(
    c("hh_size", "cases", "observed", "predicted", "cri_lower", "cri_upper",
      "outside")
  )
  cases$cases |> expect_all_lte(cases$hh_size)
  as.vector(tapply(cases$observed, cases$hh_size, sum)) |>
    expect_identical(counts$total[match(unique(cases$hh_size), counts$hh_size)])

  # rates share the columns, `outside` flags observed values beyond the CrI
  colnames(ar) |> expect_identical(
    c("hh_size", "observed", "predicted", "cri_lower", "cri_upper", "outside",
      "counts")
  )
  colnames(sar) |> expect_identical(colnames(ar))
  for (x in list(ar, sar, cases)) {
    x$outside |>
      expect_identical(x$observed < x$cri_lower | x$observed > x$cri_upper)
  }

  # subclass active fields and inherited dynamic measure fields are read-only
  ppc |>
    expect_read_only_field(
      "attack_rate", "sec_attack_rate", "cases_distr"
    )

  # a bound set to NULL or NA is removed
  hh_size <- ppc$format_sec_attack_rate()$hh_size
  ppc$set_max_hh_size(3)
  ppc$format_sec_attack_rate()$hh_size |> expect_all_lte(3)
  ppc$set_max_hh_size(NULL)
  ppc$format_sec_attack_rate()$hh_size |> expect_identical(hh_size)
  ppc$set_min_hh_size(3)
  ppc$format_sec_attack_rate()$hh_size |> expect_all_gte(3)
  ppc$format_attack_rate()$hh_size |> expect_all_gte(3)
  ppc$set_min_hh_size(NA)
  ppc$format_sec_attack_rate()$hh_size |> expect_identical(hh_size)

  # a clone formats with its own bounds
  ppc_2 <- ppc$clone()$set_max_hh_size(3)
  ppc_2$format_sec_attack_rate()$hh_size |> expect_all_lte(3)
  ppc$format_sec_attack_rate()$hh_size |> expect_identical(hh_size)
})

test_that("A case-ascertained check leaves the attack rate out", {
  # households enter the study through a case
  hh_model <- sarscov2_2021_icofs |>
    dplyr::filter(hh_id %in% utils::head(unique(hh_id), 20)) |>
    dplyr::mutate(end_date = 21) |>
    new_housepi_model("case-ascertained") |>
    add_infectivity_profile(import = "template")
  param_table <- data.frame(alpha = 1e-3, beta = 0.5, gamma = 1,
                            hh_size_ref = 4)
  ppc <- hh_model |>
    check_predictive(param_table, nrep = 2L, progress = FALSE) |>
    suppressWarnings()
  ppc$attack_rate |> expect_null()
  ppc$format_attack_rate |> expect_null()
  ppc$format_sec_attack_rate() |> expect_s3_class("data.frame")
  # every household holds a case: no row with zero cases
  ppc$format_cases_distr()$cases |> expect_all_gt(0L)

  # the name stays reserved
  hh_model |>
    check_predictive(param_table, attack_rate = \(x) 1, nrep = 2L) |>
    expect_error("[Aa]ttack_rate.*already used")
})

test_that("Epi PPC honours min/max_hh_size arguments", {
  ppc <- mk_ppc$model |>
    check_predictive(mk_ppc$fit, min_hh_size = 2L, max_hh_size = 3L) |>
    suppressWarnings()
  ppc$format_attack_rate()$hh_size |> expect_in(2:3)
  ppc$format_sec_attack_rate()$hh_size |> expect_in(2:3)
  ppc$format_cases_distr()$hh_size |> expect_in(2:3)
})

test_that("Epi PPC observed values are stable across calls and status codes", {
  dataset <- get_dataset(mk_ppc$model)
  run_on <- function(hh_model = mk_ppc$model, ...) {
    hh_model |>
      check_predictive(mk_ppc$fit, thin_ess = FALSE, nrep = 5, ...) |>
      suppressWarnings()
  }
  ppc <- run_on()
  # repeated PPC calls must not mutate the model's dataset
  get_dataset(mk_ppc$model) |> expect_identical(dataset)

  # without a reporting process the replicates stay engine-coded
  observed <- measure_ar_by_hh_size(dataset, mk_ppc$model$escaped)
  ppc$format_attack_rate()$observed |> expect_equal(observed$ar)
  # with one, they are relabelled before the measure sees them
  ppc_labelled <- mk_ppc$model |>
    add_reporting_process(dataset = identity) |>
    run_on()
  ppc_labelled$format_attack_rate()$observed |> expect_equal(observed$ar)

  # observed values agree across thread counts
  ppc_2 <- run_on(threads = 2L)
  ppc_2$cases_distr$observed |> expect_identical(ppc$cases_distr$observed)
  ppc_2$sec_attack_rate$observed |>
    expect_identical(ppc$sec_attack_rate$observed)
})

test_that("Custom PPC warns on observed values outside the interval", {
  # parameter set deliberately so that SAR~1
  param_table <- data.frame(alpha = 50, beta = 50, gamma = 1, hh_size_ref = 4)
  mk_ppc$model |>
    check_predictive(param_table, sar = measure_sar_by_hh_size,
                     epi = FALSE, thin_ess = FALSE, nrep = 20L) |>
    expect_warning(".*outside.*CrI")
})
