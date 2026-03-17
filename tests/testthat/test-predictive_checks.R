library(testthat)

test_that("measure_sar_by_house_size computes SAR correctly", {
  # household of size 2: 2 cases => SAR = (2-1)/(2-1) = 1
  # household of size 3: 1 case  => excluded (0 cases after filtering)
  # household of size 2: 1 case  => SAR = (1-1)/(2-1) = 0
  data_wide <- data.frame(
    id = paste0("p", 1:7),
    house_id = c("hh1", "hh1", "hh2", "hh2", "hh2", "hh3", "hh3"),
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
  data_wide <- data.frame(
    id = paste0("p", 1:5),
    house_id = c("hh1", "hh1", "hh2", "hh2", "hh2"),
    status = c(1, 1, 1, 0, 0)
  )

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
  expect_equal(housepi:::ppc_set_cri_level(0.95), 0.95)
  expect_equal(housepi:::ppc_set_cri_level(0.8), 0.8)
  # boundary and invalid values
  expect_error(housepi:::ppc_set_cri_level(0.5), regexp = "between 0.5 and 1")
  expect_error(housepi:::ppc_set_cri_level(1.0), regexp = "between 0.5 and 1")
  expect_error(housepi:::ppc_set_cri_level(0.3), regexp = "between 0.5 and 1")
  expect_error(housepi:::ppc_set_cri_level(NA), regexp = "cannot be NA")
})

test_that("ppc_frame_apply computes element-wise operations", {
  # list of data.frames
  df1 <- data.frame(a = c(1, 2), b = c(3, 4))
  df2 <- data.frame(a = c(5, 6), b = c(7, 8))
  result <- housepi:::ppc_frame_apply(list(df1, df2), mean)
  expect_equal(result$a, c(3, 4))
  expect_equal(result$b, c(5, 6))

  # named vectors
  v1 <- c(x = 1, y = 2)
  v2 <- c(x = 3, y = 4)
  result <- housepi:::ppc_frame_apply(list(v1, v2), mean)
  expect_equal(result, c(x = 2, y = 3))
})

test_that("ppc_frame_apply_mean handles NA and non-numeric", {
  # multi-row data.frames required by internal apply() on matrix path
  df1 <- data.frame(a = c(1, 2), b = c(NA_real_, NA_real_))
  df2 <- data.frame(a = c(3, 4), b = c(NA_real_, NA_real_))

  result <- housepi:::ppc_frame_apply_mean(list(df1, df2))
  expect_equal(result$a, c(2, 3))
  expect_true(all(is.na(result$b)))

  # named vector path: non-numeric uniform values
  v1 <- c(x = 1, y = NA_real_)
  v2 <- c(x = 3, y = NA_real_)
  result <- housepi:::ppc_frame_apply_mean(list(v1, v2))
  expect_equal(result[["x"]], 2)
  expect_true(is.na(result[["y"]]))
})

test_that("ppc_frame_apply_cri computes quantiles", {
  # multi-row data.frames
  df1 <- data.frame(a = c(1, 10))
  df2 <- data.frame(a = c(2, 20))
  df3 <- data.frame(a = c(3, 30))

  upper <- housepi:::ppc_frame_apply_cri(list(df1, df2, df3),
                                         cri_level = 0.95, upper = TRUE)
  lower <- housepi:::ppc_frame_apply_cri(list(df1, df2, df3),
                                         cri_level = 0.95, upper = FALSE)

  expect_true(all(upper$a >= lower$a))

  # named vector path: all-NA returns NA
  v1 <- c(x = NA_real_)
  v2 <- c(x = NA_real_)
  result <- housepi:::ppc_frame_apply_cri(list(v1, v2))
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
})
