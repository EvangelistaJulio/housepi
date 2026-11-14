library(testthat)

####  Helpers  ####

# one short inference run, quiet
mk_fit <- local({
  hh_fit <- flu_1918_frederick |>
    new_housepi_model("final-size") |>
    add_parameter(import = "template") |>
    run_inference(steps = 200L, chains = 1L, adapt = 200L, thin = 1L,
                  warn = FALSE)
  list(fit = hh_fit)
})

####  Definitions  ####

test_that("Credible intervals are formatted on vectors and frames", {
  # basic formatting
  format_estimates(1.23, 0.9, 1.5) |> expect_match("\\[")

  # as_integer
  format_estimates(2.1, 1.7, 2.3, as_integer = TRUE,
                   fmt = "%s-%s-%s") |>
    expect_identical("2-2-2")

  # NA mean returns NA_character_
  format_estimates(NA, NA, NA) |> expect_identical(NA_character_)

  # equal bounds
  format_estimates(2, 2, 2) |> expect_type("character")

  # two significant digits of the estimate, one of the width, none above units
  format_estimates(c(1.97, 9.9, 1234, 0.52013),
                   c(0.673, 2.6, 1180, 0.51904),
                   c(3.24, 23.7, 1290, 0.52097)) |>
    expect_identical(c("2.0 [0.7, 3.2]", "9.9 [2.6, 23.7]",
                       "1234 [1180, 1290]", "0.520 [0.519, 0.521]"))

  # interval covering zero -> two significant digits of the width
  format_estimates(c(0.0012, 0), c(-0.5, 0), c(0.5, 0.3)) |>
    expect_identical(c("0.0 [-0.5, 0.5]", "0.00 [0.00, 0.30]"))

  # inconsistent interval warns
  format_estimates(5, 10, 20) |>
    expect_warning("inconsistent")

  # mismatched lengths error
  format_estimates(c(1, 2), c(1), c(2)) |>
    expect_error("same length")

  # wrong str_format error
  format_estimates(1, 0, 2, fmt = "%s only one") |>
    expect_error("3.*placeholder")

  # data.frame input
  mean_frame <- data.frame(a = c(1, 2), b = c(3, 4))
  upper_frame <- data.frame(a = c(2, 3), b = c(4, 5))
  result <- format_estimates(
    mean_frame, data.frame(a = c(0, 1), b = c(2, 3)), upper_frame
  )
  result |> expect_s3_class("data.frame")
  result |> expect_dim(dim(mean_frame))

  # mismatched data.frame dimensions
  format_estimates(mean_frame, data.frame(a = c(0)), upper_frame) |>
    expect_error("same dimensions")
})

test_that("Point estimates alone are formatted to two significant digits", {
  format_estimates(c(0.0123, 45.6, 1234, NA)) |>
    expect_identical(c("0.012", "46", "1234", NA))

  # as_integer
  format_estimates(2.6, as_integer = TRUE) |> expect_identical("3")

  # a matrix keeps its shape
  format_estimates(matrix(c(1.23, 45.6), 1)) |>
    expect_dim(c(1L, 2L))

  # an interval placeholder without bounds errors
  format_estimates(1.5, fmt = "%s [%s, %s]") |>
    expect_error("exactly 1.*placeholder")

  # one bound alone errors
  format_estimates(1.5, cri_lower = 1) |> expect_error("pairs")
})

test_that("A posterior summary is formatted from the global summary", {
  old <- set_housepi_summary(point_estimate = "median", cri_level = 0.95)
  withr::defer(do.call(set_housepi_summary, old))
  summary_95 <- summarise_draws(mk_fit$fit)
  estimates <- format_estimates(summary_95)
  # every column kept, `formatted` appended
  estimates |> expect_named(c(names(summary_95), "formatted"))
  estimates[names(summary_95)] |> expect_identical(summary_95)
  estimates$formatted |>
    expect_identical(format_estimates(summary_95$median, summary_95$q2.5,
                                      summary_95$q97.5))

  # columns follow the global point estimate and level
  set_housepi_summary(point_estimate = "mean", cri_level = 0.8)
  summary_80 <- summarise_draws(mk_fit$fit)
  format_estimates(summary_80, fmt = "%s (%s, %s)")$formatted |>
    expect_identical(format_estimates(summary_80$mean, summary_80$q10,
                                      summary_80$q90, fmt = "%s (%s, %s)"))

  # a summary taken at another level lacks the columns
  set_housepi_summary(cri_level = 0.95)
  format_estimates(summary_80) |> expect_error("no column.*q2.5")
})

test_that("A table with credible interval columns is formatted at its bounds", {
  # shaped as a predictive check measure
  hh_data <- data.frame(hh_size = 2:3, observed = c(0.4, 0.3),
                        predicted = c(0.42, 0.28), cri_lower = c(0.3, 0.2),
                        cri_upper = c(0.55, 0.37), outside = FALSE)
  formatted <- format_estimates(hh_data, as_integer = TRUE)
  formatted |> expect_named(c(names(hh_data), "formatted"))
  formatted[names(hh_data)] |> expect_identical(hh_data)
  formatted$formatted |>
    expect_identical(format_estimates(hh_data$predicted, hh_data$cri_lower,
                                      hh_data$cri_upper, as_integer = TRUE))

  # point estimate named as in the performance checks
  hh_data_2 <- data.frame(variable = "rho", estimate = 1.97,
                          cri_lower = 0.673, cri_upper = 3.24)
  format_estimates(hh_data_2)$formatted |> expect_identical("2.0 [0.7, 3.2]")

  # two point estimate columns, or none, error
  cbind(hh_data_2, predicted = 2) |>
    format_estimates() |>
    expect_error("one point estimate column")
  hh_data_2[c("variable", "cri_lower", "cri_upper")] |>
    format_estimates() |>
    expect_error("one point estimate column")

  # one bound alone errors
  hh_data_2[c("estimate", "cri_lower")] |>
    format_estimates() |>
    expect_error("no column.*cri_upper")
})
