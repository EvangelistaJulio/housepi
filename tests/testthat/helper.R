# Helper to test read-only active fields.
expect_read_only_field <- function(obj, ...) {
  for (f in list(...)) {
    (obj[[f]] <- NULL) |> testthat::expect_error(regexp = "read only")
  }
}

# Mask the durations a print reports, which vary from run to run.
helper_scrub_secs <- function(lines) sub("(time[^:]*: ).*$", "\\1<secs>", lines)

# Lines a print emits as messages, stripped of their indentation.
helper_print_lines <- function(expr) {
  trimws(unlist(strsplit(testthat::capture_messages(expr), "\n")))
}

# Give the latents of a model a sink, so inference samples them at all.
# A uniform observation over a window nothing can leave contributes a constant,
# leaving the posterior of what it reads equal to what the model declares.
# `data` must already carry the column the observation reads its datum from.
helper_flat_sink <- function(hh_model, ..., column = "obs") {
  hh_model |>
    add_observation(column, model = "uniform", min = -1e6, max = 1e6, ...)
}

# Case-ascertained household of three cases and four escapees followed for 21
# days, the introduction at day 0, and the transmission parameters it is
# evaluated at.
helper_onset_house <- list(
  data = data.frame(id = sprintf("m%s", seq(7)), hh_id = "h1",
                    status = c(rep(1, 3), rep(0, 4)),
                    infec_time = c(0, 3, 4, rep(NA, 4)),
                    onset_date = c(2, 6, 6, rep(NA, 4)),
                    end_date = 21),
  param = list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
               hh_size_ref = 4, shape = 2, scale = 1.1)
)

# Expect every element to be finite.
expect_all_finite <- function(object) {
  testthat::expect_all_true(is.finite(object))
}

# Expect every element to be missing, or none of them; an empty vector passes.
# The failure names the offending positions.
expect_all_na <- function(object) {
  expect_na_share_(testthat::quasi_label(rlang::enquo(object)), TRUE)
}
expect_no_na <- function(object) {
  expect_na_share_(testthat::quasi_label(rlang::enquo(object)), FALSE)
}
expect_na_share_ <- function(act, missing) {
  off <- which(is.na(act$val) != missing)
  testthat::expect(length(off) == 0, sprintf(
    "Expected %s to be %s.\nOffending position(s): %s", act$lab,
    if (missing) "all missing" else "free of missing values",
    paste(utils::head(off, 10), collapse = ", ")
  ))
  invisible(act$val)
}

# Expect number of row of `object` equal to `n`.
expect_nrow <- function(object, n) {
  testthat::expect_identical(nrow(object), n)
}

# convert draws to data.frame
draws_to_frame <- function(x) {
  as.data.frame(posterior::as_draws_matrix(x))
}
