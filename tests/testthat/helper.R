# Helper to test read-only active fields.
helper_test_read_only <- function(obj, ...) {
  for (f in list(...)) {
    testthat::expect_error(obj[[f]] <- NULL, regexp = "read only")
  }
}

# Give the latents of a model a sink, so inference samples them at all.
# A uniform observation over a window nothing can leave contributes a constant,
# leaving the posterior of what it reads equal to what the model declares.
# `data` must already carry the column the observation reads its datum from.
helper_flat_sink <- function(hhmodel, ..., column = "obs") {
  hhmodel |>
    add_observation(column, model = "uniform", min = -1e6, max = 1e6, ...)
}

# Case-ascertained household of three cases and four escapees followed for 21
# days, the introduction at day 0, and the transmission parameters it is
# evaluated at.
helper_onset_house <- list(
  data = data.frame(id = sprintf("m%s", seq(7)), hhid = "h1",
                    status = c(rep(1, 3), rep(0, 4)),
                    infec_time = c(0, 3, 4, rep(NA, 4)),
                    onset_date = c(2, 6, 6, rep(NA, 4)),
                    start_date = 0, end_date = 21),
  param = list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
               hhsize_ref = 4, shape = 2, scale = 1.1)
)

# convert draws to data.frame
draws_to_frame <- function(x) {
  as.data.frame(posterior::as_draws_matrix(x))
}
