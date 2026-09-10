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
  add_observation(hhmodel, column, model = "uniform",
                  min = -1e6, max = 1e6, ...)
}
