# Helper to load definition lists.
helper_test_aux_load <- function(self, ..., append = FALSE) {
  self$add_quo_list(rlang::enquos(...), append = append)
  invisible(self)
}

# Helper to test read-only active fields.
helper_test_read_only <- function(obj, ...) {
  for (f in list(...)) {
    testthat::expect_error(obj[[f]] <- NULL, regexp = "read only")
  }
}
