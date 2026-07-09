# Helper to load definition lists.
helper_test_aux_load <- function(self, ..., append = TRUE,
                                 .kind = NULL) {
  # only the augmented loader accepts '.obs' / '.house'; keep the parameter
  # loader signature untouched otherwise
  if (is.null(.kind)) {
    self$add_quo_list(rlang::enquos(...), append = append)
  } else {
    self$add_quo_list(rlang::enquos(...), append = append, .kind = .kind)
  }
  invisible(self)
}

# Helper to test read-only active fields.
helper_test_read_only <- function(obj, ...) {
  for (f in list(...)) {
    testthat::expect_error(obj[[f]] <- NULL, regexp = "read only")
  }
}
