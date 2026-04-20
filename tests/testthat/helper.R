# Helper to load definition lists.
helper_test_aux_load <- function(self, ..., append = FALSE) {
  self$load_quosures(rlang::enquos(...), append = append)
  invisible(self)
}
