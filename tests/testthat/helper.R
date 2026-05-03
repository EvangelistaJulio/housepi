# Helper to load definition lists.
helper_test_aux_load <- function(self, ..., append = FALSE) {
  self$add_quo_list(rlang::enquos(...), append = append)
  invisible(self)
}
