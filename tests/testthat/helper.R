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

# Reproduce the scalar log-density (transmission + augmented prior) from the
# named evaluate_log_posterior() terms, for tests asserting a single value.
eval_log_density <- function(x, house_id = NA_character_) {
  terms <- evaluate_log_posterior(x, house_id = house_id)
  terms[["transmission"]] + terms[["augm_lprior"]]
}
