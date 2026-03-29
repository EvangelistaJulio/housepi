
# Helper to load an associated dataset internally.
helper_test_get_data <- function(data_name) {
  env <- new.env(parent = emptyenv())
  utils::data(list = data_name, package = "housepi", envir = env)
  env[[data_name]]
}

# Helper to load definition lists.
helper_test_aux_load <- function(self, ..., append = FALSE) {
  self$load_quosures(rlang::enquos(...), append = append)
  invisible(self)
}
