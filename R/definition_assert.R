## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

family2model_col <- function(family) {
  if (family == "augmented") c("model", "func")
  else if (family == "parameter") "prior"
  else if (family == "infectivity profile") "func"
  else if (family == "group rate") "prior"
  else stop_h(sprintf(
    "[Internal] Unrecognized module family '%s'.", family
  ))
}

helper_format_error <- function(err) {
  err <- lapply(err, function(x) utils::tail(strsplit(x, "\n")[[1]], 1))
  if (length(err) == 1) {
    sms <- err[[1]]
  } else {
    unused <- stringr::str_starts(err, "Unused|Missing")
    if (any(unused)) sms <- err[[which(unused)[[1]]]]
    else if ("func" %in% err) sms <- err["func"]
    else sms <- err[[1]]
  }
  sms
}

# Trigger C++ errors
# UPDATE: I should find an hollow export mechanism instead
raise_def_error <- function(func_name, arg_names, def_name,
                            family, model_col) {
  # WARNING: I am assuming `1` is a safe argument value...
  arg_names <- arg_names[arg_names != "x"]
  func_args <- stats::setNames(rep(1, length(arg_names)), arg_names)

  is_infec_prof <- family == "infectivity profile"
  is_augm_func <- family == "augmented" && model_col == "func"

  err <- list()
  if (is_augm_func) tryCatch(
    expose_function(func_name, func_args),
    error = function(e) err$func <<- e$message
  )
  if (!is_infec_prof) tryCatch(
    expose_distribution(func_name, func_args),
    error = function(e) err$dist <<- e$message
  )
  if (is_augm_func || is_infec_prof) tryCatch(
    expose_infec_profile(func_name, func_args),
    error = function(e) err$infec <<- e$message
  )

  if ((is_augm_func && length(err) == 3) ||
        (!is_augm_func && length(err) == 1)) stop_h(
    sprintf("Wrong definition found for %s.",
            ifelse(is.na(def_name), family, sprintf("'%s'", def_name))),
    helper_format_error(err)
  )
}

# Assert that arguments provided to model are legitimate.
# UPDATE: at the moment augmented functions are not testable
# because my registries are exposed separately
housepi_assert_args <- function(gear, last = FALSE) {
  defs <- gear$defs
  if (isTRUE(last)) defs <- utils::tail(defs, 1)
  family <- gear$family

  # load requested faimily definitions
  for (col in family2model_col(family)) {
    col_args <- ifelse(col == "func", "func_args", "args")

    readable_defs <- housepi_defin_showcaser(defs)
    func_names <- readable_defs[[col]]

    for (r in which(!is.na(func_names))) raise_def_error(
      func_names[[r]],
      names(readable_defs[r, col_args][[1]]),
      if ("name" %in% colnames(defs)) defs[r, "name"] else NA,
      family, col
    )
  }
}

# Assert that given definition name is valid.
assert_valid_definition_name <- function(x, family, is_implicit = FALSE) {
  sms <- c()
  if (length(x) != 1 || is.na(x) || !is.character(x) || nchar(x) == 0) {
    sms <- c(
      sprintf("Definition name of %s must be a non-empty character(1).", family)
    )
    if (isTRUE(is_implicit)) sms <- c(sms,
      "Did you provide a list and forgot rlang injection operator `!!!`?"
    )
  } else if (!grepl("^[A-Za-z0-9_.]+$", x)) {
    sms <- c(
      sprintf("Definition name of %s '%s' contained invalid characters.",
              family, x),
      "Only letters, numbers, underscores '_' and dots '.' are allowed."
    )
  } else if (x %in% .reserved_math_consts) {
    sms <- sprintf("Definition name of %s '%s' is reserved for constants.",
                   family, x)
  }
  if (!is.null(sms)) stop_h(sms)
}

# Assert that given status(es) name is(are) valid.
assert_valid_status <- function(status) {
  status <- status[!is.na(status)]
  status_as_int <- as.integer(status)
  if (any(status_as_int <= 0) || !all(status == status_as_int)) stop_h(
    "Status specific definition is valid only for positive integer values."
  )
}
