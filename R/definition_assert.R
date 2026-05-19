## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Function columns associated to module family.
family2model_col <- function(family) {
  if (family == "augmented") "model"
  else "prior" # family is "parameter"
}

# Function columns associated to module family.
family2function_col <- function(family) {
  if (family == "augmented") c("model", "func")
  else if (family == "infectivity profile") "func"
  else "prior" # family is "parameter" or "group rate"
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
  for (col in family2function_col(family)) {
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

# Assert that character respect convention.
assert_valid_character <- function(x, family, arg_name) {
  if (nchar(x) == 0) {
    stop_h(sprintf("%s of %s must be a non-empty string.", arg_name, family))
  } else if (!grepl("^[A-Za-z0-9_.]+$", x)) {
    stop_h(
      sprintf("%s of %s '%s' contained invalid characters.",
              arg_name, family, x),
      "Only letters, numbers, underscores '_' and dots '.' are allowed."
    )
  }
}

# Assert that x can be treated as a feature value.
assert_valid_feature <- function(x, family, arg_name = "Feature value") {
  if (length(x) != 1 || is.na(x)) {
    stop_h(sprintf("%s of %s must be a scalar.", arg_name, family))
  } else if (!any(c(is.factor(x), is.logical(x),
                    is.numeric(x), is.character(x)))) {
    # UPDATE: is there something more robust than checking types?
    stop_h(
      sprintf("%s of %s must be a unidimensional type.", arg_name, family),
      format_type_error(x)
    )
  } else if (is.factor(x) || is.character(x)) {
    assert_valid_character(as.character(x), family, arg_name)
  }
}

# Assert that given definition name is valid.
assert_valid_definition_name <- function(x, family, is_implicit = FALSE) {
  if (length(x) != 1 || is.na(x) || !is.character(x)) {
    stop_h(
      sprintf("Name of %s must be a non-empty string.", family),
      ifelse(
        isTRUE(is_implicit),
        "Did you provide a list and forgot rlang injection operator '!!!'?",
        format_type_error(x)
      )
    )
  } else if (is.character(x)) {
    assert_valid_character(x, family, "Name")
    if (x %in% .reserved_math_consts) stop_h(
      sprintf("Name '%s' for %s is reserved for constants.", x, family)
    )
  }
}

# Assert that given status(es) name is(are) valid.
assert_valid_status <- function(status) {
  status <- status[!is.na(status)]
  status_as_int <- as.integer(status)
  if (any(status_as_int <= 0) || !all(status == status_as_int)) stop_h(
    "Status specific definition is valid only for positive integer values."
  )
}
