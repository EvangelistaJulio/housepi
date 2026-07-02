## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Augmented stores three display families.
is_augmented_family <- function(family) {
  family %in% c("augmented", "observation", "household observation")
}

# Function columns associated to holder family.
family2function_col <- function(family) {
  if (is_augmented_family(family)) c("prior", "func")
  else if (family == "infectivity profile") "func"
  else "prior"
}

# Select which candidate exposure error to surface, keeping its full
# (multi-line header + hints) text so the cli emitter can split it. Selection
# keys on the salient last line (the cause), not the "Exposing ..." header.
helper_format_error <- function(err) {
  if (length(err) == 1) return(err[[1]])
  last <- vapply(err, function(x) utils::tail(strsplit(x, "\n")[[1]], 1),
                 character(1))
  unused <- stringr::str_starts(last, "Unused|Missing")
  if (any(unused)) err[[which(unused)[[1]]]]
  else if ("func" %in% err) err["func"]
  else err[[1]]
}

# Probe a C++ exposer through the raw Rcpp module, bypassing the `.engine`
# re-styling proxy. A failure then yields the plain C++ text (glue markup, but
# no cli bullets) that `helper_format_error` classifies on; success returns
# NULL. Assigning NULL to a list slot leaves it absent, so `err` below holds
# only the failed probes.
probe_exposer_ <- function(fn, name, args) {
  member <- cpp_member_(.subset2(.engine, "module"), fn)
  tryCatch({
    member(name, unlist(args), NA_real_, NA_real_, FALSE)
    NULL
  }, error = function(e) conditionMessage(e))
}

# Trigger C++ errors
raise_def_error <- function(func_name, arg_names, def_name,
                            family, model_col) {
  # WARNING: I am assuming `1` is a safe argument value...
  arg_names <- arg_names[arg_names != "x"]
  func_args <- stats::setNames(rep(1, length(arg_names)), arg_names)

  is_infec_prof <- family == "infectivity profile"
  is_augm_func <- is_augmented_family(family) && model_col == "func"

  err <- list()
  if (is_augm_func)
    err$func <- probe_exposer_("expose_function", func_name, func_args)
  if (!is_infec_prof)
    err$dist <- probe_exposer_("expose_distribution", func_name, func_args)
  if (is_augm_func || is_infec_prof)
    err$infec <- probe_exposer_("expose_infec_profile", func_name, func_args)

  if ((is_augm_func && length(err) == 3) ||
        (!is_augm_func && length(err) == 1)) stop_h(
    sprintf("Wrong definition found for %s.",
            ifelse(is.na(def_name), family, sprintf("{.val %s}", def_name))),
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
    stop_h("{arg_name} of {family} must be a non-empty string.")
  } else if (!grepl("^[A-Za-z0-9_.]+$", x)) {
    stop_h(
      "{arg_name} of {family} {.val {x}} contained invalid characters.",
      paste("Only letters, numbers, underscores {.code _}",
            "and dots {.code .} are allowed.")
    )
  }
}

# Assert that x can be treated as a feature value.
assert_valid_feature <- function(x, family, arg_name = "Feature value") {
  if (length(x) != 1 || is.na(x)) {
    stop_h("{arg_name} of {family} must be a scalar.")
  } else if (!any(c(is.factor(x), is.logical(x),
                    is.numeric(x), is.character(x)))) {
    # UPDATE: is there something more robust than checking types?
    stop_h(
      "{arg_name} of {family} must be a unidimensional type.",
      format_type_error(x)
    )
  } else if (is.factor(x) || is.character(x)) {
    assert_valid_character(as.character(x), family, arg_name)
  }
}

# Assert that given definition name is valid.
assert_valid_definition_name <- function(x, family) {
  if (length(x) != 1 || is.na(x) || !is.character(x)) {
    stop_h(
      "Name of {family} must be a non-empty string.",
      format_type_error(x)
    )
  } else {
    assert_valid_character(x, family, "Name")
    if (x %in% .reserved_math_consts) stop_h(
      "Name {.val {x}} for {family} is reserved for constants."
    )
  }
}

# Assert that given status(es) name is(are) valid.
# Numeric statuses must be positive integers; character labels are allowed and
# resolved to integer codes against the status dictionary at build time.
assert_valid_status <- function(status) {
  status <- status[!is.na(status)]
  if (length(status) == 0) return(invisible())
  if (is.numeric(status) &&
        (any(status <= 0) || any(status != as.integer(status)))) stop_h(
    "Status specific definition is valid only for positive integer values."
  )
  invisible()
}
