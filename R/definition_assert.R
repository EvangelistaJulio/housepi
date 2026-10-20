## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Kinds of definition the augmented store holds.
definition_kinds_ <- c("latent", "member", "house")

# Display family of a definition kind.
kind2family <- function(kind) {
  switch(kind,
    "member" = "observation",
    "house" = "household observation",
    "latent"
  )
}

# Reserved data columns, each with its possible space-free spelled mistakes.
data_reserved_aliases <- list(
  id = outer(c("id", "ind", "indiv", "individual"),
             c("", "id"), paste0),
  hh_id = outer(c("h", "hh", "house", "household"),
                c("", "id"), paste0),
  status = outer(c("status", "inf", "infec", "infection"),
                 c("", "status"), paste0)
)
# Reserved column names (non-covariate) in a logic table.
data_reserved_cols <- names(data_reserved_aliases)

# Augmented stores one display family per kind.
is_latent_family <- function(family) {
  family %in% vapply(definition_kinds_, kind2family, character(1))
}

# Column carrying the defining distribution of the holder family.
family2model_col <- function(family) {
  if (is_latent_family(family)) "model" else "prior"
}

# Function columns associated to holder family.
family2function_col <- function(family) {
  if (is_latent_family(family)) c("model", "func")
  else if (family == "infectivity profile") "func"
  else if (family == "parameter") c("prior", "func")
  else "prior"
}

# Select the exposure error to surface, keyed on its last line, the cause.
# The text is kept whole, header and hints, for the cli emitter to split.
helper_format_error <- function(err) {
  if (length(err) == 1) return(err[[1]])
  last <- vapply(err, \(x) utils::tail(strsplit(x, "\n")[[1]], 1),
                 character(1))
  unused <- grepl(paste0("^(", bullet_tag_, ".)?(Unused|Missing)"), last)
  if (any(unused)) err[[which(unused)[[1]]]]
  else if ("func" %in% err) err["func"]
  else err[[1]]
}

# Probe a C++ exposer through the raw Rcpp module, bypassing `.engine`.
# * failure -> plain C++ text, glue markup without cli bullets,
#   which helper_format_error() classifies;
# * success -> the exposer.
probe_exposer_ <- function(fn, name, args, min = NA_real_, max = NA_real_) {
  member <- cpp_member_(.subset2(.engine, "module"), fn)
  tryCatch(member(name, unlist(args), min, max, FALSE),
           error = \(e) conditionMessage(e))
}

# Value a readable definition entry is probed at:
# * constant -> its own value;
# * sourced or computed, an encoding carrying '|' -> `placeholder`.
probe_value_ <- function(x, placeholder) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  out[is.na(x) | grepl("|", x, fixed = TRUE) | is.na(out)] <- placeholder
  out
}

# Build the definition in every registry the engine queries for this column
# and family, raising the C++ error when none accepts it.
# * one registry accepts it -> the exposer, else NULL;
# * constant argument or bound -> probed at its value, rejected at definition;
# * other argument -> probed at `1`;
# * bounds -> left out unless every given one is constant.
probe_definition_ <- function(func_name, args, def_name, family, model_col,
                              min = NA_real_, max = NA_real_) {
  # a definition without arguments reads back as an unnamed NA
  args <- if (is.null(names(args))) numeric() else args[names(args) != "x"]
  func_args <- stats::setNames(probe_value_(args, 1), names(args))
  # a constant bound beside a sourced one would meet the default in its place
  bounds <- probe_value_(c(min, max), NA_real_)
  if (any(is.na(bounds) & !is.na(c(min, max)))) bounds[] <- NA_real_

  probed <- if (!identical(model_col, "func")) "dist"
  else if (family == "infectivity profile") "infec"
  else if (is_latent_family(family)) c("func", "dist", "infec")
  else "func"

  exposers <- c(func = "expose_function",
                dist = "expose_distribution",
                infec = "expose_infec_profile")
  out <- lapply(stats::setNames(nm = probed), function(k) {
    probe_exposer_(exposers[[k]], func_name, func_args, bounds[1], bounds[2])
  })
  err <- Filter(is.character, out)

  # every registry refused the name -> the definition has no reading
  if (length(err) == length(probed)) stop_h(
    if (is.na(def_name)) "Wrong definition found for {family}."
    else "Wrong definition found for {.val {def_name}}.",
    helper_format_error(err)
  )
  if (length(probed) - length(err) == 1) Filter(Negate(is.character), out)[[1]]
}

# Append named numeric arguments to an encoded call: "f" or "f|a|[..]".
append_call_args_ <- function(call, args) {
  vals <- stats::setNames(sprintf("%s|N", as.character(args)), names(args))
  if (!housepi_is_call(call)) return(housepi_combine_to_call(call, vals))
  parts <- regmatches(call, regexec("^([^|]*)\\|([^|]*)\\|(.*)$", call))[[1]]
  sprintf("%s|%s:%s|%s:%s", parts[2], parts[3],
          paste(names(args), collapse = ":"), parts[4],
          paste(sprintf("[%s]", vals), collapse = ":"))
}

# Write into row `r` the defaults the user left out:
# * arguments of the parametrization in use;
# * bounds holding a value of their own, not those reading an argument,
#   which passing would detach from it.
fill_defaults_ <- function(gear, r, col, obj, user_args) {
  args <- cpp_member_(obj, "args")()
  args <- args[!is.na(args) & !names(args) %in% user_args]
  if (length(args))
    gear$defs[r, col] <- append_call_args_(gear$defs[r, col], args)

  domain <- cpp_member_(obj, "domain")()
  aliases <- cpp_member_(obj, "bound_aliases")()
  for (key in intersect(c("min", "max"), colnames(gear$defs))) {
    bound <- domain[[key]]
    if (is.na(gear$defs[r, key]) && is.finite(bound) && !aliases[[key]])
      gear$defs[r, key] <- sprintf("%s|N", as.character(bound))
  }
}

# Assert the arguments given to a model are legitimate, and write the defaults
# the user left out into the definitions.
# * `last` -> only that many trailing rows;
# * `last` NA -> every row.
# @todo latent functions are not testable, their registries are exposed apart
housepi_assert_args <- function(gear, last = NA_integer_) {
  rows <- seq_len(NROW(gear$defs))
  if (!is.na(last)) rows <- utils::tail(rows, last)
  family <- gear$family

  # load the definitions of the requested family
  for (col in family2function_col(family)) {
    col_args <- if (col == "func") "func_args" else "args"

    readable_defs <- housepi_defin_showcaser(gear$defs[rows, , drop = FALSE])
    func_names <- readable_defs[[col]]

    # an alias reads a definition, not a function
    is_alias <- vapply(gear$defs[rows, col], housepi_is_sourced_value,
                       logical(1))
    for (i in which(!is.na(func_names) & !is_alias)) {
      r <- rows[[i]]
      # the bounds truncate the defining distribution, never the function
      bound <- function(key) {
        if (col == "func" || !key %in% colnames(readable_defs)) NA_real_
        else readable_defs[[key]][[i]]
      }
      user_args <- readable_defs[i, col_args][[1]]
      obj <- probe_definition_(
        func_names[[i]], user_args,
        if ("name" %in% colnames(gear$defs)) gear$defs[r, "name"] else NA,
        family, col, min = bound("min"), max = bound("max")
      )
      if (!is.null(obj))
        fill_defaults_(gear, r, col, obj, names(user_args))
    }
  }
}

# Assert no argument name is given twice to one definition.
assert_unique_arg_names <- function(nms, family, name = NULL) {
  dups <- unique(nms[nzchar(nms) & duplicated(nms)])
  if (length(dups)) stop_h(
    paste("{cli::qty(dups)}Repeated argument{?s} {.arg {dups}} for",
          if (is.null(name)) "{family}." else "{family} {.val {name}}.")
  )
  invisible()
}

# Assert a name carries no declaration of another kind.
# Latents and observations share one store, where a name is one quantity.
assert_no_kind_clash <- function(defs, name, kind) {
  if (is.null(defs)) return(invisible())
  other <- defs[["kind"]][defs[["name"]] == name & defs[["kind"]] != kind]
  if (length(other)) stop_h(
    "Name {.val {name}} is already declared as {kind2family(other[1])}.",
    "Each name stands for one quantity: choose a different name."
  )
  invisible()
}

# Assert a string argument respects the naming convention.
# @param what `character(1)` Argument name in the message.
assert_valid_character <- function(x, what) {
  if (!nzchar(x)) stop_h("{what} must be a non-empty string.")
  else if (!grepl("^[A-Za-z0-9_.]+$", x)) stop_h(
    "{what} {.val {x}} contained invalid characters.",
    "Only letters, numbers, underscores and dots are allowed."
  )
  invisible()
}

# Assert a distribution name respects the naming convention.
assert_valid_law_name <- function(x, what) {
  if (is.character(x) && length(x) == 1 && !is.na(x))
    assert_valid_character(x, what)
  invisible()
}

# Assert `x` can be a covariate value.
assert_valid_covariate <- function(x, family) {
  arg_name <- "Covariate value"
  if (length(x) != 1 || is.na(x)) {
    stop_h("{arg_name} of {family} must be a scalar.")
  } else if (is.factor(x) || is.character(x)) {
    assert_valid_character(as.character(x), paste(arg_name, "of", family))
  } else if (!is.logical(x) && !is.numeric(x)) {
    stop_h("{arg_name} of {family} must be a unidimensional type.",
           format_type_error(x))
  }
  invisible()
}

# Assert a definition name is valid.
# * leading dot -> reserved to definitions the engine generates.
assert_valid_definition_name <- function(x, family) {
  if (length(x) != 1 || is.na(x) || !is.character(x)) {
    stop_h("Name of {family} must be a non-empty string.",
           format_type_error(x))
  } else {
    assert_valid_character(x, paste("Name of", family))
    if (x %in% .reserved_math_consts) stop_h(
      "Name {.val {x}} for {family} is reserved for constants."
    )
    if (startsWith(x, ".")) stop_h(
      "Name {.val {x}} of {family} cannot start with a dot.",
      "A leading dot is reserved to the definitions housepi generates itself."
    )
  }
  invisible()
}

# Assert data columns have valid names and can be sourced by definitions.
# * column resembling a missing reserved one -> warning, it stays a covariate.
assert_valid_data_columns <- function(cols, targets = data_reserved_cols) {
  forbidden <- cols[!grepl("^[A-Za-z0-9_.]+$", cols)]
  if (length(forbidden)) stop_h(
    "Column {cli::qty(forbidden)}name{?s} are forbidden: {.field {forbidden}}.",
    "Use only letters, digits and/or underscores {.code _}."
  )
  # warn about lookalikes
  # nolint next: object_usage_linter. read by the cli glue below
  what <- c(id = "individual identifier",
            hh_id = "household identifier",
            status = "infection status")
  flat <- gsub(" ", "", spaced_choice(cols), fixed = TRUE)
  for (target in setdiff(targets, cols)) {
    hits <- cols[flat %in% data_reserved_aliases[[target]]]
    if (length(hits)) warn_h(
      paste("{cli::qty(hits)}Found column{?s} {.field {hits}} in data",
            "recalling a {what[[target]]}."),
      paste("If {cli::qty(hits)}{?it/one of them} is meant to be used as",
            "a {what[[target]]}, rename it to {.field {target}}.")
    )
  }
  invisible()
}

# Highest infected status level the engine stores.
# StatusLevel is an unsigned 8-bit integer (see inst/include/data/status.h).
max_status_level <- 254L

# Assert statuses are valid.
# * numeric -> positive integers, at most `max_status_level`;
# * character -> labels, resolved to codes at build time.
assert_valid_status <- function(status) {
  status <- status[!is.na(status)]
  if (!length(status)) return(invisible())
  if (!is.numeric(status)) return(invisible())
  if (any(status <= 0) || any(status != as.integer(status))) stop_h(
    "Status specific definition is valid only for positive integer values."
  )
  if (any(status > max_status_level)) stop_h(
    paste("Status specific definition exceeded the maximum of",
          "{max_status_level} infected statuses.")
  )
  invisible()
}
