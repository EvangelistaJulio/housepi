## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
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

# Select which candidate exposure error to surface, keeping its full
# (multi-line header + hints) text so the cli emitter can split it. Selection
# keys on the salient last line (the cause), not the "Exposing ..." header.
helper_format_error <- function(err) {
  if (length(err) == 1) return(err[[1]])
  last <- vapply(err, \(x) utils::tail(strsplit(x, "\n")[[1]], 1),
                 character(1))
  unused <- grepl(paste0("^(", bullet_tag_, ".)?(Unused|Missing)"), last)
  if (any(unused)) err[[which(unused)[[1]]]]
  else if ("func" %in% err) err["func"]
  else err[[1]]
}

# Probe a C++ exposer through the raw Rcpp module, bypassing the `.engine`
# re-styling proxy. A failure then yields the plain C++ text (glue markup, but
# no cli bullets) that `helper_format_error` classifies on; success returns
# the exposer.
probe_exposer_ <- function(fn, name, args, min = NA_real_, max = NA_real_) {
  member <- cpp_member_(.subset2(.engine, "module"), fn)
  tryCatch(member(name, unlist(args), min, max, FALSE),
           error = \(e) conditionMessage(e))
}

# Value a readable definition entry is probed at: a constant at its own value,
# anything sourced or computed (an encoding carrying '|') at `placeholder`.
probe_value_ <- function(x, placeholder) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(x))
  out[is.na(x) | grepl("|", x, fixed = TRUE) | is.na(out)] <- placeholder
  out
}

# Build the definition in every registry the engine queries for this column
# and family, and raise the C++ error when none accepts it. Return the exposer
# when a single registry accepts it, NULL otherwise.
# Constant arguments and bounds are probed at their value, so that a value the
# engine rejects errors at definition. Other arguments are probed at `1`; the
# bounds are left out unless all given ones are constant.
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
  out <- lapply(stats::setNames(nm = probed), \(k) {
    probe_exposer_(exposers[[k]], func_name, func_args, bounds[1], bounds[2])
  })
  err <- Filter(is.character, out)

  # every registry refused the name: the definition has no reading
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

# Write into row `r` of the definitions the defaults the user left out: the
# arguments of the parametrization in use, and the bounds that hold a value of
# their own. A bound reading an argument is left out, as passing it would
# detach the bound from that argument.
fill_defaults_ <- function(gear, r, col, obj, user_args) {
  args <- cpp_member_(obj, "args")()
  args <- args[!is.na(args) & !names(args) %in% user_args]
  if (length(args)) {
    gear$defs[r, col] <- append_call_args_(gear$defs[r, col], args)
  }

  domain <- cpp_member_(obj, "domain")()
  aliases <- cpp_member_(obj, "bound_aliases")()
  for (key in intersect(c("min", "max"), colnames(gear$defs))) {
    bound <- domain[[key]]
    if (is.na(gear$defs[r, key]) && is.finite(bound) && !aliases[[key]]) {
      gear$defs[r, key] <- sprintf("%s|N", as.character(bound))
    }
  }
}

# Assert that arguments provided to model are legitimate, and write the
# defaults the user left out into the definitions. `last` restricts both to
# that many trailing rows, all rows when NA.
# UPDATE: at the moment latent functions are not testable
# because my registries are exposed separately
housepi_assert_args <- function(gear, last = NA_integer_) {
  rows <- seq_len(NROW(gear$defs))
  if (!is.na(last)) rows <- utils::tail(rows, last)
  family <- gear$family

  # load requested faimily definitions
  for (col in family2function_col(family)) {
    col_args <- if (col == "func") "func_args" else "args"

    readable_defs <- housepi_defin_showcaser(gear$defs[rows, , drop = FALSE])
    func_names <- readable_defs[[col]]

    for (i in which(!is.na(func_names))) {
      r <- rows[[i]]
      # the bounds truncate the defining distribution, never the function
      bound <- \(key) {
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
  if (!length(dups)) return(invisible())
  stop_h(paste("{cli::qty(dups)}Repeated argument{?s} {.arg {dups}} for",
               if (is.null(name)) "{family}." else "{family} {.val {name}}."))
}

# Assert a name carries no declaration of another kind:
# latent variables and observations share one definition store,
# and each name stands for a single quantity of the model.
assert_no_kind_clash <- function(defs, name, kind) {
  if (is.null(defs)) return(invisible())
  other <- defs[["kind"]][defs[["name"]] == name & defs[["kind"]] != kind]
  if (!length(other)) return(invisible())
  stop_h(
    "Name {.val {name}} is already declared as {kind2family(other[1])}.",
    "Each name stands for one quantity: choose a different name."
  )
}

# Assert that a string argument respects the naming convention, `what` naming
# it in the message.
assert_valid_character <- function(x, what) {
  if (!nzchar(x)) stop_h("{what} must be a non-empty string.")
  else if (!grepl("^[A-Za-z0-9_.]+$", x)) stop_h(
    "{what} {.val {x}} contained invalid characters.",
    "Only letters, numbers, underscores and dots are allowed."
  )
}

# Assert a given law name respects the naming convention.
assert_valid_law_name <- function(x, what) {
  if (is.character(x) && length(x) == 1 && !is.na(x)) {
    assert_valid_character(x, what)
  }
}

# Assert that x can be treated as a covariate value.
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
}

# Assert that given definition name is valid.
# A leading dot is reserved to the definitions the engine generates itself.
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
}

# Assert data columns have valid names and can be sourced by definitions.
# Eventually warn about columns similar to missing reserved columns.
# The engine reads reserved columns by exact name: lookalike stays a covariate.
assert_valid_data_columns <- function(cols, targets = data_reserved_cols) {
  forbidden <- cols[!grepl("^[A-Za-z0-9_.]+$", cols)]
  if (length(forbidden)) stop_h(
    "Column {cli::qty(forbidden)}name{?s} are forbidden: {.field {forbidden}}.",
    "Use only letters, digits and/or underscores {.code _}."
  )
  # warn about lookalikes
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

# Assert that given status(es) name is(are) valid.
# Numeric statuses must be positive integers; character labels are allowed and
# resolved to integer codes against the status dictionary at build time.
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
