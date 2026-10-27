## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Unprintable tag leading a C++ message line, followed by its bullet name.
bullet_tag_ <- "\x1F"

# Render a `cli_bullets()` block as one message.
print_bullets <- function(x, .envir = parent.frame()) {
  message(paste(cli::cli_fmt(cli::cli_bullets(x, .envir = .envir)),
                collapse = "\n"))
}

# Render `header` as an info line,
# then one indented `Label: value` line per element of `rows`,
# its name the label and its value cli markup.
# * `warnings` with rows -> closing `Warning` line pointing to `$warnings`.
print_rows_ <- function(header, rows = character(), warnings = NULL,
                        .envir = parent.frame()) {
  lines <- paste0(names(rows), ": ", rows)
  n_warn <- NROW(warnings)
  warn <- if (n_warn) c(`!` = paste0(
    cli::col_yellow("Warning"), ": ",
    inline_cli_("{n_warn} issue{?s} flagged, see {.code $warnings}.")
  ))
  print_bullets(c(i = header, stats::setNames(lines, rep(" ", length(lines))),
                  warn), .envir = .envir)
}

# Duration of `secs` seconds, as the progress bar prints it.
format_secs_ <- function(secs) .engine$format_secs(secs)

# Flatten `...` and vectors into one line per element, splitting on newlines.
# Names are kept.
# * segment produced by a split -> empty name,
#   unless it opens with `bullet_tag_`, whose next character names it.
flatten_cli_parts_ <- function(dots) {
  parts <- unlist(dots, use.names = TRUE)
  nms <- names(parts)
  if (is.null(nms)) nms <- rep("", length(parts))
  out <- character()
  out_nms <- character()
  for (k in seq_along(parts)) {
    segs <- strsplit(parts[[k]], "\n", fixed = TRUE)[[1]]
    seg_nms <- c(nms[[k]], rep("", max(0, length(segs) - 1)))[seq_along(segs)]
    tagged <- startsWith(segs, bullet_tag_)
    seg_nms[tagged] <- substr(segs[tagged], 2, 2)
    segs[tagged] <- substring(segs[tagged], 3)
    out <- c(out, segs)
    out_nms <- c(out_nms, seg_nms)
  }
  stats::setNames(out, out_nms)
}

# Housepi warnings are held back while the switch is on.
quiet_warnings_ <- new.env(parent = emptyenv())
quiet_warnings_$on <- FALSE

#' Value of `code(bar)`, `bar` a progress bar over `total` units
#' printed every `progress` seconds.
#' * `progress` FALSE or 0 -> `bar` is `NULL`;
#' * `progress` TRUE -> one print per second;
#' * knitted document -> closing line alone.
#' Leaving `code` closes the bar:
#' * returned -> with its closing line;
#' * failed -> its line ended, before the error prints;
#' * interrupted -> as left, R ends the line.
#' @param label `character(1)` Name heading the bar.
#' @param code `function(bar)` Code filling the bar.
#' @noRd
with_progress_ <- function(progress, label, total, code) {
  pace <- as.double(progress[[1L]])
  if (is.na(pace) || pace <= 0) return(code(NULL))
  # a knitted document ignores carriage returns, so its prints pile up
  document <- isTRUE(getOption("knitr.in.progress"))
  # half the console, at most 39 symbols: squares on UTF-8 output, "#" otherwise
  bar <- .engine$progress_open(label, total, pace,
                               min(getOption("width"), 78L) %/% 2L,
                               interactive(), document,
                               cli::num_ansi_colors(stderr()) > 1L,
                               if (cli::is_utf8_output()) "\u25a0" else "#")
  out <- withCallingHandlers(code(bar),
                             error = \(e) .engine$progress_crash(bar))
  .engine$progress_close(bar)
  out
}

#' Hold back housepi warnings for the rest of the calling frame.
#' @param quiet `logical(1)` Whether to hold the switch, `FALSE` doing nothing.
#' @param .envir `environment` Frame the switch is released at the end of.
#' @noRd
local_quiet_warnings_ <- function(quiet, .envir = parent.frame()) {
  if (!isTRUE(quiet)) return(invisible())
  held <- quiet_warnings_$on
  quiet_warnings_$on <- TRUE
  withr::defer(quiet_warnings_$on <- held, envir = .envir)
  invisible()
}

#' Emit a harmonized cli error or warning condition.
#' * first element -> header, `x` for errors, unnamed for warnings;
#' * other element without a bullet name -> hint.
#' @param parts `character()` Message lines, optionally named with bullet keys.
#' @param level `character(1)` Either `"stop"` or `"warn"`.
#' @param trusted `logical(1)` Whether the text is trusted.
#' If `FALSE`, braces are escaped so foreign text cannot trip cli glue.
#' @param .envir `environment` Frame used to resolve `{ }` interpolations.
#' @noRd
emit_cli_condition <- function(parts, level = c("stop", "warn"),
                               trusted = TRUE, .envir = parent.frame()) {
  parts <- flatten_cli_parts_(parts)
  level <- match.arg(level)

  # foreign text carries no markup: doubling renders its braces literally
  if (!trusted) parts[] <- escape_cli_(parts)

  # capitalize the leading letter of each line
  parts[] <- sub("^([a-z])", "\\U\\1", parts, perl = TRUE)

  # un-named continuation lines become info hints
  nms <- names(parts)
  nms[nms == ""] <- "i"
  if (length(nms)) nms[[1]] <- if (level == "stop") "x" else ""
  bullets <- stats::setNames(parts, nms)

  if (level == "stop") cli::cli_abort(bullets, call = NULL, .envir = .envir)
  # skip warnings if disabled
  else if (quiet_warnings_$on) invisible()
  # display at emission time, so warnings are not deferred past the errors
  else withr::with_options(list(warn = 1), cli::cli_warn(
    bullets, call = NULL, .envir = .envir
  ))
}

#' Stop execution with a harmonized, cli-rendered error.
#' @param ... Message lines: the header, then hints `i` unless named otherwise.
#' @noRd
stop_h <- function(...) {
  emit_cli_condition(list(...), "stop", .envir = parent.frame())
}

#' Emit a harmonized, cli-rendered warning.
#' @param ... Message lines, with the bullet convention of stop_h().
#' @noRd
warn_h <- function(...) {
  emit_cli_condition(list(...), "warn", .envir = parent.frame())
}

#' Restyle any condition raised beneath `expr` through the cli emitter,
#' split on newlines into a header and hints.
#' * C++ condition -> recognised by class, rendered as cli markup;
#' * rlang condition, raised by an R callback -> passed on, its bullets kept;
#' * foreign condition -> shown literally.
#' @noRd
with_cpp_conditions <- function(expr) {
  reemit <- function(cnd, level, trusted) {
    msg <- conditionMessage(cnd)
    emit_cli_condition(list(msg), level, trusted = trusted)
  }
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      if (inherits(e, "rlang_error")) stop(e)
      reemit(e, "stop", inherits(e, "housepi::TrustedException"))
    }),
    warning = function(w) {
      if (inherits(w, "rlang_warning")) return()
      reemit(w, "warn", inherits(w, "housepi_cpp_warning"))
      invokeRestart("muffleWarning")
    }
  )
}

#' Stop when a value is assigned to a read-only active field.
#' @param value The `value` argument of the active binding, missing on read.
#' @param active `character(1)` Name of the active method.
#' @noRd
assert_read_only_ <- function(value, active) {
  if (!missing(value)) stop_h("Active field {.field {active}} is read only.")
  invisible()
}

# Unprintable separator hashing several values into one key.
key_sep_ <- "\x1F"

#' Hash the `cols` of each row of `x` into one key.
#' @noRd
row_keys_ <- function(x, cols) {
  # no column: every row shares one key
  if (!length(cols)) return(rep(key_sep_, nrow(x)))
  do.call(paste, c(lapply(x[cols], as.character), list(sep = key_sep_)))
}

#' Stop on a non-empty message when `wall` is set, else return it.
#' @noRd
wall_or_return_ <- function(sms, wall) {
  if (isTRUE(wall) && length(sms)) stop_h(sms) else sms
}

#' Thin draws to `n` evenly spaced elements, or rows of a `data.frame`,
#' first and last kept.
#' Spacing them out keeps autocorrelated draws as far apart as they can be.
#' * `n` NA or at least `NROW(x)` -> `x` unchanged.
#' @param x Vector or `data.frame` of draws.
#' @param n `integer(1)` Number of draws kept.
#' @noRd
helper_thin_even_ <- function(x, n) {
  if (is.na(n) || NROW(x) <= n) return(x)
  keep <- round(seq.int(1L, NROW(x), length.out = n))
  if (is.data.frame(x)) x[keep, , drop = FALSE] else x[keep]
}

#' Seeds of the engine random number generators, drawn from the R stream.
#' @param n `integer(1)` Number of seeds.
#' @noRd
draw_seed_ <- function(n = 1L) {
  sample.int(.Machine$integer.max, n)
}

#' A count argument as a whole positive number.
#' @param what `character(1)` What the count is, opening the error.
#' @param allow_na `logical(1)` Let a single `NA` through, as `NA_integer_`.
#' @noRd
assert_count_ <- function(x, what, allow_na = FALSE) {
  if (allow_na && length(x) == 1L && is.na(x)) return(NA_integer_)
  if (!rlang::is_scalar_integerish(x) || is.na(x) || x < 1) stop_h(
    "{what} should be a positive integer."
  )
  as.integer(x)
}

#' Read a parameter table from draws or a data.frame, thinned to the ESS.
#' @param thin_ess `logical(1)` Whether to thin to the ESS.
#' @noRd
helper_read_param_table_ <- function(ptab, thin_ess = TRUE) {
  # extract parameter draws from table
  if (inherits(ptab, "HousepiPosteriorDraws")) {
    ptab <- draws_block_(ptab, "param")
  }

  if (posterior::is_draws(ptab)) {
    ptab <- posterior::as_draws_matrix(ptab)
    # reduce sample to effective sample size (ESS)
    if (isTRUE(thin_ess)) {
      # an ESS above the draw count keeps every draw, capped or not
      # degenerate draws give NA, asserted below
      ess_bulk <- suppressWarnings(posterior::summarise_draws(
        ptab, ess_bulk = posterior::ess_bulk
      )$ess_bulk)
      # check no NA values in ESS
      if (anyNA(ess_bulk)) stop_h(
        "ESS values return NA.",
        "Check your draws for degenerate parameters."
      )
      # keep at most ESS draws, evenly spaced over the rows
      ptab <- helper_thin_even_(as.data.frame(ptab), floor(min(ess_bulk)))
    }
  }

  tryCatch({
    return(as.data.frame(ptab))
  }, error = \(e) stop_h("Type of {.arg draws} not supported.")
  )
}

#' Split the `id_draws` selection into individual and household identifiers.
#' * `NULL` -> nothing stored;
#' * `NA` or empty -> the default household is followed.
#' Status filtering and the default household pick are left to the engine.
#' @param id_draws Selection given to [run_inference()].
#' @noRd
resolve_id_draws_ <- function(id_draws) {
  out <- list(id = character(), hh_id = character())
  # explicit opt-out: forward the sentinel so the engine stores nothing
  if (is.null(id_draws)) {
    out$id <- ".none"
    return(out)
  }
  # unspecified: let the engine follow one random household
  if (!length(id_draws) ||
        (!is.list(id_draws) && length(id_draws) == 1L && is.na(id_draws))) {
    return(out)
  }
  # named list: split individual and household identifiers
  if (is.list(id_draws)) {
    extra <- setdiff(names(id_draws), c("id", "hh_id"))
    if (length(extra)) stop_h(
      "Unknown {.arg id_draws} name{?s}: {.val {extra}}."
    )
    if (!is.null(id_draws$id)) out$id <- as.character(id_draws$id)
    if (!is.null(id_draws$hh_id))
      out$hh_id <- as.character(id_draws$hh_id)
    return(out)
  }
  # bare vector: individual identifiers
  out$id <- as.character(id_draws)
  out
}

#' Double the braces of plain text, so a cli pass prints it verbatim.
#' @param x `character()` Text carrying no markup.
#' @noRd
escape_cli_ <- function(x) {
  gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
}

#' Render cli markup now, into text that shows verbatim when emitted later.
#' Values are interpolated from `.envir` as data, never parsed as markup.
#' @param ... Message text, pasted together.
#' @param .envir `environment` Frame used to resolve `{ }` interpolations.
#' @return `character(1)`
#' @noRd
inline_cli_ <- function(..., .envir = parent.frame()) {
  escape_cli_(cli::format_inline(..., .envir = .envir))
}

# Lower-case name with `-`, `_`, `.` and runs of spaces folded to one space.
spaced_choice <- function(x) {
  trimws(tolower(gsub("\\s+", " ", gsub("[-_. ]", " ", x))))
}

#' Suggest the option closest to `input`, as a hint.
#' The accepted distance is capped by half the length of the input.
#' @param min_dist `integer(1)` Threshold for minimal distance.
#' @param otherwise `character()` Hint when no option is close,
#'   naming where the options are shown.
#' @noRd
find_closest_match <- function(input, options, min_dist = 3,
                               otherwise = NULL) {
  if (!length(options) || min_dist <= 0) return(otherwise)
  dists <- utils::adist(input, options)
  # shorter inputs tolerate fewer edits
  idxs <- apply(dists, 1, min) <= pmin(min_dist, nchar(input) %/% 2)
  if (!any(idxs)) return(otherwise)
  closest <- options[apply(dists[idxs, , drop = FALSE], 1, which.min)]
  vapply(closest, \(x) inline_cli_("Did you mean {.val {x}}?"),
         character(1), USE.NAMES = FALSE)
}

# Format an invalid-type error.
format_type_error <- function(x) {
  inline_cli_("Got an object of class {.cls {rev(class(x))}} instead.")
}

raise_permissions_error <- function() {
  stop_h("Failed to generate necessary files.",
         "Check your user permissions.")
}

# file.copy() with an overwrite option.
overwrite <- function(from, to, ...) {
  if (from == to) return(invisible())
  # user permission error, typically
  if (!suppressWarnings(file.copy(from, to, overwrite = TRUE, ...))) {
    raise_permissions_error()
  }
}

########  Files  ########

# File shipped by a model for each templated definition family.
template_def_files <- c(parameter = "param.R",
                        latent = "latent.R",
                        `infectivity profile` = "infecprof.R",
                        `status process` = "status.R")

# Deep clone R6 objects.
helper_clone <- function(x) {
  if (is.data.frame(x)) x[TRUE, , drop = FALSE]
  else if (is.list(x)) lapply(x, helper_clone)
  else if (inherits(x, "R6")) x$clone(deep = TRUE)
  else x
}

########  Settings  ########

#' Assert undocumented advanced settings overrides from `...`.
#' @param allowed `character()` Accepted argument names.
#' @param ... Advanced settings overrides.
#' @return `list` Non-null overrides keyed by engine settings key.
#' @noRd
assert_advanced_settings_ <- function(allowed, ...) {
  dots <- list(...)
  for (nm in names(dots)) {
    if (!nm %in% allowed) stop_h(
      "Unexpected argument {.arg {nm}} passed through {.code ...}.",
      find_closest_match(nm, setdiff(allowed, names(dots)))
    )
  }
  invisible()
}

#' Match `value` to one of `choices`, exactly.
#' * `value` identical to `choices` -> the first, as for an unset formal.
#' @param arg `character(1)` Argument name, as printed.
#' @param help `character(1)` Help topic documenting the choices.
#' @noRd
match_choice_ <- function(value, choices, arg, help) {
  if (identical(value, choices)) return(choices[[1]])
  if (rlang::is_string(value) && value %in% choices) return(value)
  pointer <- paste0("See {.run [?{help}](utils::help(\"{help}\", ",
                    "\"housepi\"))} for the accepted values.")
  stop_h(
    "Unknown {.arg {arg}} {.val {value}}.",
    if (rlang::is_string(value))
      find_closest_match(value, choices, otherwise = pointer)
    else pointer
  )
}

#' Format settings for the import_settings Rcpp method.
#' @noRd
char_settings_input_ <- function(defs) {
  # change scientific notation locally
  withr::local_options(list(scipen = 500))
  # correct id / house trace vectors
  for (key in intersect(c("id_draws", "house_draws"), names(defs))) {
    y <- unlist(defs[[key]])
    if (!length(y)) defs[[key]] <- ""
    else defs[[key]] <- paste0(unique(y), collapse = ":")
  }
  data.frame(lapply(defs, as.character))
}

#' Assert, then resolve the sampler settings of an inference run:
#' * `steps` -> extended to store the last thinned draw;
#' * sample size -> per chain and in total.
#' Counts are whole numbers, `thin` up to `steps`.
#' @noRd
resolve_settings_ <- function(steps, thin, adapt, chains) {
  # a whole number past the integer range stays a double
  assert_whole <- function(x, what, min = 1) {
    if (rlang::is_scalar_integerish(x) && !is.na(x) && x >= min) return()
    # nolint next: object_usage_linter. read by the cli glue below
    kind <- if (min > 0) "positive" else "non-negative"
    stop_h("Wrong settings: {what} should be a {kind} integer.")
  }
  assert_whole(steps, "number of steps")
  assert_whole(thin, "thinning interval")
  assert_whole(adapt, "number of adaptive steps", min = 0)
  assert_whole(chains, "number of chains")
  if (thin > steps) stop_h(
    "Wrong settings: {.arg thin} cannot be greater than {.arg steps}."
  )
  # round the sample size up, then extend steps to fill it
  chain <- steps %/% thin + (steps %% thin > 0)
  list(steps = chain * thin, chain = chain, total = as.double(chain) * chains)
}

#' Format data, wide or long, for the import Rcpp method.
#' Column names cross verbatim: they are the keys the engine registers,
#' and the definitions sourcing them are written against the user's spelling.
#' @noRd
char_data_input_ <- function(df) {
  if (is.null(df)) NULL
  else data.frame(lapply(df, as.character), check.names = FALSE)
}

#' Combine R6 `objects` into a clone of the first one.
#' Functions compare by arguments, body and enclosing environment.
#' * classes differ -> error;
#' * private field outside `pooled` differs -> error naming it,
#'   down to the entry of a named list;
#' * `secs_` pooled -> run times add up, other pooled fields left to the caller.
#' @param pooled `character()` Private fields the combination pools.
#' @param what `character(1)` Plural noun of the objects, for the error.
#' @return `list` The clone `out`, its private environment `pvt`,
#' and the private environment of each object `privs`.
#' @noRd
combine_r6_ <- function(objects, pooled, what) {
  classes <- unique(vapply(objects, \(x) class(x)[[1L]], character(1L)))
  if (length(classes) > 1L) stop_h(
    "Objects of class{?es} {.cls {classes}} cannot be combined.",
    "Combine {what} of one kind at a time."
  )
  privs <- lapply(objects, \(x) x$.__enclos_env__$private)
  first <- privs[[1L]]
  fields <- setdiff(names(Filter(Negate(is.function), as.list(first))), pooled)
  # name of the field, or of its entries, `x` and `y` differ in
  differ_in <- function(field, x, y) {
    a <- x[[field]]
    b <- y[[field]]
    entries <- is.list(a) && !is.data.frame(a) && !is.null(names(a)) &&
      identical(names(a), names(b))
    if (entries) return(names(a)[!mapply(identical, a, b)])
    if (identical(a, b)) return()
    if (field == "name_") "model" else sub("_$", "", field)
  }
  differ <- unique(unlist(lapply(privs[-1L], function(other) {
    lapply(fields, differ_in, x = first, y = other)
  })))
  if (length(differ)) stop_h(
    "The {what} differ in {.field {differ}}.",
    "Only {what} of the same model, data and settings combine."
  )
  out <- objects[[1L]]$clone()
  pvt <- out$.__enclos_env__$private
  if ("secs_" %in% pooled)
    pvt$secs_ <- sum(vapply(privs, `[[`, numeric(1L), "secs_"))
  list(out = out, pvt = pvt, privs = privs)
}

#' Resolve the character `status =` labels of latent definitions
#' to codes 1..K against the model status dictionary.
#' * integer code -> passes through;
#' * integer and character declarations mixed -> error;
#' * label absent from the dictionary -> error.
#' @param dictionary `character()` Ordered infected-status labels.
#' @return `data.frame` | `NULL` Definitions, `status` coded as integers.
#' @noRd
resolve_status_codes_ <- function(defs, dictionary) {
  if (is.null(defs)) return(defs)
  st <- as.character(unlist(defs[["status"]], use.names = FALSE))
  known <- !is.na(st)
  as_int <- suppressWarnings(as.integer(st))
  is_int <- known & !is.na(as_int) & as.character(as_int) == st
  is_lab <- known & !is_int
  # guard: a model is all-integer or all-character on 'status ='
  if (any(is_int) && any(is_lab)) stop_h(
    "Mixed integer and character {.code status =} declarations in one model.",
    "Use either integer codes or character labels throughout."
  )
  if (any(is_lab)) {
    code <- match(st[is_lab], dictionary)
    if (anyNA(code)) stop_h(paste(
      "Status label(s) {.val {unique(st[is_lab][is.na(code)])}}",
      "not among the declared statuses {.val {dictionary}}."
    ))
    st[is_lab] <- as.character(code)
    defs[["status"]] <- st
  }
  defs
}

#' Resolve unsettled sources against the model's own declarations.
#' The distribution and proposal columns carry encoded arguments too,
#' so they are rewritten alongside the value-bearing ones.
#' @param dict `character()` Source letter, named by declared name.
#' @noRd
resolve_source_types_ <- function(defs, dict) {
  if (!length(defs)) return(defs)
  cols <- c("const", "func", "prior", "model", "min", "max", "proposal")
  for (col in intersect(names(defs), cols)) {
    if (is.character(defs[[col]]))
      defs[[col]] <- housepi_resolve_definition(defs[[col]], dict)
  }
  defs
}

#' Resolve unspecified `discr` declarations against model and data.
#' * definition neither carrying a model/prior nor sampled -> `discr` FALSE.
#' @param discr_req_defs `character()` Definitions the model requires discrete.
#' @param discr_col_types `logical()` Storage type per data column, by name.
#' @param discr_req_cols `character()` Data columns the model requires discrete.
#' @param is_prior `logical(1)` Whether prior, otherwise model.
#' @return `data.frame` Definitions with a logical `discr` column.
#' @noRd
resolve_discr_flags_ <- function(defs, discr_req_defs,
                                 discr_col_types = NULL,
                                 discr_req_cols = character(),
                                 is_prior = FALSE) {
  if (is.null(defs)) return(defs)
  discr <- as.logical(defs[["discr"]])
  unresolved <- is.na(discr)
  if (!any(unresolved)) return(defs)

  # whether the transmission model requires the definition's name discrete
  by_model <- defs[["name"]] %in% discr_req_defs

  # whether the datum column the definition reads is discrete
  by_data <- rep(FALSE, length(discr))
  if (!is.null(discr_col_types)) {
    # name of the data column a definition is pinned to
    col <- as.character(defs[["const"]])
    hit <- grepl("^[A-Za-z0-9_.]+\\|D$", col)
    # NA if it reads none
    col <- ifelse(hit, sub("\\|D$", "", col), NA_character_)
    # the column storage type, or the model requiring that column discrete
    by_data <- unname(discr_col_types[col]) | col %in% discr_req_cols
    by_data[is.na(by_data)] <- FALSE
  }

  distr_col <- if (is_prior) "prior" else "model"
  # sampled: a prior, a model, or bare (no constant, no function)
  has_distr <- !is.na(defs[[distr_col]]) |
    (is.na(defs[["const"]]) & is.na(defs[["func"]]))
  mask <- has_distr & (by_model | by_data)
  discr[unresolved] <- mask[unresolved]
  defs[["discr"]] <- discr
  defs
}

#' Code a character or factor status column into engine status codes.
#' * escaped label -> `0`;
#' * every other label -> `1..K` in level or appearance order,
#'   so a simulated factor round-trips to the same codes;
#' * escaped label unmatched -> every individual counted as infected.
#' @param escaped `character(1)` Label of status code 0.
#' @return `list` Entries `codes`, `integer()`, and `labels`, `character()`.
#' @noRd
status_to_codes_ <- function(x, escaped = "escaped") {
  lv <- if (is.factor(x)) levels(x) else unique(as.character(x))
  lv <- lv[!is.na(lv)]
  # an all-missing column carries no outcome at all, so it stays silent
  if (length(lv) && !escaped %in% lv) {
    guess <- lv[tolower(lv) %in% c("susceptible", "negative", "control",
                                   "no", "false", "uninfected", "healthy", "0")]
    warn_h(
      "No {.field status} value matches the escaped label {.val {escaped}}.",
      `!` = "Every individual counts as infected. Values found: {.val {lv}}.",
      if (length(guess) == 1) "Did you mean {.code escaped = \"{guess}\"}?"
      else paste("Name the value carried by individuals who were not infected",
                 "with {.code escaped =}.")
    )
  }
  labels <- as.character(setdiff(lv, escaped))
  code <- stats::setNames(c(0L, seq_along(labels)), c(escaped, labels))
  list(codes = unname(code[as.character(x)]), labels = labels)
}

#' Factorize the dates and non-numeric columns of `df` for the engine.
#' @param escaped `character(1)` Label of status code 0,
#' in a character or factor `status` column.
#' @noRd
cpp_data_factorize <- function(df, escaped = "escaped") {
  skip <- c("id", "hh_id", "status")
  cols <- colnames(df)
  info <- list()

  # character/factor status: factorize to integer codes (NA -> missing)
  status_labels <- character()
  if ("status" %in% cols && (is.character(df[["status"]]) ||
                               is.factor(df[["status"]]))) {
    coded <- status_to_codes_(df[["status"]], escaped)
    df[["status"]] <- coded$codes
    status_labels <- coded$labels
  }

  # extract all infectious status codes (now integer-coded)
  status <- integer()
  if ("status" %in% cols) {
    s <- suppressWarnings(as.integer(df[["status"]]))
    status <- sort(unique(s[!is.na(s) & s > 0]))
  }

  discr <- stats::setNames(logical(length(cols)), cols)

  for (col in setdiff(cols, skip)) {
    x <- df[[col]]
    discr[[col]] <- inherits(x, "Date") || is.integer(x) || is.logical(x)
    if (inherits(x, "Date")) {
      info[[col]] <- list(class = "Date")
      df[[col]] <- as.numeric(x)
    } else if (is.logical(x)) {
      # FALSE reads 0 and TRUE reads 1
      info[[col]] <- list(class = "logical")
      df[[col]] <- as.integer(x)
    } else if (is.character(x) || is.factor(x)) {
      cls <- if (is.character(x)) "character" else "factor"
      # a missing cell is missing data, not a category: it carries no level
      # and codes to NA, which the engine reads as an empty cell
      y <- factor(x, levels = unique(x))
      info[[col]] <- list(class = cls, levels = levels(y))
      df[[col]] <- as.integer(y)
    }
  }

  list(data = char_data_input_(df), info = info, cols = cols, discr = discr,
       status = status, status_labels = status_labels)
}

#' Restore the factor columns of data returned by the engine.
#' @param info `list` of `list`s Class and levels of the data columns.
#' @param headers `character()` Header columns.
#' @param status_levels `character()` Status levels.
#' @param escaped `character(1)` Label of status code 0.
#' @noRd
cpp_data_refactor <- function(data, info, headers,
                              status_levels = character(),
                              escaped = "escaped") {
  # re-factor
  for (col in names(info)) {
    cls <- info[[col]]$class
    if (identical(cls, "Date")) {
      data[[col]] <- as.Date(data[[col]])
    } else if (identical(cls, "logical")) {
      data[[col]] <- as.logical(data[[col]])
    } else {
      lv <- info[[col]]$levels
      data[[col]] <- factor(data[[col]],
                            levels = seq_along(lv),
                            labels = lv)
      if (identical(cls, "character")) {
        data[[col]] <- as.character(data[[col]])
      }
    }
  }

  if ("status" %in% colnames(data)) {
    s <- as.integer(data[["status"]])
    if (length(status_levels)) {
      # code 0 -> escaped label, code k -> k-th infected label
      lv <- c(escaped, status_levels)
      data[["status"]] <- factor(s + 1L, levels = seq_along(lv), labels = lv)
    } else {
      data[["status"]] <- s
    }
  }
  # re-order data and append any extra columns
  data[c(headers, setdiff(colnames(data), headers))]
}

#' @title Set a warm default theme for housepi visualizations
#' @param base_size `numeric(1)` Base font size, in points.
#' @param base_family `character(1)` Base font family, `""` for the device's.
#' @param transparent `logical(1)`
#' Whether to leave the plot background transparent.
#' @description
#' Set the housepi vignette themes of `ggplot2`, `bayesplot` and `flextable`,
#' for each package installed.
#' @return `NULL`, invisibly.
#' The function sets the themes.
#' @examples
#' \donttest{
#' # apply housepi defaults to ggplot2, bayesplot and flextable
#' set_housepi_theme()
#' }
#' @export
set_housepi_theme <- function(base_size = 11,
                              base_family = "",
                              transparent = FALSE) {
  # ggplot theme
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    my_elem_text <- ggplot2::element_text(size = base_size,
                                          color = "black",
                                          face = "plain")
    my_elem_tag <- ggplot2::element_text(size = base_size + 1,
                                         color = "black",
                                         face = "bold")
    my_panel <- ggplot2::element_rect(fill = "#FEF9EC", color = NA)
    my_strip <- ggplot2::element_rect(fill = "#FBE09B", color = NA)
    my_grid <- ggplot2::element_line(color = "#FDFEFE", linewidth = 1)
    my_plot <- ggplot2::element_rect(
      color = NA, fill = if (isTRUE(transparent)) "transparent" else "white"
    )
    ggplot2::theme_set(
      ggplot2::theme_minimal(base_size = base_size,
                             base_family = base_family) +
        ggplot2::theme(text = my_elem_text,
                       axis.text = my_elem_text,
                       strip.text = my_elem_text,
                       plot.tag = my_elem_tag,
                       panel.background = my_panel,
                       strip.background = my_strip,
                       plot.background = my_plot,
                       legend.position = "bottom",
                       legend.text = my_elem_text,
                       panel.grid.major = my_grid,
                       panel.grid.minor = ggplot2::element_blank(),
                       legend.title = ggplot2::element_blank())
    )

    # mimic element_geom for older ggplot2 versions
    my_text_feats <- list(size = (base_size - 1) / ggplot2::.pt)
    my_point_feats <- list(color = "#FF7F50", size = 2)
    my_line_feats <- list(color = "#FF7F50", size = 1)
    my_errbar_feats <- c(my_line_feats, width = 0.1)
    # labels on a borderless box of the panel colour, hiding the grid lines
    my_label_feats <- c(my_text_feats, fill = my_panel$fill, linewidth = 0)
    ggplot2::update_geom_defaults("text", my_text_feats)
    ggplot2::update_geom_defaults("label", my_label_feats)
    ggplot2::update_geom_defaults("line", my_line_feats)
    ggplot2::update_geom_defaults("point", my_point_feats)
    ggplot2::update_geom_defaults("vline", my_line_feats)
    ggplot2::update_geom_defaults("hline", my_line_feats)
    ggplot2::update_geom_defaults("abline", my_line_feats)
    ggplot2::update_geom_defaults("segment", my_line_feats)
    ggplot2::update_geom_defaults("errorbar", my_errbar_feats)

    # bayesplot theme
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      bayesplot::bayesplot_theme_set(ggplot2::theme_get())
      bayesplot::color_scheme_set(c("#6ec5ff", "#3db1ff", "#0d8fff",
                                    "#0070dd", "#0055aa", "#003d77"))
    }
  }

  # flextable theme
  if (requireNamespace("flextable", quietly = TRUE))
    flextable::set_flextable_defaults(font.family = base_family,
                                      font.size = base_size,
                                      font.color = "black",
                                      background.color = "white",
                                      padding.left = 2,
                                      padding.right = 2,
                                      theme_fun = "theme_vanilla")
  invisible()
}
