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
