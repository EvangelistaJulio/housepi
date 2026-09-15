## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

# Render a `cli_bullets()` block to stdout as a single block.
print_bullets <- function(x, .envir = parent.frame()) {
  cat(cli::cli_fmt(cli::cli_bullets(x, .envir = .envir)), sep = "\n")
}

# Flatten `...`/vectors into one line per element, splitting on newlines.
# Names are preserved; segments produced by a split inherit an empty name.
flatten_cli_parts_ <- function(dots) {
  parts <- unlist(dots, use.names = TRUE)
  nms <- names(parts)
  if (is.null(nms)) nms <- rep("", length(parts))
  out <- character()
  out_nms <- character()
  for (k in seq_along(parts)) {
    segs <- strsplit(parts[[k]], "\n", fixed = TRUE)[[1]]
    out <- c(out, segs)
    out_nms <- c(out_nms, nms[[k]], rep("", max(0, length(segs) - 1)))
  }
  stats::setNames(out, out_nms)
}

#' Emit a harmonized cli error or warning condition.
#' @description
#' First element becomes the header (`x` for errors, unnamed for warnings).
#' Every other element without an explicit bullet name becomes a hint.
#' @param parts `character()` Message lines, optionally named with bullet keys.
#' @param level `character(1)` Either `"stop"` or `"warn"`.
#' @param trusted `logical(1)` When `FALSE`, braces are escaped so foreign text
#'   cannot trip cli glue.
#' @param .envir `environment` Frame used to resolve `{ }` interpolations.
#' @noRd
emit_cli_condition <- function(parts, level = c("stop", "warn"),
                               trusted = TRUE, .envir = parent.frame()) {
  parts <- flatten_cli_parts_(parts)
  level <- match.arg(level)

  # foreign text carries no markup: doubling renders its braces literally
  if (!trusted) parts[] <- gsub("}", "}}", gsub("{", "{{", parts, fixed = TRUE),
                                fixed = TRUE)

  # capitalize the leading letter of each line
  parts[] <- sub("^([a-z])", "\\U\\1", parts, perl = TRUE)

  # un-named continuation lines become info hints
  nms <- names(parts)
  nms[nms == ""] <- "i"
  if (length(nms) > 0) nms[[1]] <- ifelse(level == "stop", "x", "")
  bullets <- stats::setNames(parts, nms)

  if (level == "stop") cli::cli_abort(bullets, call = NULL, .envir = .envir)
  # display at emission time, so warning are not deferred past the errors
  else withr::with_options(list(warn = 1), cli::cli_warn(
    bullets, call = NULL, .envir = .envir
  ))
}

#' Stop execution with a harmonized, cli-rendered error.
#' @param ... Message lines.
#' First is the header; further lines become `i`, unless specified otherwise.
#' @noRd
stop_h <- function(...) {
  emit_cli_condition(list(...), "stop", .envir = parent.frame())
}

#' Emit a harmonized, cli-rendered warning.
#' @param ... Message lines.default.
#' Following the same bullet convention as `stop_h()`.
#' @noRd
warn_h <- function(...) {
  emit_cli_condition(list(...), "warn", .envir = parent.frame())
}

#' Restyle any condition raised beneath `expr` through the cli emitter.
#' @description
#' Condition raised are split on newlines into a header plus hints.
#' C++ conditions are recognised by class and rendered as cli markup.
#' Foreign conditions are shown literally.
#' @param expr Expression to evaluate.
#' @noRd
with_cpp_conditions <- function(expr) {
  reemit <- function(cnd, level, trusted) {
    msg <- conditionMessage(cnd)
    emit_cli_condition(list(msg), level, trusted = trusted)
  }
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      reemit(e, "stop", inherits(e, "housepi::CppError"))
    }),
    warning = function(w) {
      reemit(w, "warn", inherits(w, "housepi_cpp_warning"))
      invokeRestart("muffleWarning")
    }
  )
}

#' Standard error message for read-only fields.
#' @param active `character(1)` Name of the active method.
#' @noRd
stop_read_only <- function(active) {
  stop_h("Active field {.field {active}} is read only.")
}

#' Subsample a data.frame to ~n evenly spaced rows.
#' Returns the input unchanged when n is NA or n >= nrow(df).
#' @param df `data.frame` input table.
#' @param n `integer(1)` Target number of rows.
#' @noRd
helper_thin_table_ <- function(df, n) {
  if (is.null(df) || is.na(n) || nrow(df) <= n) return(df)
  df[round(seq.int(1L, nrow(df), length.out = n)), , drop = FALSE]
}

#' Return 's' if length of input is larger than 1.
#'
#' Thin wrapper over cli pluralization: the marker `{?false/true}` selects
#' `false` for a singular quantity and `true` for a plural one. cli treats a
#' quantity of 0 as plural (unlike the old `len > 1` test), but every call site
#' guards on a non-empty list, so only the 1-vs-many split is ever exercised.
#' @param txt `character(1)` Leading text.
#' @param len_list `vector(1)` Reference list length.
#' @param true `character(1)` What to print if length > 1.
#' @param false `character(1)` What to print otherwise.
#' @noRd
as_plural <- function(txt, len_list, true = "s", false = "") {
  form <- paste0("{cli::qty(", len_list, ")}{?", false, "/", true, "}")
  paste0(txt, vapply(form, cli::pluralize, character(1), USE.NAMES = FALSE))
}

#' Helper for "there were/was  `1`/`n`/`more than n1`" sentence.
#' @param len_list `vector(1)` Reference list length.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @noRd
format_there_were_was <- function(len_list, max_print = 4) {
  paste0(as_plural("There ", len_list, "were ", "was "),
         ifelse(len_list > max_print,
                paste0("more than ", max_print), len_list))
}

#' Collapse a character vector into a phrase.
#' @param mylist `character()` Vector to print.
#' @param quote `logical(1)` Whether to add quotes.
#' @param and `logical(1)` Whether to include a 'and'.
#' @param style `character(1)` cli markup class wrapping each item.
#' @return `character(1)`
#' @noRd
character_to_phrase <- function(mylist, quote = FALSE, and = TRUE,
                                style = "val") {
  if (length(mylist) == 0) return(NULL)
  # wrap in cli markup on request (styled downstream at emit time)
  if (isTRUE(quote)) mylist <- paste0("{.", style, " ", mylist, "}")
  # let cli collapse; " and " before the last item, or plain comma without it
  cli::ansi_collapse(mylist, sep = ", ",
                     last = if (isTRUE(and)) " and " else ", ", trunc = Inf)
}

#' Collapse a character vector into a leading-spaced string of quoted words.
#' The leading space separates the list from the text that introduces it.
#' @param mylist `character()` List to be printed.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @param style `character(1)` cli markup class wrapping each item.
#' @return `character(1)`
#' @noRd
character_to_string_list <- function(mylist, max_print = 4, style = "val") {
  len_list <- length(mylist)
  max_print <- min(max_print, len_list, na.rm = TRUE)
  overflow <- len_list > max_print
  sms <- character_to_phrase(mylist[seq_len(max_print)],
                             quote = TRUE, and = !overflow, style = style)
  paste0(sms, ifelse(overflow, ", ...", "."))
}

#' Find closest match across list of options.
#' Return a string with helpful comments.
#' @param input `character(1)` String to be inspected.
#' @param options `character()` Options.
#' @param min_dist `integer(1)` Threshold for minimal distance.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @return `character()`
#' @keywords internal
#' @noRd
find_closest_match <- function(input, options, min_dist = 3, max_print = 4) {
  output <- c()

  len_opts <- length(options)
  if (len_opts > 0 && min_dist > 0) {
    dists <- utils::adist(input, options)
    idxs <- apply(dists, 1, min) <= min_dist

    if (any(idxs)) {
      closest <- options[apply(dists[idxs, , drop = FALSE], 1, which.min)]
      output <- sprintf("Did you mean {.val %s}?", closest)
    }

    if (!all(idxs) && max_print > 0L) output <- c(output,
      sprintf("Available %s %s",
              as_plural("option", len_opts, "s are:", " is"),
              character_to_string_list(options, max_print = max_print))
    )
  }
  output
}

# Helper to format invalid type errors.
format_type_error <- function(x) {
  sprintf("Got an object of class {.cls %s} instead.",
          paste(rev(class(x)), collapse = "/"))
}

raise_permissions_error <- function() {
  stop_h("Failed to generate necessary files.",
         "Check your user permissions.")
}

# wrapper to file.copy with overwrite option
overwrite <- function(from, to, ...) {
  if (from == to) return(invisible())
  # user permission error (tipically)
  if (!suppressWarnings(file.copy(from, to, overwrite = TRUE, ...))) {
    raise_permissions_error()
  }
}

############################################################################
############################################################################
############################################################################

# return time of modification truncated
trunc_mtime <- function(path) {
  max(trunc(file.info(path)$mtime))
}

#' @title Detect untracked changes in template source files
#' @description
#' Whether source files have been modified or a more recent installation exists.
#' @param source_dir `character(1)` Path to source directory.
#' @param dll `character(1)` Expected dynamic library file name.
#' @return `bool`
#' @keywords internal
#' @noRd
untracked_changes <- function(source_dir, dll) {
  last_mod <- trunc_mtime(file.path(source_dir, "model.cpp"))

  # last changed should be older than installation time
  housepi_dir <- system.file("template", package = "housepi")
  check_src_files <- list.files(housepi_dir, full.names = TRUE)
  install_time <- trunc_mtime(check_src_files)

  any(trunc_mtime(dll) < c(last_mod, install_time))
}

# File shipped by a model for each templated definition family.
template_def_files <- c(parameter = "param.R",
                        latent = "latent.R",
                        `infectivity profile` = "infecprof.R")

# helper function to deep clone R6 objects
helper_clone <- function(x) {
  if (is.data.frame(x)) x[TRUE, , drop = FALSE]
  else if (is.list(x)) lapply(x, helper_clone)
  else if (inherits(x, "R6")) x$clone(deep = TRUE)
  else x
}

############################################################################
############################################################################
############################################################################

#' Assert undocumented advanced setup overrides from `...`.
#' @param allowed `character()` accepted argument names.
#' @param ... `r rd_dots_config_`
#' @return `list` Non-null overrides keyed by engine setup key.
#' @keywords internal
#' @noRd
assert_advanced_setup_ <- function(allowed, ...) {
  dots <- list(...)
  for (nm in names(dots)) {
    if (!(nm %in% allowed)) stop_h(
      "Unexpected argument {.arg {nm}} passed through {.code ...}.",
      # NOTE: do not print all advanced options, just suggest
      find_closest_match(nm, setdiff(allowed, names(dots)), max_print = 0)
    )
  }
  invisible()
}

#' Format elements of a list to pass to import_setup Rcpp method.
#' @param defs `list` Setup options.
#' @return `character()`
#' @keywords internal
#' @noRd
char_setup_input_ <- function(defs) {
  # change scientific notation locally
  withr::local_options(list(scipen = 500))
  # correct logical flags
  for (flag in intersect(c("progress", "follow_default"), names(defs))) {
    y <- unlist(defs[[flag]])
    if (is.logical(y)) defs[[flag]] <- as.integer(y)
  }
  # correct id / house trace vectors
  for (key in intersect(c("id_draws", "house_draws"), names(defs))) {
    y <- unlist(defs[[key]])
    if (length(y) == 0) defs[[key]] <- NA_character_
    else defs[[key]] <- paste0(unique(y), collapse = ":")
  }
  data.frame(lapply(defs, as.character))
}

# Derive (per-chain, total) sample sizes from a setup list, mirroring the
# ceiling-division formula in Setup::assert_mcmc_setup.
setup_sample_size_ <- function(steps, burn, thin, chains) {
  diff <- steps - burn
  c_size <- suppressWarnings(diff %/% thin + (diff %% thin > 0))
  list(chain = c_size, total = c_size * chains)
}

#' Format elements of data (wide or long) to pass to import Rcpp method.
#' Column names cross verbatim: they are the keys the engine registers, and
#' the definitions that source them are written against the user's spelling.
#' @param df `data.frame` Definition table or list.
#' @return `character()`
#' @keywords internal
#' @noRd
char_data_input_ <- function(df) {
  if (is.null(df)) NULL
  else data.frame(lapply(df, as.character), check.names = FALSE)
}

#' Resolve character `status =` labels in a latent defs frame to integer
#' codes (1..K) against the model status dictionary. Integer codes pass through.
#' Guards against mixing integer and character status declarations and against
#' labels absent from the dictionary.
#' @param defs `data.frame` or `NULL`. Latent definitions.
#' @param dictionary `character()`. Ordered infected-status labels.
#' @return `data.frame` with an integer-coded `status` column.
#' @keywords internal
#' @noRd
resolve_status_codes_ <- function(defs, dictionary) {
  if (is.null(defs) || !"status" %in% names(defs)) return(defs)
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

#' Resolve unsettled source sigils against the model's own declarations.
#' The distribution and proposal columns carry encoded arguments too, so they
#' are rewritten alongside the value-bearing ones.
#' @param defs `data.frame` Definitions table.
#' @param dict named `character()` Declared name to source letter.
#' @return `data.frame`
#' @keywords internal
#' @noRd
resolve_source_types_ <- function(defs, dict) {
  if (is.null(defs) || length(defs) == 0) return(defs)
  cols <- c("value", "func", "prior", "model", "min", "max", "proposal")
  for (col in intersect(names(defs), cols)) {
    if (is.character(defs[[col]])) {
      defs[[col]] <- housepi_resolve_definition(defs[[col]], dict)
    }
  }
  defs
}

#' Resolve unspecified `discr` declarations against model and data.
#' A definition without model/prior keeps discr = FALSE.
#' @param defs `data.frame` Definitions table.
#' @param discr_req_defs `character()` Definitions the model requires discrete.
#' @param discr_col_types named `logical()` Storage type per data column.
#' @param discr_req_cols `character()` Data columns the model requires discrete.
#' @param is_prior `logical(1)` Whether prior, otherwise model.
#' @return `data.frame` with a logical `discr` column.
#' @keywords internal
#' @noRd
resolve_discr_flags_ <- function(defs, discr_req_defs,
                                 discr_col_types = NULL,
                                 discr_req_cols = character(),
                                 is_prior = FALSE) {
  if (is.null(defs) || !"discr" %in% names(defs)) return(defs)
  discr <- as.logical(defs[["discr"]])
  unresolved <- is.na(discr)
  if (!any(unresolved)) return(defs)

  # whether the transmission model requires the definition's name discrete
  by_model <- defs[["name"]] %in% discr_req_defs

  # whether the datum column the definition reads is discrete
  by_data <- rep(FALSE, length(discr))
  if (!is.null(discr_col_types)) {
    # name of the data column a definition value reads
    col <- as.character(defs[["value"]])
    hit <- grepl("^[A-Za-z0-9_.]+\\|D$", col)
    # NA if it reads none
    col <- ifelse(hit, sub("\\|D$", "", col), NA_character_)
    # the column storage type, or the model requiring that column discrete
    by_data <- unname(discr_col_types[col]) | col %in% discr_req_cols
    by_data[is.na(by_data)] <- FALSE
  }

  distr_col <- ifelse(is_prior, "prior", "model")
  mask <- !is.na(defs[[distr_col]]) & (by_model | by_data)
  discr[unresolved] <- mask[unresolved]
  defs[["discr"]] <- discr
  defs
}

#' Code a character/factor status column into engine status codes.
#'
#' The escaped label maps to `0` and every other label to `1..K` in
#' level/appearance order, so a simulated factor round-trips to the same codes.
#' An unmatched escaped label leaves every individual counted as infected.
#' @param x `character()` | `factor()` Status column.
#' @param escaped `character(1)` Label for status code 0.
#' @return `list` with `codes` (`integer()`) and `labels` (`character()`).
#' @keywords internal
#' @noRd
status_to_codes_ <- function(x, escaped = "escaped") {
  lv <- if (is.factor(x)) levels(x) else unique(as.character(x))
  lv <- lv[!is.na(lv)]
  # an all-missing column carries no outcome at all, so it stays silent
  if (length(lv) > 0 && !escaped %in% lv) {
    guess <- lv[tolower(lv) %in% c("susceptible", "negative", "control",
                                   "no", "false", "uninfected", "healthy", "0")]
    warn_h(
      "No {.field status} value matches the escaped label {.val {escaped}}.",
      "Every individual counts as infected. Values found: {.val {lv}}.",
      if (length(guess) == 1) "Did you mean {.code escaped = \"{guess}\"}?"
      else paste("Name the value carried by individuals who were not infected",
                 "with {.code escaped =}.")
    )
  }
  labels <- as.character(setdiff(lv, escaped))
  code <- stats::setNames(c(0L, seq_along(labels)), c(escaped, labels))
  list(codes = unname(code[as.character(x)]), labels = labels)
}

#' Factorize dates and non-numerical quantities from data to pass to C++.
#' @param df `data.frame` Dataset.
#' @param escaped `character(1)` Label for status code 0 in character/factor
#' `status` columns.
#' @return `list`
#' @keywords internal
#' @noRd
cpp_data_factorize <- function(df, escaped = "escaped") {
  skip <- c("id", "house_id", "status")
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
    discr[[col]] <- inherits(x, "Date") || is.integer(x)
    if (inherits(x, "Date")) {
      info[[col]] <- list(class = "Date")
      df[[col]] <- as.numeric(x)
    } else if (is.character(x) || is.factor(x) || is.logical(x)) {
      cls <- ifelse(is.character(x), "character",
                    ifelse(is.logical(x), "logical", "factor"))
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

#' Restore factored quantities on data returned from the C++ implementation.
#' @param data `data.frame` Output of the Rcpp model object.
#' @param info `list` of `list`s. Info about class and levels of data columns.
#' @param headers `character()` Header columns.
#' @param status_levels `character()` Status levels.
#' @param escaped `character(1)` Label for status code 0.
#' @return `data.frame`
#' @keywords internal
#' @noRd
cpp_data_refactor <- function(data, info, headers,
                              status_levels = character(),
                              escaped = "escaped") {
  # re-factor
  for (col in names(info)) {
    cls <- info[[col]]$class
    if (identical(cls, "Date")) {
      data[[col]] <- as.Date(data[[col]])
    } else {
      lv <- info[[col]]$levels
      data[[col]] <- factor(data[[col]],
                            levels = seq_along(lv),
                            labels = lv)
      if (identical(cls, "character")) {
        data[[col]] <- as.character(data[[col]])
      } else if (identical(cls, "logical")) {
        data[[col]] <- as.logical(data[[col]])
      }
    }
  }

  if ("status" %in% colnames(data)) {
    s <- as.integer(data[["status"]])
    if (length(status_levels) > 0) {
      # code 0 -> escaped label, code k -> k-th infected label
      lv <- c(escaped, status_levels)
      data[["status"]] <- factor(s + 1L, levels = seq_along(lv), labels = lv)
    } else {
      data[["status"]] <- s
    }
  }
  # re-order data and append eventual extra columns
  data[c(headers, setdiff(colnames(data), headers))]
}

#' Evaluate the symbolic chain-binomial seropositive-contribution.
#'
#' @param community_escape `numeric` Community escape probabilities.
#' One entry for each of the $n_\mathrm{infec}$ seropositives.
#' @param p2p_escape `matrix` Person-to-person escape-probabilities.
#' Square matrix $n_\mathrm{infec} \times n_\mathrm{infec}$.
#' Convention `[i][k]` of infectee `k` from infector `i`.
#' @return `numeric(1)` Positive-contribution probability.
#' @keywords internal
#' @noRd
eval_chain_binomial_debugger_ <- function(community_escape, p2p_escape) {
  n_infec <- length(community_escape)
  stopifnot(is.matrix(p2p_escape),
            identical(dim(p2p_escape), c(n_infec, n_infec)))

  formula <- .engine$aux_chain_binomial_debugger(n_infec)

  for (r in seq_len(n_infec)) for (c in seq_len(n_infec)) {
    formula <- gsub(paste0("\\bQ", r, "_", c, "\\b"),
                    format(p2p_escape[r, c], digits = 17), formula)
  }
  for (i in seq_len(n_infec)) {
    formula <- gsub(paste0("\\bB", i, "\\b"),
                    format(community_escape[i], digits = 17), formula)
  }
  eval(str2lang(formula))
}

#' @title Set a warm default theme for housepi visualizations
#' @param base_size Base font size (pts). Default `11`.
#' @param base_family Base font family. Default `""` (device default).
#' @param transparent `logical(1)` Leave the plot background transparent.
#' @description
#' If suggested packages are installed, sets corresponding `ggplot2`,
#' `bayesplot` and `flextable` themes as used in housepi vignettes.
#' @examples
#' # apply housepi defaults to ggplot2, bayesplot and flextable
#' housepi_theme_set()
#' @export
housepi_theme_set <- function(base_size = 11,
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
      color = NA, fill = ifelse(isTRUE(transparent), "transparent", "white")
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
    ggplot2::update_geom_defaults("text", my_text_feats)
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
  if (requireNamespace("flextable", quietly = TRUE)) {
    flextable::set_flextable_defaults(font.family = base_family,
                                      font.size = base_size,
                                      font.color = "black",
                                      background.color = "white",
                                      theme_fun = "theme_vanilla")
  }
}
