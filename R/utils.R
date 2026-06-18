## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' Stop execution with a multi-line error message.
#' @param ... Messages printed in sequence on different lines.
#' @noRd
stop_h <- function(...) {
  # FIXME: let's just expose rlang
  rlang::abort(paste(unlist(list(...)), collapse = "\n"),
               call = NULL)
}

#' Emit a multi-line warning message.
#' @param ... Messages printed in sequence on different lines.
#' @noRd
warn_h <- function(...) {
  # FIXME: let's just expose rlang
  rlang::warn(paste(unlist(list(...)), collapse = "\n"),
              call = NULL)
}

#' Standard error message for read-only fields.
#' @param active `character(1)` Name of the active method.
#' @noRd
stop_read_only <- function(active) {
  stop_h(sprintf("'%s' is read only.", active))
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
#' @param txt `character(1)` Leading text.
#' @param len_list `vector(1)` Reference list length.
#' @param true `character(1)` What to print if length > 1.
#' @param false `character(1)` What to print otherwise.
#' @noRd
as_plural <- function(txt, len_list, true = "s", false = "") {
  paste0(txt, ifelse(len_list > 1, true, false))
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
#' @return `character(1)`
#' @noRd
character_to_phrase <- function(mylist, quote = FALSE, and = TRUE) {
  len_list <- length(mylist)
  if (len_list == 0) return(NULL)
  if (len_list == 1) and <- FALSE

  # eventually add quotes
  if (isTRUE(quote)) mylist <- paste0("'", mylist, "'")
  # collapse w\o 'and'
  if (isFALSE(and)) return(paste(mylist, collapse = ", "))
  # collapse with 'and'
  sms <- paste(mylist[1:(len_list - 1)], collapse = ", ")
  paste(sms, "and", mylist[len_list])
}

#' Collapse a character vector into string of quoted words.
#' @param mylist `character()` List to be printed.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @param fmt `logical(1)` Format as in message.
#' @return `character(1)`
#' @noRd
character_to_string_list <- function(mylist, max_print = 4, fmt = TRUE) {
  len_list <- length(mylist)
  max_print <- min(max_print, len_list, na.rm = TRUE)
  overflow <- len_list > max_print
  sms <- character_to_phrase(mylist[seq_len(max_print)],
                             quote = TRUE, and = !overflow)
  out <- paste0(sms, ifelse(overflow, ", ...", "."))
  if (isTRUE(fmt)) out <- paste0(as_plural("", len_list, "\n", " "), out)
  out
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
      output <- sprintf("Did you mean '%s'?", closest)
    }

    if (!all(idxs)) output <- c(output,
      sprintf("Available %s %s",
              as_plural("option", len_opts, "s are:", " is"),
              character_to_string_list(options, fmt = FALSE))
    )
  }
  output
}

# Helper to format invalid type errors.
format_type_error <- function(x) {
  sprintf("Got an object of class %s instead.",
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
  if (!file.copy(from, to, overwrite = TRUE, ...)) raise_permissions_error()
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

# Optional helper files
helper_files <- function(path) {
  helper_names <- c(paste0("template_", c("augm.R", "param.R")), "info.Rd")
  file.path(path, helper_names)
}

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

#' Format elements of a list to pass to import_setup Rcpp method.
#' @param defs `list` Setup options.
#' @return `character()`
#' @keywords internal
#' @noRd
char_setup_input_ <- function(defs) {
  # change scientific notation locally
  withr::local_options(list(scipen = 500))
  # correct progress
  if ("progress" %in% names(defs)) {
    y <- unlist(defs$progress)
    if (is.logical(y)) defs$progress <- as.integer(y)
  }
  # correct id trace vector
  if ("id_draws" %in% names(defs)) {
    y <- unlist(defs$id_draws)
    if (length(y) == 0) defs$id_draws <- NA_character_
    else defs$id_draws <- paste0(unique(y), collapse = ":")
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
#' @param df `data.frame` Definition table or list.
#' @return `character()`
#' @keywords internal
#' @noRd
char_data_input_ <- function(df) {
  if (is.null(df)) NULL
  else data.frame(lapply(df, as.character))
}

#' Factorize dates and non-numerical quantities from data to pass to C++.
#' @param df `data.frame` Dataset.
#' @return `list`
#' @keywords internal
#' @noRd
cpp_data_factorize <- function(df) {
  skip <- c("id", "house_id", "status")
  cols <- colnames(df)
  info <- list()
  # extract all infectious status
  status <- integer()
  if ("status" %in% cols) {
    status <- unique(df[["status"]])
    status <- sort(status[!is.na(status) & status > 0])
  }

  for (col in setdiff(cols, skip)) {
    x <- df[[col]]
    if (inherits(x, "Date")) {
      info[[col]] <- list(class = "Date")
      df[[col]] <- as.numeric(x)
    } else if (is.character(x) || is.factor(x) || is.logical(x)) {
      cls <- ifelse(is.character(x), "character",
                    ifelse(is.logical(x), "logical", "factor"))
      y <- factor(x, levels = unique(x), exclude = NULL)
      info[[col]] <- list(class = cls, levels = levels(y))
      df[[col]] <- as.integer(y)
    }
  }

  list(data = char_data_input_(df), info = info, cols = cols, status = status)
}

#' Restore factored quantities on data returned from the C++ implementation.
#' @param data `data.frame` Output of the Rcpp model object.
#' @param info `list` of `list`s. Info about class and levels of data columns.
#' @param headers `character()` Header columns.
#' @return `data.frame`
#' @keywords internal
#' @noRd
cpp_data_refactor <- function(data, info, headers) {
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
    data[["status"]] <- replace(data[["status"]], data[["status"]] < 0, NA)
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

  formula <- .engine$auxiliary_chain_binomial_debugger(n_infec)

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
#' @description
#' If suggested packages are installed, sets corresponding `ggplot2`,
#' `bayesplot` and `flextable` themes as used in housepi vignettes.
#' @examples
#' # apply housepi defaults to ggplot2, bayesplot and flextable
#' housepi_theme_set()
#' @export
housepi_theme_set <- function(base_size = 11, base_family = "") {
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
    ggplot2::theme_set(
      ggplot2::theme_minimal(base_size = base_size,
                             base_family = base_family) +
        ggplot2::theme(text = my_elem_text,
                       axis.text = my_elem_text,
                       strip.text = my_elem_text,
                       plot.tag = my_elem_tag,
                       panel.background = my_panel,
                       strip.background = my_strip,
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
