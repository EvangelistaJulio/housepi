#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# utils.R
#

#' @title
#' Housepi stop utility
#' @description
#' Utility functions for error messages.
#'
#' @param ...
#' Lists of messages printed in sequence on different lines.
#'
#' @noRd
stop_h <- function(...) {
  base::stop(paste(unlist(list(...)), collapse = "\n"),
             call. = FALSE)
}

#' @title
#' Housepi warning utility
#' @description
#' Utility functions for warning messages.
#'
#' @param ...
#' Lists of messages printed in sequence on different lines.
#'
#' @noRd
warning_h <- function(...) {
  base::warning(paste(unlist(list(...)), collapse = "\n"),
                call. = FALSE)
}

#' @title
#' Character to string list
#' @description
#' Collapse a character vector into string of quoted words.
#'
#' @param mylist
#' A `character()` vector to be print.
#' @param max_print
#' An `integer(1)` with the maximum number of elements to print.
#'
#' @return
#' String formatted as a list of quoted words.
#'
#' @noRd
character_to_string_list <- function(mylist, max_print = 4) {
  max_print <- min(max_print, length(mylist), na.rm = TRUE)
  sms <- paste(lapply(mylist[seq_len(max_print)],
                      function(x) paste0("'", x, "'")),
               collapse = ", ")
  if (length(mylist) > max_print) paste0(sms, ", ...")
  else paste0(sms, ".")
}

#' @title
#' Character to phrase
#' @description
#' Collapse a character vector into a phrase.
#'
#' @param mylist
#' A `character()` vector.
#'
#' @return
#' String formatted as a phrase.
#'
#' @noRd
character_to_phrase <- function(mylist) {
  len_list <- length(mylist)
  if (len_list == 0) {
    return(NULL)
  } else if (len_list == 1) {
    sms <- mylist
  } else {
    sms <- paste(mylist[1:(len_list - 1)], collapse = ", ")
    sms <- paste(sms, "and", mylist[len_list])
  }
  sms
}

#' @title
#' Find closest match
#' @description
#' Find the closest match among a list of options.
#'
#' @param input
#' A `character(1)` with the string to be inspected.
#' @param options
#' A `character()` with the options.
#' @param min_dist
#' An `integer(1)` with the threshold for minimal distance.
#'
#' @return `character()`\cr
#' String with helpful comments
#'
#' @keywords internal
#' @noRd
find_closest_match <- function(input, options, min_dist = 1) {
  output <- c()

  if ((length(options) > 0) && (min_dist > 0)) {
    dists <- utils::adist(input, options)
    idxs <- apply(dists, 1, min) <= min_dist

    if (any(idxs)) {
      closest <- options[apply(dists[idxs, , drop = FALSE], 1, which.min)]
      output <- sprintf("Did you mean '%s'?", closest)
    }

    if (!all(idxs)) output <- c(output,
      sprintf("Available option%s %s",
              ifelse(length(options) > 1, "s are:", " is"),
              character_to_string_list(options))
    )
  }
  output
}

# wrapper to file.copy with overwrite option
overwrite <- function(from, to, ...) {
  if (from != to) {
    # user permission error (tipically)
    if (!file.copy(from, to, overwrite = TRUE, ...)) {
      print(from) # FLAG
      print(to)
      stop_h("Failed to setup necessary files.",
             "Check your user permissions.")
    }
  }
}

# Define a custom get function (like in Python dictionaries)
get_val_ <- function(lst, key, default = NA) {
  if (key %in% names(lst)) lst[[key]]
  else default
}

############################################################################
############################################################################
############################################################################


# Generate a random string
rnd_string <- function(size = 8) {
  paste0(sample(c(letters, LETTERS, 0:9), size, replace = TRUE),
         collapse = "")
}

# Source code files
source_files <- function(path) {
  file.path(path, paste0("model.", c("cpp", "h")))
}

# return time of modification truncated
trunc_mtime <- function(path) {
  max(trunc(file.info(path)$mtime))
}

#' @title
#' Untracjed Changes
#' @description
#' Tell if source files was modified or if a more recent installation exists.
#'
#' @param source_dir
#' A `character(1)` with the path to source directory.
#' @param dll
#' A `character(1)` with the expected dynamic library file name.
#'
#' @return
#' `bool`
#'
#' @keywords internal
#' @noRd
untracked_changes <- function(source_dir, dll) {
  last_mod <- trunc_mtime(source_files(source_dir))

  # last changed should be older than installation time
  housepi_dir <- system.file("template", package = "housepi")
  check_src_files <- list.files(housepi_dir, full.names = TRUE)
  install_time <- trunc_mtime(check_src_files)

  any(trunc_mtime(dll) < c(last_mod, install_time))
}

# Optional helper files
helper_files <- function(path) {
  helper_names <- paste0("template_", c("augm.R", "param.R", "info.Rd"))
  file.path(path, helper_names)
}

# Associate dynamic library extension to name
dll_file_ext <- function(name) {
  paste0(name, .Platform$dynlib.ext)
}

# Templated name convention
template_name <- function(path) {
  # remove punctuation from string
  name <- gsub("[[:punct:]]", "", basename(path))
  paste0("housepi", name)
}

# Add random string to templated name
random_name <- function(path) {
  paste0(template_name(path), rnd_string())
}

# DLL fyle with templated name convention
template_dll_file <- function(path) {
  file.path(path, dll_file_ext(template_name(path)))
}

# flags used to set the compilation environment
get_cpp_environment_flags <- function(paths = NULL,
                                      flags = NULL,
                                      linkers = NULL) {
  # path to headers
  template_path <- system.file("template/include", package = "housepi")
  include_paths <- paste0("-I", c(template_path, paths), collapse = " ")
  include_flags <- paste0("-D", c("_HOUSEPI_RCPP", flags), collapse = " ")
  include_linkers <- paste0(linkers, collapse = " ")

  # C++ flags for housepi Rcpp and OpenMP
  # + linker flags for OpenMP
  c(PKG_CXXFLAGS = sprintf("$(SHLIB_OPENMP_CXXFLAGS) %s %s",
                           include_flags, include_paths),
    PKG_LIBS = sprintf("$(SHLIB_OPENMP_CXXFLAGS) %s",
                       include_linkers))
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

#' @title
#' Character Setup Input
#' @description
#' Format elements of a list to pass to import_setup Rcpp method.
#'
#' @param defs
#' A `list` of setup command
#'
#' @return `character()`
#'
#' @keywords internal
#' @noRd
char_setup_input_ <- function(defs) {
  # change scientific notation locally
  withr::local_options(list(scipen = 500))
  # correct id trace vector
  if ("id_draws" %in% names(defs)) {
    y <- unlist(defs$id_draws)
    if (length(y) == 0) defs$id_draws <- NA_character_
    else defs$id_draws <- paste0(unique(y), collapse = ":")
  }
  data.frame(lapply(defs, as.character))
}

#' @title
#' Character Data Input
#' @description
#' Format elements of data (wide or long) to pass to import Rcpp method.
#'
#' @param df
#' A `data.frame` with definition table or list.
#'
#' @return `character()`
#'
#' @keywords internal
#' @noRd
char_data_input_ <- function(df) {
  if (is.null(df)) NULL
  else data.frame(lapply(df, as.character))
}

#' @title
#' Cpp Data Factorize
#' @description
#' Factorize dates and non numerical quanties from data to pass Cpp.
#'
#' @param df
#' A `data.frame` with dataset.
#'
#' @return `list`
#'
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

#' @title
#' Cpp Data Refactor
#' @description
#' Restore factored quanties on data return from Cpp implementation.
#'
#' @param data
#' A `data.frame` as output from the Rcpp model object.
#' @param info
#' A `list` of `list`s with info about class and levels of data columns.
#' @param headers
#' A `character()` with the column headers.
#'
#' @return `data.frame`
#'
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
  # re-order data and append eventual extra columns
  data[c(headers, setdiff(colnames(data), headers))]
}

#' @title Set the theme for housepi visualizations
#' @description
#' If suggested packages are installad, sets the theme used
#' in housepi vignettes for `ggplot2`, `bayesplot` and `flextable`.
#' @export
housepi_theme_set <- function() {
  # ggplot theme
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    my_elem_text <- ggplot2::element_text(size = 12, face = "plain")
    ggplot2::theme_set(
      ggplot2::theme_gray() +
        ggplot2::theme(text = my_elem_text,
                       axis.text = my_elem_text,
                       strip.text = my_elem_text,
                       legend.text = my_elem_text,
                       legend.title = ggplot2::element_blank())
    )
    # mimic element_geom for older ggplot2 versions
    ggplot2::update_geom_defaults("line", list(color = "#FF7F50", size = 1))
    ggplot2::update_geom_defaults("point", list(color = "#FF7F50", size = 3))
    ggplot2::update_geom_defaults("vline", list(color = "#FF7F50", size = 1))
    ggplot2::update_geom_defaults("hline", list(color = "#FF7F50", size = 1))
    ggplot2::update_geom_defaults("abline", list(color = "#FF7F50", size = 1))
    ggplot2::update_geom_defaults("segment", list(color = "#FF7F50", size = 1))
    ggplot2::update_geom_defaults("errorbar", list(color = "#FF7F50", size = 1,
                                                   width = 0.1))

    # bayesplot theme
    if (requireNamespace("bayesplot", quietly = TRUE)) {
      bayesplot::bayesplot_theme_set(ggplot2::theme_get())
      bayesplot::color_scheme_set(c("#6ec5ff", "#3db1ff",
                                    "#0d8fff", "#0070dd",
                                    "#0055aa", "#003d77"))
    }
  }

  # flextable theme
  if (requireNamespace("flextable", quietly = TRUE)) {
    flextable::set_flextable_defaults(font.color = "#111111",
                                      background.color = "white",
                                      theme_fun = "theme_vanilla")
  }
}
