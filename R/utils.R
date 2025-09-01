#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# utils.R
#

#' @title
#' Housepi stop utility
#' @description
#' Utility functions for error messages.
#'
#' @param ...\cr
#' Lists of messages printed in sequence on different lines.
#'
#' @noRd
stop_h <- function(...) {
  base::stop(paste(unlist(list(...)), collapse = "\n"),
             call. = FALSE)
  invisible()
}

#' @title
#' Housepi warning utility
#' @description
#' Utility functions for warning messages.
#'
#' @param ...\cr
#' Lists of messages printed in sequence on different lines.
#'
#' @noRd
warning_h <- function(...) {
  base::warning(paste(unlist(list(...)), collapse = "\n"),
                call. = FALSE)
  invisible()
}

# check if suggested library is installed
require_suggested_library <- function(libname) {
  if (!requireNamespace(libname, quietly = TRUE)) stop_h(
    sprintf("Attempt to call missing library \"%s\".", libname),
    sprintf("Install it running: install.packages(\"%s\")", libname)
  )
}

#' @title
#' Character to string list
#' @description
#' Collapse a character vector into string of quoted words.
#'
#' @param mylist`character()`\cr
#' Vector of characters.
#' @param max_print `integer(1)`\cr
#' Maximum number of elements to print (default is 4).
#'
#' @return
#' String formatted as a list of quoted words.
#'
#' @noRd
character_to_string_list <- function(mylist, max_print = 4) {
  sms <- paste(lapply(mylist[seq_len(min(max_print, length(mylist)))],
                      function(x) paste0("\"", x, "\"")),
               collapse = ", ")
  if (length(mylist) > max_print) return(paste0(sms, ", ..."))
  else return(paste0(sms, "."))
}

#' @title
#' Character to phrase
#' @description
#' Collapse a character vector into a phrase.
#'
#' @param mylist`character()`\cr
#' Vector of characters.
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
  paste0(sms, ".")
}

#' @title
#' Find closest match
#' @description
#' Find the closest match among a list of options.
#'
#' @param input `character(1)`\cr
#' Vector of characters.
#' @param options `character()`\cr
#' Vector of options.
#' @param min_dist `integer(1)`\cr
#' Threshold for minimal distance.
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
      output <- sprintf("Did you mean \"%s\"?", closest)
    }

    if (!all(idxs)) output <- c(output,
      sprintf("Available option%s %s",
              ifelse(length(options) > 1, "s are:", " is"),
              character_to_string_list(options))
    )

    return(output)
  }
}

# wrapper to file.copy with overwrite option
overwrite <- function(from, to, ...) {
  if (from != to) {
    # user permission error (tipically)
    if (!file.copy(from, to, overwrite = TRUE, ...)) stop_h(
      "Failed to setup necessary files.",
      "Check your user permissions."
    )
  }
}

# Define a custom get function (like in Python dictionaries)
get_val_ <- function(lst, key, default = NA) {
  if (key %in% names(lst)) return(lst[[key]])
  else return(default)
}

#' @title
#' Load a List of Definitions
#' @description
#' Load definitions either as direct arguments or from predefined list.
#'
#' @param ...\cr
#' Named arguments with respective list of definition details.
#' @param from named `list`\cr
#' Named list of definitions (alternative to named arguments).
#'
#' @return
#' `list` of `list`
#'
#' @keywords internal
#' @noRd
load_definition_list <- function(..., from = NULL) {
  definitions <- list(...)

  if (length(definitions) == 0) {
    if (length(from) == 0) stop_h("No definitions provided.")
    else if (is.list(from)) definitions <- from
    else stop_h("Argument 'from' should be a list.")
  } else if (any(names(definitions) == "")) {
    stop_h("Found definition without name.")
  }

  return(definitions)
}

############################################################################
############################################################################
############################################################################

# remove punctuation from string
format_string <- function(name) {
  gsub("[[:punct:]]", "", name)
}

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
#' @param source_dir `character(1)`\cr
#' Path to source directory.
#' @param dll `character(1)`\cr
#' Expected dynamic library file name.
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
  helper_names <- paste0("template_",
                         c("param.R", "info.Rd", "augm.R"))
  file.path(path, helper_names)
}

# Associate dynamic library extension to name
dll_file_ext <- function(name = "") {
  return(paste0(name, .Platform$dynlib.ext))
}

# Templated name convention
template_name <- function(path) {
  paste0("housepi", format_string(basename(path)))
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
#' Default types
#' @description
#' Check that types are set correctly.
#' Format NA types with defaults according to values.
#'
#' @param types `character()`\cr
#' Vector of characters.
#' @param values `list`\cr
#' List of values.
#'
#' @return `character()`\cr
#' Vector of characters.
#'
#' @keywords internal
#' @noRd
default_types_ <- function(types, values) {
  valid_types <- list("S" = "scalar",
                      "P" = "parameter",
                      "D" = "dataset",
                      "A" = "augmented")

  for (i in seq_along(types)) {
    k <- types[i]
    y <- values[[i]]

    if (is.na(k)) {
      # no explicit type provided, then format for default types
      if (is.na(y)) types[i] <- NA
      else if (is.numeric(y)) types[i] <- "S"
      else if (is.character(y)) types[i] <- ifelse(
        is.na(suppressWarnings(as.numeric(y))), "A", "S"
      )
      else types[i] <- "A"

    } else {
      # explicit type checking
      k <- toupper(k)
      if (!(k %in% names(valid_types))) stop_h(
        paste0("Invalid type \"", k, "\" for value \"", y, "\"."),
        paste(
          sprintf("\"%s\"", names(valid_types)),
          sprintf("(%s)", unlist(valid_types)),
          collapse = ", "
        )
      )
    }
  }

  types
}

# @title
#' Type split
#' @description
#' Split a typed string or character vector
#'
#' @param lst `character()`\cr
#' Vector of characters.
#'
#' @return `list`\cr
#' List with two elements: value and type.
#'
#' @keywords internal
#' @noRd
type_split_ <- function(lst) {
  if (all(is.na(lst))) {
    list(val = NA, at = NA)

  } else {
    splitted_lst <- strsplit(as.character(lst), "\\@")
    values <- sapply(splitted_lst, function(x) x[1])
    types <- sapply(splitted_lst, function(x) x[2])
    list(
      val = values,
      at = default_types_(types, values)
    )
  }
}

#' @title
#' Type join
#' @description
#' Join elements of a dictionary to pass to csv
#'
#' @param x `character()`\cr
#' Vector of characters.
#'
#' @return `character(1)`\cr
#' String formatted as a list of quoted words.
#'
#' @keywords internal
#' @noRd
type_join_ <- function(x) {
  paste(x, collapse = "/")
}

#' @title
#' Character Setup Input
#' @description
#' Format elements of a list to pass to import_setup Rcpp method.
#'
#' @param defs `list`\cr
#' List of setup command
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
    if (length(y) == 0) defs$id_draws <- NA
    else defs$id_draws <- type_join_(unique(y))
  }
  data.frame(lapply(defs, as.character))
}

#' @title
#' Character Data Input
#' @description
#' Format elements of dataset (or datalist) to pass to import Rcpp method.
#'
#' @param df `data.frame()`\cr
#' Dataest or datalist.
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
#' Factorize dates and non numerical quanties from dataset to pass Cpp.
#'
#' @param df `data.frame`\cr
#' Dataset.
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
      info[[col]] <- list(class = cls,
                          levels = levels(y))
      df[[col]] <- as.integer(y)
    }
  }
  list(data = char_data_input_(df),
       info = info,
       cols = cols,
       status = status)
}

#' @title
#' Cpp Data Refactor
#' @description
#' Restore factored quanties on dataset return from Cpp implementation.
#'
#' @param data `data.frame`\cr
#' Data as output from rcpp model object.
#' @param info `list` of `list\cr
#' Infor about class and levels of dataset columns.
#' @param headers `character()`\cr
#' Headers of the columns.
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

############################################################################
############################################################################
############################################################################

#' @title
#' Count unique function calls
#' @description
#' Count how many base functions are invoked to avoid double calls
#'
#' @param s `character()`\cr
#' Vector of characters.
#'
#' @return `list`
#'
#' @keywords internal
#' @noRd
count_unique_function_calls_ <- function(s) {
  forbid <- c("-", "+", "=", "*", "/", "^")
  reg_expr <- paste0(c("[", forbid, "]"), collapse = "")
  sum(forbid %in% unlist(regmatches(s, gregexpr(reg_expr, s))))
}

#' @title
#' Extract augmented function equality
#' @description
#' Extract equality augmented function declaration
#'
#' @param s `character()`\cr
#' Vector of characters.
#'
#' @return `list`
#'
#' @keywords internal
#' @noRd
extract_func_equality_ <- function(s) {
  # split according to the regular expression
  y <- strsplit(s, "=")[[1]]
  if (length(y) != 2) stop_h(
    "Invalid number of arguments for equality.",
    "Expected 1 argument, found ", length(y), "."
  )
  y <- type_split_(y[-1])
  list(
    func = "eq", key = "x",
    val = type_join_(y$val), at = type_join_(y$at)
  )
}

#' @title
#' Extract augmented function basic
#' @description
#' Extract basic augmented function declaration
#'
#' @param s `character(1)`\cr
#' Vector of characters.
#' @param reg_expr `character(1)`\cr
#' Regular expression to split the string.
#' @param func_name `character(1)`\cr
#' Function name.
#'
#' @return `list`
#'
#' @keywords internal
#' @noRd
extract_func_arithmetic_ <- function(s, reg_expr, func_name) {
  # split according to the regular expression
  y <- strsplit(s, reg_expr)[[1]]
  if (length(y) != 2) stop_h(
    "Invalid number of arguments for ", func_name, ".",
    "Expected 2 argument, found ", length(y), "."
  )
  y <- type_split_(y)
  list(
    func = func_name, key = type_join_(c("x", "x2")),
    val = type_join_(y$val), at = type_join_(y$at)
  )
}

#' @title
#' Extract function name and argument from a function call string
#'
#' @description
#' Parses a string and returns a list with the function name and argument.
#' Handles common mathematical functions and strips whitespace.
#' Raises an error if parentheses are unmatched.
#'
#' @param s `character(1)`\cr
#' String representing a function call.
#' @return
#' `list`
#'
#' @keywords internal
#' @noRd
extract_func_call_ <- function(s) {
  options <- c("sqrt", "exp", "log", "sin", "cos", "tan")

  if (count_unique_function_calls_(s) > 0) stop_h(
    sprintf("Multiple functions invoked in definition \"%s\"", s)
  )

  # check parentheses are matched
  if (nchar(gsub("[^\\(]", "", s)) != nchar(gsub("[^\\)]", "", s))) stop(
    "Unmatched parentheses in input string"
  )

  # extract function name and argument using regex
  m <- regexec("^([a-zA-Z_][a-zA-Z0-9_]*)\\((.*)\\)$", s)
  parts <- regmatches(s, m)[[1]]

  if (length(parts) != 3) stop_h(
    "Input string must be a function call like func(x)"
  )

  func_name <- parts[2]
  if (!(func_name %in% options)) stop_h(
    sprintf("Function '%s' is not supported.", func_name),
    find_closest_match(func_name, options, 1)
  )

  y <- type_split_(parts[3])
  if (length(y$val) != 1) stop_h(
    sprintf("Function '%s' expects exactly one argument.", func_name)
  )
  list(func = func_name, key = "x", val = y$val, at = y$at)
}

#' @title
#' Extract augmented function
#' @description
#' Extract augmented function declaration.
#'
#' @param `character()`\cr
#' Vector of characters.
#'
#' @return `list`cr
#' function name, key, variable and type
#'
#' @keywords internal
#' @noRd
extract_func_ <- function(s) {
  # UPDATE: use rlang instead !!!
  # remove spaces
  s <- gsub(" ", "", s)

  if (count_unique_function_calls_(s) > 1) stop_h(
    sprintf("Multiple functions invoked in definition \"%s\"", s)
  )
  else if (grepl("[()]", s)) return(extract_func_call_(s))

  else if (grepl("=", s)) return(extract_func_equality_(s))

  else if (grepl("\\+", s)) return(extract_func_arithmetic_(s, "\\+", "add"))

  else if (grepl("\\*", s)) return(extract_func_arithmetic_(s, "\\*", "mult"))

  else if (grepl("-", s)) return(extract_func_arithmetic_(s, "-", "sub"))

  else if (grepl("/", s)) return(extract_func_arithmetic_(s, "/", "div"))

  else if (grepl("\\^", s)) return(extract_func_arithmetic_(s, "\\^", "pow"))

  else if (!is.na(s)) {
    # use the convention of equality
    y <- type_split_(s)
    if (length(y$val) == 1) { 
      return(list(func = "eq", key = "x", val = y$val, at = y$at))
    }
  }

  # else it's not a call
  list(func = NA, key = NA, val = NA, at = NA)
}
