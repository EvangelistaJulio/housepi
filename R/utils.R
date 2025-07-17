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
#' @keywords internal
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
#' @keywords internal
warning_h <- function(...) {
  base::warning(paste(unlist(list(...)), collapse = "\n"),
                call. = FALSE)
  invisible()
}

internal_setup_error <- function() {
  stop_h(
    "Failed to setup necessary files.",
    "Check your user permissions."
  )
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
  if (len_list == 0) return(NULL)
  else if (len_list == 1) sms <- mylist
  else {
    sms <- paste(mylist[1:(len_list-1)], collapse = ", ")
    sms <- paste(sms, "and", mylist[len_list])
  }
  return(paste0(sms, "."))
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
    if (!file.copy(from, to, overwrite = TRUE, ...)) {
      # user permission error (tipically)
      internal_setup_error()
    }
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
#' @param from (named `list()`\cr
#' Named list of definitions (alternative to named arguments).
#'
#' @return
#' `list() of list()`
#'
#' @keywords internal
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
############################################################################
############################################################################

# remove punctuation from string
format_string <- function(name) {
  return(gsub("[[:punct:]]", "", name))
}

# Generate a random string
rnd_string <- function(size = 8) {
  return(paste0(sample(c(letters, LETTERS, 0:9), size, replace = TRUE),
                collapse = ""))
}

# Source code files
source_files <- function(path) {
  return(file.path(path, paste0("model.", c("cpp", "h"))))
}

# return time of modification truncated
trunc_mtime <- function(path) {
  return(max(trunc(file.info(path)$mtime)))
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
untracked_changes <- function(source_dir, dll) {
  last_mod <- trunc_mtime(source_files(source_dir))

  # last changed should be older than installation time
  housepi_dir <- file.path(system.file(package = "housepi"), "template")
  check_src_files <- list.files(housepi_dir, full.names = TRUE)
  install_time <- trunc_mtime(check_src_files)

  return(any(trunc_mtime(dll) < c(last_mod, install_time)))
}

# Optional helper files
helper_files <- function(path) {
  helper_names <- paste0("template_",
                         c("param.R", "info.Rd", "augm.R"))
  return(file.path(path, helper_names))
}

# Associate dynamic library extension to name
dll_file_ext <- function(name = "") {
  return(paste0(name, .Platform$dynlib.ext))
}

# Templated name convention
template_name <- function(path) {
  return(paste0("housepi", format_string(basename(path))))
}

# Add random string to templated name
random_name <- function(path) {
  return(paste0(template_name(path), rnd_string()))
}

# DLL fyle with templated name convention
template_dll_file <- function(path) {
  return(file.path(path, dll_file_ext(template_name(path))))
}

# Get template name from DLL file
get_template_from_dll <- function(dll) {
  if (is.null(dll) || is.na(dll) || !file.exists(dll)) return(NA)
  else return(tools::file_path_sans_ext(basename(dll)))
}

# Get Rcpp module from DLL file
get_rcpp_module <- function(module_name, model = NULL) {
  dll <- ifelse(is.null(model), NA, model$path_to_dll)
  if (is.na(dll)) {
    housepi_dir <- file.path(system.file(package = "housepi"), "template")
    dll <- file.path(housepi_dir, "src", dll_file_ext("housepienvironment"))
  }
  dyn.load(dll)
  package <- get_template_from_dll(dll)
  return(Rcpp::Module(module_name, PACKAGE = package))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# split into a list of matrix
split_traces_by_chain_ <- function(x) {
  if (length(x) == 0) stop_h("Traces are empty.")
  return(lapply(
    split(seq_len(nrow(x)), x[, "_chain_"]),
    function(idx) x[idx, colnames(x) != "_chain_", drop = FALSE]
  ))
}

#' @title
#' Convert `housepi` Rcpp traces to coda.
#' @description
#' Convert traces export from a Rcpp `housepi` model to coda.
#'
#' @param trs
#' Traces as exported from Rcpp `housepi` model.
#'
#' @return
#' `posterior::draws`
#'
#' @keywords internal
cpp_traces_to_posterior <- function(trs) {
  if (methods::is(trs, "draws")) {
    return(trs)

  } else if (is.matrix(trs) || (is.data.frame(trs))) {
    list_trs <- split_traces_by_chain_(trs)

    # create empty array
    n_chain <- length(list_trs)
    par_names <- colnames(list_trs[[1]])
    n_iter <- nrow(list_trs[[1]])
    iter_names <- rownames(list_trs[[1]])
    array_trs <- array(NA_real_, dim = c(n_iter, n_chain, length(par_names)),
                       dimnames = list(iter_names, seq_len(n_chain), par_names))
    # associate traces from Cpp
    for (i in seq_along(list_trs)) array_trs[, i, ] <- list_trs[[i]]
    output <- posterior::as_draws_df(array_trs)
    output$.iteration <- as.integer(rownames(trs))
    return (output)

  } else stop_h(
    "Unrecognized type for 'housepi' Rcpp traces export."
  )
}

#' @title
#' Convert `housepi` Rcpp traces to coda.
#' @description
#' Convert traces export from a Rcpp `housepi` model to coda.
#'
#' @param trs
#' Traces as exported from Rcpp `housepi` model.
#'
#' @return
#' `coda::mcmc` | `coda::mcmc.list`
#'
#' @keywords internal
cpp_traces_to_coda <- function(trs) {
  require_suggested_library("coda")
  if ((methods::is(trs, "mcmc") || methods::is(trs, "mcmc.list"))) {
    return(trs)

  } else if (is.matrix(trs) || (is.data.frame(trs))) {
    list_trs <- split_traces_by_chain_(trs)
    sweeps <- as.integer(rownames(list_trs[[1]]))
    list_coda_mcmc <- lapply(list_trs, function(x) coda::mcmc(
      data = x,
      start = sweeps[1],
      end = sweeps[length(sweeps)],
      thin = sweeps[2] - sweeps[1]
    ))
    if (length(list_coda_mcmc) == 1) return(list_coda_mcmc[[1]])
    return(coda::as.mcmc.list(list_coda_mcmc))

  } else stop_h(
    "Unrecognized type for 'housepi' Rcpp traces export."
  )
}

############################################################################
############################################################################
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
#' @param values `list()`\cr
#' List of values.
#'
#' @return `character()`\cr
#' Vector of characters.
#'
#' @keywords internal
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

  return(types)
}
# @title
#' Type split
#' @description
#' Split a typed string or character vector
#'
#' @param lst `character()`\cr
#' Vector of characters.
#'
#' @return `list()`\cr
#' List with two elements: value and type.
#'
#' @keywords internal
type_split_ <- function(lst) {
  if (all(is.na(lst))) {
    return(list(val = NA, at = NA))

  } else {
    splitted_lst <- strsplit(as.character(lst), "\\@")
    values <- sapply(splitted_lst, function(x) x[1])
    types <- sapply(splitted_lst, function(x) x[2])
    return(list(
      val = values,
      at = default_types_(types, values)
    ))
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
type_join_ <- function(x) {
  return(paste(x, collapse = "/"))
}

#' @title
#' Format setup input
#' @description
#' Format elements of a list to pass to import_setup Rcpp method.
#'
#' @param defs `list()`\cr
#' List of setup command
#'
#' @return `character()`
#'
#' @keywords internal
fmt_setup_input_ <- function(defs) {
  # change scientific notation locally
  withr::local_options(list(scipen = 500))
  # correct id trace vector
  if ("id_trace" %in% names(defs)) {
    y <- unlist(defs$id_trace)
    if (length(y) == 0) defs$id_trace <- NA
    else defs$id_trace <- type_join_(unique(y))
  }
  return(data.frame(lapply(defs, as.character)))
}

#' @title
#' Format data input
#' @description
#' Format elements of dataset (or datalist) to pass to import Rcpp method.
#'
#' @param df `data.frame()`\cr
#' Dataest or datalist.
#'
#' @return `character()`
#'
#' @keywords internal
fmt_data_input_ <- function(df) {
  if (is.null(df)) return(NULL)
  else return(data.frame(lapply(df, as.character)))
}

#' @title
#' Cpp Data Refactor
#' @description
#' Restore factored quanties on dataset return from Cpp implementation.
#'
#' @param df `data.frame`\cr
#' Dataset as output from rcpp model object.
#' @param levels `list()`\cr
#' Levels of the factors in the dataset.
#' @param headers `character()`\cr
#' Headers of the columns.
#'
#' @return `data.frame`
#'
#' @keywords internal
cpp_data_refactor = function(df, levels, headers) {
  if (length(df) > 0) {
    # re-factor
    if (length(levels) > 0) {
      for (col in names(levels)) df[[col]] <- factor(
        df[[col]],
        levels = seq_along(levels[[col]]),
        labels = levels[[col]]
      )
    }
    # re-order data and eventually extra columns
    df <- df[c(headers, setdiff(colnames(df), headers))]
  }
  return(df)
}

############################################################################
############################################################################
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
#' @return `list()`
#'
#' @keywords internal
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
#' @param `character()`\cr
#' Vector of characters.
#'
#' @return `list()`
#'
#' @keywords internal
extract_func_equality_ <- function(s) {
  # split according to the regular expression
  y <- strsplit(s, "=")[[1]]
  if (length(y) != 2) stop_h(
    "Invalid number of arguments for equality.",
    "Expected 1 argument, found ", length(y), "."
  )
  
  y <- type_split_(y[-1])

  return(list(
    func = "eq", key = "x",
    val = type_join_(y$val), at = type_join_(y$at)
  ))
}

#' @title
#' Extract augmented function basic
#' @description
#' Extract basic augmented function declaration
#'
#' @param `character(1)`\cr
#' Vector of characters.
#' @param `character(1)`\cr
#' Regular expression to split the string.
#' @param `character(1)`\cr
#' Function name.
#'
#' @return `list()`
#'
#' @keywords internal
extract_func_arithmetic_ <- function(s, reg_expr, func_name) {
  # split according to the regular expression
  y <- strsplit(s, reg_expr)[[1]]
  if (length(y) != 2) stop_h(
    "Invalid number of arguments for ", func_name, ".",
    "Expected 2 argument, found ", length(y), "."
  )
  y <- type_split_(y)

  return(list(
    func = func_name, key = type_join_(c("x", "y")),
    val = type_join_(y$val), at = type_join_(y$at)
  ))
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
#' `list()`
#'
#' @keywords internal
extract_func_call_ <- function(s) {
  options <- c("pow", "sqrt", "exp", "log", "sin", "cos", "tan")

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
  return(list(func = func_name, key = "x", val = y$val, at = y$at))
}

#' @title
#' Extract augmented function
#' @description
#' Extract augmented function declaration.
#'
#' @param `character()`\cr
#' Vector of characters.
#'
#' @return `list()`cr
#' function name, key, variable and type
#'
#' @keywords internal
extract_func_ <- function(s) {
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

  else return(list(func = NA, key = NA, val = NA, at = NA))
}


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @title
#' Test Housepi Function
#' @description
#' Call housepi registered function.
#' This function internally calls a large C++ class, 
#' thus it is not optimized for repeated  operations.
#'
#' @param name `character(1)`\cr
#' Name of the function.
#' @param x `numeric()`\cr
#' Input vector.
#' @param args `numeric()`\cr
#' Additional arguments.
#' @param min `double(1)`\cr
#' Value of the minimum.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param model `HouseholdModel`\cr
#' Custom household model to load from.
#'
#' @return
#' `numeric()`
#'
#' @export
test_function <- function(name, x, args = numeric(), 
                          min = NA, max = NA, model = NULL) {
  # note: for simple function no reshaping -> discr = FALSE
  get_rcpp_module("housepi_function", model)$FunctionHousepi$new(
    name, FALSE, unlist(args), c(min = min, max = max)
  )$call(x)
}

#' @title
#' Test Housepi Distribution
#' @description
#' Call housepi registered distribution.
#' This function internally calls a large C++ class, 
#' thus it is not optimized for repeated  operations.
#'
#' @param name `character(1)`\cr
#' Name of the distribution.
#' @param x `numeric()`\cr
#' Input vector.
#' @param discr `logical(1)`\cr
#' If `TRUE` discretize the distribution.
#' @param args `numeric()`\cr
#' Additional arguments.
#' @param min `double(1)`\cr
#' Value of the minimum.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param model `HouseholdModel`\cr
#' Custom household model to load from.
#'
#' @return
#' `numeric()`
#'
#' @export
test_distribution <- function(name, x, discr = FALSE, args = numeric(), 
                              min = NA, max = NA, model = NULL) {    
  get_rcpp_module("housepi_distribution", model)$DistributionHousepi$new(
    name, discr, unlist(args), c(min = min, max = max)
  )$call(x)
}

#' @title
#' Test Housepi Infectivity Profile
#' @description
#' Call housepi registered infectivity profile.
#' This function internally calls a large C++ class, 
#' thus it is not optimized for repeated  operations.
#'
#' @param name `character(1)`\cr
#' Name of the infectivity profile
#' @param x `numeric()`\cr
#' Input vector.
#' @param discr `logical(1)`\cr
#' If `TRUE` discretize the infectivity profile.
#' @param args `numeric()`\cr
#' Additional arguments.
#' @param max `double(1)`\cr
#' Value of the maximum.
#' @param model `HouseholdModel`\cr
#' Custom household model to load from.
#'
#' @return
#' `numeric()`
#'
#' @export
test_infec_profile <- function(name, x, discr = FALSE, args = numeric(), 
                               max = NA, model = NULL) {          
  # UPDATE: the minimum should be imposed from within         
  get_rcpp_module("housepi_infec_profile", model)$InfecProfileHousepi$new(
    name, discr, unlist(args), c(min = ifelse(discr, 1, 0), max = max)
  )$call(unlist(x))
}
