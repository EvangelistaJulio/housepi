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
  if (requireNamespace(libname, quietly = TRUE)) {
    return(TRUE)
  } else {
    stop_h(
      sprintf("Attempt to call missing library \"%s\".", libname),
      sprintf("Install it running: install.packages(\"%s\")", libname)
    )
    return(FALSE)
  }
}

#' @title
#' Character to string list
#' @description
#' Collapse a character vector into string of quoted words.
#'
#' @param mylist`character()`\cr
#' Vector of characters.
#' @param max_print `integer(1)`\cr
#' Maximum number of elements to print (default is 5).
#'
#' @return
#' String formatted as a list of quoted words.
#'
#' @keywords internal
character_to_string_list <- function(mylist, max_print = 5) {
  sms <- paste(lapply(mylist[seq_len(min(max_print, length(mylist)))],
                      function(x) paste0("\"", x, "\"")),
               collapse = ", ")
  if (length(mylist) > max_print) {
    sms <- paste0(sms, ", ...")
  }
  return(sms)
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
find_closest_match <- function(input,
                               options,
                               min_dist = 2) {
  output <- c()

  if ((length(options) > 0) && (min_dist > 0)) {
    dists <- utils::adist(input, options)
    idxs <- apply(dists, 1, min) <= min_dist

    if (any(idxs)) {
      closest <- options[apply(dists[idxs, , drop = FALSE], 1, which.min)]
      output <- sprintf("Did you mean %s?",
                        character_to_string_list(closest))
    }

    if (!all(idxs)) {
      output <- c(output,
        sprintf("Available option%s %s",
                ifelse(length(options) > 1, "s are:", " is"),
                character_to_string_list(options))
      )
    }

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
  if (key %in% names(lst)) {
    return(lst[[key]])
  } else {
    return(default)
  }
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
    if (length(from) == 0) {
      stop_h("No definitions provided.")
    } else if (is.list(from)) {
      definitions <- from
    } else {
      stop_h("Argument 'from' should be a list.")
    }
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

# return time of creation truncated
trunc_ctime <- function(path) {
  return(max(trunc(file.info(path)$ctime)))
}

# say if source files modified after dll or new installation
untracked_changes <- function(model_choice, dll) {
  last_mod <- trunc_ctime(source_files(model_choice))

  # last changed should be older than installation time
  housepi_dir <- file.path(system.file(package = "housepi"), "template")
  check_src_files <- list.files(housepi_dir, full.names = TRUE)
  install_time <- trunc_ctime(check_src_files)

  return(any(trunc_ctime(dll) < c(last_mod, install_time)))
}

# Optional helper files
helper_files <- function(path) {
  helper_names <- paste0("template_",
                         c("params.R", "info.Rd", "augm.R"))
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

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @title
#' Coerce traces to coda.
#' @description
#' Convert HouseholdModel traces to coda.
#'
#' @param trs
#' Traces as exported from HouseholdModel.
#'
#' @return
#' `coda::mcmc` | `coda::mcmc.list`
#'
#' @keywords internal
coerce_traces_to_coda <- function(trs) {
  if (require_suggested_library("coda")) {
    if ((methods::is(trs, "mcmc") || methods::is(trs, "mcmc.list"))) {
      return(trs)

    } else if (is.matrix(trs) || (is.data.frame(trs))) {
      sweeps <- as.integer(rownames(df))
      return(coda::mcmc(data = trs,
                        start = sweeps[1],
                        end = sweeps[length(sweeps)],
                        thin = sweeps[2] - sweeps[1]))

    } else if (is.list(trs)) {
      sweeps <- as.integer(rownames(trs[[1]]))
      df_list <- lapply(trs, function(x) {
        coda::mcmc(
          data = x,
          start = sweeps[1],
          end = sweeps[length(sweeps)],
          thin = sweeps[2] - sweeps[1]
        )
      })
      return(coda::as.mcmc.list(df_list))
    }
  }
}

#' @title
#' Extract results from mcmc traces.
#' @description
#' Collects relevant statistics from the input mcmc traces.
#'
#' @param trs
#' Traces as exported from HouseholdModel.
#' @param cri_level `double(1)`\cr
#' Level for the credible interval.
#' Value between 0 and 1 (default is 0.95).
#'
#' @return
#' `data.frame`
#'
#' @keywords internal
retrieve_traces_results <- function(trs,
                                    cri_level = 0.95) {
  if ((cri_level <= 0) || (cri_level >= 1)) {
    stop_h(
      "Credible interval level must be a value between 0 and 1."
    )
  }

  # change in case of multiple chains
  if (is.list(trs)) trs <- do.call(rbind, trs)

  cri_lower <- 0.5 - 0.5 * cri_level
  cri_upper <- 0.5 + 0.5 * cri_level
  out <- data.frame(
    name = NA,
    mean = apply(trs, 2, mean, na.rm = TRUE),
    sd = apply(trs, 2, sd, na.rm = TRUE),
    cri_lower = apply(trs, 2, quantile, probs = cri_lower, na.rm = TRUE),
    cri_upper = apply(trs, 2, quantile, probs = cri_upper, na.rm = TRUE),
    check.names = FALSE
  )
  out$name <- rownames(out)

  return(out)
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
      if (is.na(y)) {
        types[i] <- NA
      } else if (is.numeric(y)) {
        types[i] <- "S"
      } else if (is.character(y)) {
        if (is.na(suppressWarnings(as.numeric(y)))) {
          types[i] <- "A"
        } else {
          types[i] <- "S"
        }
      } else {
        types[i] <- "A"
      }

    } else {
      # explicit type checking
      k <- toupper(k)
      if (!(k %in% names(valid_types))) {
        valid_message <- paste(
          sprintf("\"%s\"", names(valid_types)),
          sprintf("(%s)", unlist(valid_types)),
          collapse = ", "
        )

        stop_h(
          paste0("Invalid type \"", k, "\" for value \"", y, "\"."),
          valid_message
        )
      }
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
#' @param `character()`\cr
#' Vector of characters.
#'
#' @return `character(1)`\cr
#' String formatted as a list of quoted words.
#'
#' @keywords internal
type_join_ <- function(x) {
  return(paste(x, collapse = "/"))
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
  if (length(y) != 2) {
    stop_h(
      "Invalid number of arguments for equality.",
      "Expected 1 argument, found ", length(y), "."
    )
  }
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
  if (length(y) != 2) {
    stop_h(
      "Invalid number of arguments for ", func_name, ".",
      "Expected 2 argument, found ", length(y), "."
    )
  }
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

  if (count_unique_function_calls_(s) > 0) {
    stop_h(sprintf("Multiple functions invoked in definition \"%s\"", s))
  }

  # check parentheses are matched
  if (nchar(gsub("[^\\(]", "", s)) != nchar(gsub("[^\\)]", "", s))) {
    stop("Unmatched parentheses in input string")
  }

  # extract function name and argument using regex
  m <- regexec("^([a-zA-Z_][a-zA-Z0-9_]*)\\((.*)\\)$", s)
  parts <- regmatches(s, m)[[1]]

  if (length(parts) != 3) {
    stop_h("Input string must be a function call like func(x)")
  }

  func_name <- parts[2]
  if (!(func_name %in% options)) {
    stop_h(
      sprintf("Function '%s' is not supported.", func_name),
      find_closest_match(func_name, options, 1)
    )
  }

  y <- type_split_(parts[3])
  if (length(y$val) != 1) {
    stop_h(sprintf("Function '%s' expects exactly one argument.", func_name))
  }
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

  if (count_unique_function_calls_(s) > 1) {
    stop_h(sprintf("Multiple functions invoked in definition \"%s\"", s))

  } else if (grepl("[()]", s)) {
    return(extract_func_call_(s))

  } else if (grepl("=", s)) {
    return(extract_func_equality_(s))

  } else if (grepl("\\+", s)) {
    return(extract_func_arithmetic_(s, "\\+", "add"))

  } else if (grepl("\\*", s)) {
    return(extract_func_arithmetic_(s, "\\*", "mult"))

  } else if (grepl("-", s)) {
    return(extract_func_arithmetic_(s, "-", "sub"))

  } else if (grepl("/", s)) {
    return(extract_func_arithmetic_(s, "/", "div"))

  } else if (grepl("\\^", s)) {
    return(extract_func_arithmetic_(s, "\\^", "pow"))

  } else {
    return(list(func = NA, key = NA, val = NA, at = NA))
  }
}
