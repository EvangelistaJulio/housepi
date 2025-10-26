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

# check if suggested library is installed
require_suggested_library <- function(libname) {
  if (!requireNamespace(libname, quietly = TRUE)) stop_h(
    sprintf("Attempt to call missing library '%s'.", libname),
    sprintf("Install it running: install.packages('%s')", libname)
  )
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
  if (length(mylist) > max_print) return(paste0(sms, ", ..."))
  else return(paste0(sms, "."))
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
  paste0(sms, ".")
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
#' @param ...
#' Named arguments with respective list of definition details.
#' @param from
#' A named `list` of definitions (alternative to named arguments).
#'
#' @return
#' `list` of `list`
#'
#' @keywords internal
#' @noRd
helper_load_definition_list <- function(..., from = NULL) {
  definitions <- list(...)

  if (length(definitions) == 0) {
    if (length(from) == 0) stop_h("No definition provided while loading list.")
    else if (is.list(from)) definitions <- from
    else stop_h("Argument 'from' should be a list.")
  } else if (any(names(definitions) == "")) {
    stop_h("Found definition(s) without name.")
  }

  return(definitions)
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
      info[[col]] <- list(class = cls,
                          levels = levels(y))
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

############################################################################
############################################################################
############################# Definition Encoder ###########################
############################################################################
############################################################################

#' @title Definition Encoder
#' @description
#' Helper with definitions.
#'
#' @keywords internal
#' @noRd
DefEncoder <- R6::R6Class(
  classname = "DefEncoder",

  private = list(
    # If it's a call to "(" with a single argument that's a symbol, unwrap it.
    .strip_useless_parens = function(x) {
      while (rlang::is_call(x) &&
               rlang::call_name(x) == "(" &&
               length(rlang::call_args(x)) == 1) {
        inner_arg <- rlang::call_args(x)[[1]]
        if (rlang::is_symbol(inner_arg)) x <- inner_arg
        else if (rlang::call_name(inner_arg) == "@") x <- inner_arg
        else break
      }
      x
    },

    # Check type definitions and set defaults.
    .default_types = function(types, values) {
      valid_types <- list("A" = "augmented",
                          "D" = "data",
                          "G" = "group rate",
                          "N" = "numeric",
                          "P" = "parameter")

      for (i in seq_along(types)) {
        k <- types[i]
        y <- values[[i]]

        if (is.na(k)) {
          # no explicit type provided, then format for default types
          if (is.na(y)) types[i] <- NA_character_
          else if (is.numeric(y)) types[i] <- "N"
          else if (is.character(y)) types[i] <- ifelse(
            is.na(suppressWarnings(as.numeric(y))), "A", "N"
          )
          else types[i] <- "A"

        } else {
          # explicit type checking
          k <- toupper(k)
          if (!(k %in% names(valid_types))) stop_h(
            paste0("Invalid type '", k, "' for value '", y, "'."),
            paste(
              sprintf("'%s'", names(valid_types)),
              sprintf("(%s)", unlist(valid_types)),
              collapse = ", "
            )
          )
        }
      }

      types
    },

    # Add default types if not a call.
    .format_types = function(func_args) {
      output <- NA_character_
      if (!(is.null(func_args) || all(is.na(func_args)))) {
        splitted <- strsplit(as.character(func_args), "\\@")
        output <- sapply(splitted, function(x) {
          ifelse(
            stringr::str_count(x[1], "\\|"), x[1],
            sprintf("%s|%s", x[1], private$.default_types(x[2], x[1]))
          )
        }, USE.NAMES = FALSE)
      }
      output
    },

    # Format a call to match housepi standard typed definitions
    .format_identity = function(func_args) {
      if (is.null(func_args)) NA_character_
      else private$.format_types(func_args)
    },

    # Format a call to match housepi standard typed definitions
    .format_call = function(func_name, func_args) {
      paste0(c(
        func_name,
        paste0(names(func_args), collapse = ":"),
        paste0(
          sprintf("{%s}", private$.format_types(func_args)),
          collapse = ":"
        )
      ), collapse = "|")
    },

    # Helper to extract functin from call.
    .helper_func_from_call = function(call) {
      # strip meaningless parentheses from arguments
      args <- lapply(rlang::call_args(call), private$.strip_useless_parens)

      # list operands as literals
      func_args <- character()
      if (length(args) > 0) {
        for (i in seq_along(args)) {
          y <- args[[i]]
          func_args[[paste0("x", i)]] <- ifelse(
            rlang::is_call(y) && rlang::call_name(y) != "@",
            private$.helper_func_from_call(y),
            rlang::as_label(y)
          )
        }
        names(func_args)[1] <- "x"
      }

      # map to function name
      func_name <- switch(
        rlang::call_name(call),
        `+` = if (length(args) == 1) NA_character_ else "sum",
        `-` = if (length(args) == 1) "minus" else "sub",
        `*` = "prod",
        `/` = "div",
        `^` = "pow",
        rlang::call_name(call) # default
      )

      if (is.na(func_name)) private$.format_identity(func_args)
      else private$.format_call(func_name, func_args)
    },

    # Extract function or identiyt call from rlang expression.
    .extract = function(expr) {
      # read input expression
      if (is.null(expr)) {
        return(private$.format_identity(NULL))
      } else if (inherits(expr, "formula")) {
        call <- rlang::f_rhs(expr)
      } else if (rlang::is_call(expr)) {
        call <- expr
      } else {
        return(private$.format_identity(c(x = expr)))
      }

      if (rlang::is_symbol(call)) {
        # case of a symbol actually looking as a call
        return(private$.format_identity(c(x = call)))
      } else if (rlang::call_name(call) == "@") {
        # case of a raw typed identity
        return(private$.format_identity(c(x = call)))
      } else {
        # standard call
        return(private$.helper_func_from_call(call))
      }
      invisible()
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field x `character(1)`\cr
    #' Encoding after initialization.
    x = NA_character_,

    #' @description Flag to tell if object is a call or a value.
    is_call = function() {
      # NOTE: default function are 0 '|'
      isTRUE(stringr::str_count(self$x, "\\|") > 1)
    },

    #' @description Return the value, if applicable.
    #' @param func_name
    #' A `character(1)` with the function name.
    #'
    #' @return `character(1)` or `NA`
    combine_to_call = function(func_name) {
      if (DefEncoder$new(func_name)$is_call()) stop_h(
        "Function name is actually a call, impossible to combine."
      )
      sprintf("%s|%s|%s", func_name,
              paste(names(self$x), collapse = ":"),
              paste(sprintf("{%s}", self$x), collapse = ":"))
    },

    #' @description
    #' New object of class
    #' [HousepiAugmented][housepi::HousepiAugmented].
    initialize = function(expr) {
      if (is.list(expr)) self$x <- sapply(expr, private$.extract)
      else self$x <- private$.extract(expr)
      invisible()
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title
#' Walk Encoded Function
#' @description
#' Recursively walk an encoded function to extract all values of a given type.
#'
#' @param input
#' A `character(1)` with encoded function string.
#' @param type
#' A `character(1)` with the type to extract.
#'
#' @return `character()`
#' @keywords internal
#' @noRd
walk_encoded_func <- function(input, type) {
  if (is.na(input) || is.null(input) || input == "") return(NULL)
  blocks <- extract_blocks(input)
  if (length(blocks) > 0) {
    output <- c()
    for (block in blocks) output <- c(output, evaluate_block(block, type))
    if (length(output) == 0) NULL
    else unique(output)
  } else {
    # terminal node
    parts <- unlist(strsplit(input, "\\|"))
    if (length(parts) == 2 && parts[2] == type) parts[1]
    else NULL
  }
}

# Helper to extract balanced {...} blocks at the top level.
extract_blocks <- function(s) {
  blocks <- list()
  depth <- 0
  start <- NULL
  for (i in seq_along(strsplit(s, "")[[1]])) {
    ch <- substr(s, i, i)
    if (ch == "{") {
      if (depth == 0) start <- i
      depth <- depth + 1
    } else if (ch == "}") {
      depth <- depth - 1
      if (depth == 0 && !is.null(start)) {
        blocks <- c(blocks, substr(s, start, i))
        start <- NULL
      }
    }
  }
  blocks
}

# Recursive helper function to evaluate blocks.
evaluate_block <- function(block, type) {
  # Remove curly braces
  block <- sub("^\\{", "", block)
  block <- sub("\\}$", "", block)

  # split the block into individual conditions
  conditions <- unlist(strsplit(block, ":"))
  output <- c()

  # Evaluate each condition
  for (cond in conditions) {
    if (grepl("\\{", cond)) {
      # recurse on nested blocks
      output <- c(output, walk_encoded_func(cond, type))
    } else {
      parts <- unlist(strsplit(cond, "\\|"))
      if (length(parts) == 2 && parts[2] == type) output <- c(output, parts[1])
    }
  }
  output
}

#' @title
#' Make Definition Readable
#' @description
#' Convert encoded definition strings into a more human-readable format.
#'
#' @param defs
#' A `data.frame()` with encoded definition strings.
#'
#' @return `data.frame()`
#' @keywords internal
#' @noRd
make_def_readable <- function(defs) {
  .extract_numeric <- function(x) {
    sapply(x, function(element) {
      if (is.na(element)) {
        NA_character_
      } else if (is.character(element)) {
        z <- strsplit(element, "\\|")[[1]]
        if (length(z) == 2 && z[2] == "N") as.numeric(z[1])
        else element
      } else {
        element
      }
    }, USE.NAMES = FALSE)
  }

  .extract_args <- function(x) {
    y <- sapply(x, function(element) {
      if (is.na(element)) {
        c(NA_character_, NA_character_)
      } else if (is.character(element)) {
        z <- strsplit(element, "\\|")[[1]]
        if (length(z) > 2) c(z[2], paste(z[3:length(z)], collapse = "|"))
        else c(NA_character_, NA_character_)
      } else {
        c(element, NA_character_)
      }
    }, USE.NAMES = FALSE)

    as.data.frame(t(y))
  }

  .extract_func <- function(x) {
    sapply(x, function(element) {
      if (is.na(element)) {
        NA_character_
      } else if (is.character(element)) {
        z <- strsplit(element, "\\|")[[1]]
        if (length(z) > 2) z[1]
        else element
      } else {
        element
      }
    }, USE.NAMES = FALSE)
  }

  output <- as.data.frame(lapply(defs, .extract_numeric))
  for (col in colnames(defs)) {
    if (isTRUE(any(stringr::str_count(output[[col]], "\\|") > 1))) {
      pos <- match(col, names(output))

      output <- cbind(
        output[1:(pos - 1)],
        stats::setNames(list(.extract_func(output[[col]])), col),
        stats::setNames(
          .extract_args(output[[col]]),
          paste0(
            ifelse(col %in% c("prior", "model"),
                   "arg", paste0(col, "_arg")),
            c("_names", "_values")
          )
        ),
        output[(pos + 1):ncol(output)]
      )
    }
  }

  output
}
