#      __                                 
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ / 
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /  
#  /_/ /_/\____/\____/____/\___/ ____/_/   
#                             /_/       

# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# utils.R
#

stop_h <- function(...) {
  base::stop(paste(unlist(list(...)), collapse = "\n"),
             call. = FALSE)
}
warning_h <- function(...) {
  base::warning(paste(unlist(list(...)), collapse = "\n"),
                call. = FALSE)
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

#' @description
#' Collapse a character vector into string of quoted words.
#'
#' @param `character()`\cr
#' Vector of characters.
#'
#' @return
#' String formatted as a list of quoted words.
#'
#' @noRd
character_to_string_list <- function(mylist) {
  return(paste(lapply(mylist, function(x) paste0("\"", x, "\"")),
               collapse = ", "))
}

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
#' @noRd
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
        sprintf("Available option%s %s.", 
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

# Format NA types with defaults according to values
default_types_ <- function(typs, vals) {
  for (i in which(is.na(typs))) {
    y <- vals[[i]]
    if (is.na(y)) {
      typs[i] <- NA
    } else if (is.numeric(y)) {
      typs[i] <- "S"
    } else if (is.character(y)) {
      if (is.na(suppressWarnings(as.numeric(y)))) {
        typs[i] <- "P"
      } else {
        typs[i] <- "S"
      }
    } else {
      typs[i] <- "P"
    }
  }
  return(typs)
}

# Split a typed string or character vector
str_split_ <- function(lst) {
  values <- sapply(
    strsplit(as.character(lst), "\\@"), function(x) x[1]
  )
  types <- sapply(
    strsplit(as.character(lst), "\\@"), function(x) x[2]
  )
  return(list(var = values, typ = default_types_(types, values)))
}

# Join elements of a dictionary to pass to csv
str_join_ <- function(x) {
  return(paste(x, collapse = "/"))
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
#' @noRd
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
#' `coda::mcmc` | `coda::mcmc.list`\cr
#'
#' @noRd
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
#' `data.frame`)
#'
#' @noRd
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
