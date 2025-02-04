#
# utils.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

internal_setup_error <- function() {
  stop(cli::format_error(c(
    "Failed to setup files necessary to model.",
    i = "Check your user permissions."
  )))
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
        typs[i] <- "D"
      } else {
        typs[i] <- "S"
      }
    } else {
      typs[i] <- "D"
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

format_string <- function(name) {
  return(gsub("[[:punct:]]", "", name))
}

# Generate a random string
rnd_string <- function(size = 8) {
  return(paste0(sample(c(letters, LETTERS, 0:9), size, replace = TRUE),
                collapse = ""))
}

# Returns the path/to/models folder
template_models_dir <- function() {
  return(file.path(system.file(package = "housepi"), "models"))
}

# Source code files
source_files <- function(path) {
  return(file.path(path, paste0("model.", c("cpp", "h"))))
}

# Verify source files exist
assert_source_files <- function(path) {
  for (file_name in source_files(path)) {
    if (!file.exists(file_name)) {
      stop(cli::format_error(c(
        sprintf("File \"%s\" was not found.", file_name),
        i = "Check your model directory"
      )))
    }
  }
  invisible()
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
