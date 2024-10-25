#
# utils.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

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
  vals <- sapply(
    strsplit(as.character(lst), "\\@"), function(x) x[1]
  )
  types <- sapply(
    strsplit(as.character(lst), "\\@"), function(x) x[2]
  )
  return(list(values = vals, types = default_types_(types, vals)))
}

# Join elements of a dictionary to pass to csv
str_join_ <- function(x) {
  return(paste(x, collapse = "/"))
}