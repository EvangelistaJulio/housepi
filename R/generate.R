#
# generate.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# Returns the path/to/DLL file
executable <- function(path) {
  file.path(path, paste0("housepienvironment", .Platform$dynlib.ext))
}

# Returns the path/to/models folder
template_dir <- function() {
  file.path(system.file(package = "housepi"), "models")
}

# Source code files
source_files <- function(path) {
  file.path(path, paste0("model.", c("cpp", "h")))
}

# Optional helper files
helper_files <- function(path) {
  file.path(path, paste0("template_", c("params.R", "info.Rd", "augm.R")))
}

overwrite <- function(from, to, ...) {
  file.copy(from, to, overwrite = TRUE, ...)
}

# Create temporary folder and copy source code there
temp_src_code <- function(path) {

  #===== Temporary folder where to store sorce code =====#
  pkg_dir <- file.path(tempdir(), paste0("housepi", gsub("_", "", path)))

  #===== Copy the package to the new directory =====#
  dir.create(pkg_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(pkg_dir)) {
    stop(cli::format_error(c(
      "Impossibile to create a temporary directory for the model.",
      i = "Check your user permissions."
    )))
  }

  #===== Copy single package files =====#
  for (f in c("NAMESPACE", "DESCRIPTION")) {
    from_file <- file.path(system.file(package = "housepi"), "template", f)
    to_file <- file.path(pkg_dir, f)
    if (!overwrite(from_file, to_file)) {
      stop(cli::format_error(c(
        "Failed to setup internal files.",
        i = "Check your user permissions."
      )))
    }
  }

  #===== Copy folders from template =====#
  for (fldr in c("src", "R")) {
    from_path <- file.path(system.file(package = "housepi"), "template", fldr)
    if (!overwrite(from_path, pkg_dir, recursive = TRUE)) {
      stop(cli::format_error(c(
        "Failed to setup internal files.",
        i = "Check your user permissions."
      )))
    }
  }

  return(pkg_dir)
}

# ensure the model path is a directory or a template
check_model_path <- function(path) {
  if (dir.exists(path)) {
    #===== Folder Contains Source Code =====#
    for (file_name in source_files(path)) {
      if (!file.exists(file_name)) {
        stop(cli::format_error(c(
          sprintf("A file '%s' was not found.", file_name),
          i = "Check your model directory"
        )))
      }
    }

  } else if (path %in% list.files(template_dir())) {
    #===== Template Model Name =====#
    path <- file.path(template_dir(), path)

  } else {
    #===== ERROR =====#
    stop(cli::format_error(
      sprintf("'%s' is not a folder nor a template name.", path)
    ))
  }

  return(path)
}

# save DLL and helper files
save_environment <- function(path, pkg_dir, custom_dir = NA) {
  if (((!is.na(custom_dir))) && is.character(custom_dir)) {
    dir.create(custom_dir, recursive = TRUE, showWarnings = TRUE)
    out_dir <- custom_dir
  } else {
    out_dir <- pkg_dir
  }

  #----- Copy DLL in the output directory -----#
  if (!overwrite(executable(file.path(pkg_dir, "src")), executable(out_dir))) {
    stop(cli::format_error(c(
      "Failed to setup internal files.",
      i = "Check your user permissions."
    )))
  }

  #----- Copy helper files in the output directory -----#
  for (from_file in helper_files(path)) {
    if (file.exists(from_file)) {
      to_file <- file.path(out_dir, basename(from_file))

      if (!overwrite(from_file, to_file)) {
        stop(cli::format_error(c(
          "Failed to setup implementation files.",
          i = "Check your user permissions."
        )))
      }
    }
  }
  
  return(out_dir)
}

# @description
# Name of model dynamic library. Compile if needed.
# @param model_path (named `character(1)`)\cr
# @param store_path (named `character(1)`)\cr
# @param quiet (`logical(1)`)\cr
# @param force (`logical(1)`)\cr
# @param debug (`logical(1)`)\cr
generate_model <- function(model_path = ".", store_path = NA,
                           quiet = TRUE, force = FALSE, debug = FALSE) {

  if ((file.exists(executable(model_path))) && isFALSE(force)) {
    #===== Folder Contains Executable =====#
    out_dir <- file.path(model_path)

    # UPDATE: sugar to be tuned with pkgbuild compile_dll
    cli::cli_alert_info("Using existing {.pkg {housepienvironment}}")

  } else {
    #===== Check Model Path =====#
    model_path <- check_model_path(model_path)

    #===== Create a Temporary Model Source Code =====#
    pkg_dir <- temp_src_code(model_path)

    #===== Copy source files in the new folder =====#
    for (from_file in source_files(model_path)) {
      to_file <- file.path(pkg_dir, "src", basename(from_file))
      if (!overwrite(from_file, to_file)) {
        stop(cli::format_error(c(
          "Failed to setup implementation files.",
          i = "Check your user permissions."
        )))
      }
    }

    #===== COMPILE PACKAGE =====#
    #===== Remove object file =====#
    unlink(file.path(pkg_dir, "src", "model.o"))

    #===== Touch installation files =====#
    touch_time <- Sys.time()
    for (f in list.files(file.path(pkg_dir, "src"), full.names = TRUE)) {
      Sys.setFileTime(f, touch_time)
    }
    Sys.setFileTime(file.path(pkg_dir, "src", "model.cpp"), touch_time + 1)

    #===== Compile the dynamic library =====#
    # note: for the moment force is not within pkgbuild
    pkgbuild::compile_dll(pkg_dir, quiet = quiet, debug = debug)

    #===== SAVE OUTPUT =====#
    out_dir <- save_environment(model_path, pkg_dir, store_path)
  }

  return(executable(out_dir))
}