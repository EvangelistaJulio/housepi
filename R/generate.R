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

# user permission error (tipically)
internal_error <- function() {
  stop(cli::format_error(c(
    x = "Failed to setup internal files.",
    i = "Check your user permissions."
  )))
}

# ensure the model path is a directory or a template
check_model_path <- function(path) {
  if (dir.exists(path)) {
    #===== Folder Contains Source Code =====#
    for (file_name in source_files(path)) {
      if (!file.exists(file_name)) {
        stop(cli::format_error(c(
          x = sprintf("File '%s' was not found.", file_name),
          i = "Check your model directory"
        )))
      }
    }

  } else if (path %in% list.files(template_dir())) {
    #===== Template Model Name =====#
    path <- file.path(template_dir(), path)

  } else {
    #===== ERROR =====#
    stop(cli::format_error(c(
      x = sprintf("'%s' is not a folder nor a template name.", path),
      i = "To check available template names execute: `household_model_template()`"
    )))
  }

  return(path)
}

# @description
# Create temporary folder and copy source code there.
# @param pkg_dir
# Temporary folder where to compile package.
copy_src_code <- function(pkg_dir) { 
  # create temporary directory if necessary
  if(!dir.exists(pkg_dir)) {
    if (!dir.create(pkg_dir, recursive = TRUE, showWarnings = FALSE)) {
      stop(cli::format_error(c(
        x = "Impossibile to create temporary directory for compilation.",
        i = "Check your user permissions."
      )))
    }
  }

  # copy single package files 
  for (f in c("NAMESPACE", "DESCRIPTION")) {
    from_file <- file.path(system.file(package = "housepi"), "template", f)
    to_file <- file.path(pkg_dir, f)
    if (!overwrite(from_file, to_file)) {
      internal_error()
    }
  }

  # copy folders from template
  for (fldr in c("src", "R")) {
    from_path <- file.path(system.file(package = "housepi"), "template", fldr)
    if (!overwrite(from_path, pkg_dir, recursive = TRUE)) {
      internal_error()
    }
  }

  invisible()
}

# save DLL and helper files
save_environment <- function(path, pkg_dir, custom_dir = NA) {

  # create custom directory if necessary
  if (((!is.na(custom_dir))) && is.character(custom_dir)) {
    if(!dir.exists(custom_dir)) {
      if (!dir.create(custom_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop(cli::format_error(c(
          x = sprintf("Impossibile to create directory '%s'.", custom_dir),
          i = "Check that custom path name is valid."
        )))
      }
    }
    out_dir <- custom_dir

  } else {
    out_dir <- pkg_dir
  }

  #----- Copy DLL in the output directory -----#
  from_file <- executable(file.path(pkg_dir, "src"))
  to_file <- executable(out_dir)
  if (!overwrite(from_file, to_file)) {
    internal_error()
  }

  #----- Copy source code withot commented lines as a reference -----#
  lines <- readLines(file.path(pkg_dir, "src", "model.cpp"), warn = FALSE)
  
  # Remove single-line comments (//) 
  lines <- stringr::str_remove_all(lines, "//.*")
  # Remove inline comments (/* ... */)
  lines <- stringr::str_remove_all(lines, "/\\*.*?\\*/")
  # Filter out empty or whitespace-only lines
  cleaned_lines <- lines[stringr::str_detect(lines, "\\S")]
  
  writeLines(cleaned_lines, file.path(out_dir, ".source.txt"))

  #----- Copy helper files in the output directory -----#
  for (from_file in helper_files(path)) {
    if (file.exists(from_file)) {
      to_file <- file.path(out_dir, basename(from_file))

      if (!overwrite(from_file, to_file)) {
        internal_error()
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
generate_model <- function(model_path = ".", store_path = NA,
                           quiet = TRUE, force = FALSE) {

  if ((dir.exists(model_path)) && (file.exists(executable(model_path))) && isFALSE(force)) {
    #===== Folder Contains Executable =====#
    out_dir <- file.path(model_path)
    cli::cli_alert_info("Loading {.pkg {basename(model_path)}}")

  } else {
    #===== Check Model Path =====#
    model_path <- check_model_path(model_path)

    #===== Temporary folder where to store source code =====#
    pkg_dir <- file.path(tempdir(), paste0("housepi", gsub("_", "", model_path)))
    
    if ((file.exists(executable(pkg_dir))) && isFALSE(force)) {
      out_dir <- pkg_dir
      cli::cli_alert_info("Using chached {.pkg {basename(model_path)}}")
      
    } else {
      #===== Create a Temporary Model Source Code =====#
      copy_src_code(pkg_dir)
  
      #===== Copy source files in the new folder =====#
      for (from_file in source_files(model_path)) {
        to_file <- file.path(pkg_dir, "src", basename(from_file))
        if (!overwrite(from_file, to_file)) {
          internal_error()
        }
      }
  
      #===== COMPILE PACKAGE =====#
      # note: tricking pkgbuild::compile_dll and R CMD INSTALL
      #       to avoid recompiling immutable files

      #===== Remove object file =====#
      unlink(file.path(pkg_dir, "src", "model.o"))
  
      #===== Touch installation files =====#
      touch_time <- Sys.time()
      for (f in list.files(file.path(pkg_dir, "src"), full.names = TRUE)) {
        Sys.setFileTime(f, touch_time)
      }
      Sys.setFileTime(file.path(pkg_dir, "src", "model.cpp"), touch_time + 1)
  
      #===== Compile the dynamic library =====#
      cli::cli_alert_info("Compiling {.pkg {basename(model_path)}}")
      pkgbuild::compile_dll(pkg_dir, quiet = quiet, debug = FALSE)
  
      #===== SAVE OUTPUT =====#
      out_dir <- save_environment(model_path, pkg_dir, store_path)
    }
  }

  return(executable(out_dir))
}

# @description
# Re-compiles all required objects of the housepi library.
# Passed to .onLoad in file zzz.R
# @param pkgname
# Name of the package.
# @param force
# Force compilation even if it already exists.
re_compile_library <- function(pkgname, force = FALSE) {
  # compile the core
  template_dir <- file.path(system.file(package = pkgname), "template")
  pkgbuild::compile_dll(template_dir, debug = FALSE)
  
  # compile each model
  models_dir <- file.path(system.file(package = pkgname), "models")
  for (m in list.files(models_dir, full.names = TRUE)) {
    if ((!file.exists(executable(m))) || isTRUE(force)) {
      generate_model(m, store_path = m, quiet = FALSE, force = force)
    }
  }
}