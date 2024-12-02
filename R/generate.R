#
# generate.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

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
        x = sprintf("File '%s' was not found.", file_name),
        i = "Check your model directory"
      )))
    }
  }
  invisible()
}

# Optional helper files
helper_files <- function(path) {
  return(file.path(path, paste0("template_", c("params.R", "info.Rd", "augm.R"))))
}

format_name <- function(name) {
  fmt_name <- gsub("_", "", name)
  fmt_name <- gsub("-", "", fmt_name)
  return(fmt_name)
}

# Templated name convention
template_name <- function(path) {
  return(paste0("housepi", format_name(basename(path))))
}

# Add random string to templated name
random_name <- function(path) {
  return(paste0(template_name(path), rnd_string()))
}

# Associate dynamic library extension to name
dll_file_ext <- function(name) {
  return(paste0(name, .Platform$dynlib.ext))
}

# DLL fyle with templated name convention
template_dll_file <- function(path) {
  return(file.path(path, dll_file_ext(template_name(path))))
}

# user permission error (tipically)
internal_error <- function() {
  stop(cli::format_error(c(
    x = "Failed to setup internal files.",
    i = "Check your user permissions."
  )))
  invisible()
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# @description
# Name of model dynamic library. Compile if needed.
# @param model_path (named `character(1)`)\cr
# @param force (`logical(1)`)\cr
# @param quiet (`logical(1)`)\cr
# @param as_template (`logical(1)`)\cr
generate_dll_model <- function(model_path = ".", force = FALSE, quiet = TRUE, as_template = FALSE) {

  if ((file.exists(template_dll_file(model_path))) && isFALSE(force)) {
    #===== Folder Contains Executable =====#
    cli::cli_alert_info("Loading {.pkg {basename(model_path)}}")
    output <- template_dll_file(model_path)
  }

  else if (dir.exists(model_path)) {
    #===== Folder Contains Source Code =====#
    output <- compile_dll_model(model_path, force, quiet, as_template)

  } else if (model_path %in% list.files(template_models_dir())) {
    #===== Template Model Name =====#
    model_path <- file.path(template_models_dir(), model_path)

    if (file.exists(template_dll_file(model_path))) {
      #----- Template executable already generated -----#
      cli::cli_alert_info("Loading {.pkg {basename(model_path)}}")
      output <- template_dll_file(model_path)

    } else {
      #----- Installation call -----#
      output <- compile_dll_model(model_path, TRUE, TRUE, TRUE)
    }

  } else {
    #===== ERROR =====#
    stop(cli::format_error(c(
      x = sprintf("'%s' is not a folder nor a template name.", model_path),
      i = "To check available template names execute: `household_model_template()`"
    )))
  }
  
  return(output)
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

compile_dll_model <- function(model_path = ".", force = FALSE, quiet = TRUE, as_template = FALSE) {

    # path must contain source files                          
    assert_source_files(model_path)                
    # temporary folder where to store source code
    cache_path <- file.path(tempdir(), paste0("housepi", format_name(basename(model_path))))
    # check for previous cached dll files
    cache_dll_files <- list.files(pattern = .Platform$dynlib.ext)

    if ((length(cache_dll_files) == 1) && isFALSE(force)) {
      cli::cli_alert_info("Using chached {.pkg {basename(model_path)}}")
      output <- cache_dll_files[1]
      
    } else {
      #===== ORGANIZE SOURCE CODE =====#
  
      # clean previous dll files
      unlink(cache_dll_files)

      # create a temporary folder for model source code
      copy_src_code(cache_path)
  
      # copy model specific source files in the new folder =====#
      for (from_file in source_files(model_path)) {
        to_file <- file.path(cache_path, "src", basename(from_file))
        if (!overwrite(from_file, to_file)) {
          internal_error()
        }
      }

      # FIXME : to avoid breaking pre-compilation
      # I attempt to just change the name of the so file after...
      # let's see if it works...
      # FIXME-END

      #===== COMPILE PACKAGE =====#
      # note: tricking pkgbuild::compile_dll and R CMD INSTALL
      #       to avoid recompiling immutable files

      # remove object file
      unlink(file.path(cache_path, "src", "model.o"))
  
      # touch installation files
      touch_time <- Sys.time()
      for (f in list.files(file.path(cache_path, "src"), full.names = TRUE)) {
        Sys.setFileTime(f, touch_time)
      }
      Sys.setFileTime(file.path(cache_path, "src", "model.cpp"), touch_time + 1)
  
      # compile the dynamic library
      cli::cli_alert_info("Compiling {.pkg {basename(model_path)}}")
      pkgbuild::compile_dll(cache_path, quiet = quiet, debug = FALSE)
  
      #===== SAVE OUTPUT =====#
      output <- save_environment(model_path, cache_path, as_template)
  }

  return(output)
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# @description
# Create temporary folder and copy source code there.
# @param pkg_dir
# Temporary folder where to compile package.
copy_src_code <- function(pkg_dir) { 
  # create temporary directory if necessary
  if (dir.exists(pkg_dir)) {
    unlink(list.files(pkg_dir, full.names = TRUE), recursive = TRUE)

  } else {
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

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# save DLL and helper files
save_environment <- function(model_path, cache_path, as_template = FALSE) {

  #===== DYNAMIC LIBRARY =====#

  # move DLL outside of the source directory
  from_file <- file.path(cache_path, "src", dll_file_ext("housepienvironment"))

  if (isTRUE(as_template)) {
    output <- template_dll_file(model_path)
  } else  {
    output <- file.path(cache_path, dll_file_ext(random_name(model_path)))
  }

  if (!overwrite(from_file, output)) {
    internal_error()
  }

  #===== SOURCE CODE FOR CHECKS =====#

  out_dir <- dirname(output)

  # copy source code without commented lines as a reference
  lines <- readLines(file.path(cache_path, "src", "model.cpp"), warn = FALSE)
  
  # Remove single-line comments (//) 
  lines <- stringr::str_remove_all(lines, "//.*")
  # Remove inline comments (/* ... */)
  lines <- stringr::str_remove_all(lines, "/\\*.*?\\*/")
  # Filter out empty or whitespace-only lines
  cleaned_lines <- lines[stringr::str_detect(lines, "\\S")]
  
  writeLines(cleaned_lines, file.path(out_dir, ".source"))

  # copy helper files in the output directory
  for (from_file in helper_files(model_path)) {
    if (file.exists(from_file)) {
      to_file <- file.path(out_dir, basename(from_file))

      if (!overwrite(from_file, to_file)) {
        internal_error()
      }
    }
  }

  #===== CLEANING =====#
  unlink(file.path(cache_path, c("src", "R", "NAMESPACE", "DESCRIPTION")), recursive = TRUE)
  
  return(output)
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# @description
# Re-compiles all required objects of the housepi library.
# Passed to .onLoad in file zzz.R
# @param pkgname
# Name of the package.
# @param force
# Force compilation even if it already exists.
re_compile_library <- function(pkgname, force = FALSE) {
  # compile the core
  template_models_dir <- file.path(system.file(package = pkgname), "template")
  pkgbuild::compile_dll(template_models_dir, debug = FALSE)
  
  # compile each model
  models_dir <- file.path(system.file(package = pkgname), "models")
  for (m in list.files(models_dir, full.names = TRUE)) {
    if ((!file.exists(template_dll_file(m))) || isTRUE(force)) {
      generate_dll_model(m, force = force, quiet = FALSE, as_template = TRUE)
    }
  }
}