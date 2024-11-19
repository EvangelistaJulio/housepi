#
# generate.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# Returns the path/to/DLL file 
executable <- function(path) {
  return(file.path(path, paste0("housepienvironment", .Platform$dynlib.ext)))
}

generate_model <- function(model_path = ".", save = NA, quiet = TRUE,
                           force = FALSE, debug = FALSE) {

  if ((file.exists(executable(model_path))) && isFALSE(force)) {
    #===== Folder Contains Executable =====#
    exec_dir <- file.path(model_path)

    # UPDATE: sugar to be tuned with pkgbuild compile_dll
    cli::cli_alert_info(
      "Using existing {.pkg {pkgload::pkg_name(exec_dir)}}"
    )
  }

  else {
    # system path to housepi
    housepi_dir <- system.file(package = "housepi")

    # temporary session directory where to store model environment
    exec_dir <- file.path(tempdir(), paste0("housepi", gsub("_", "", model_path)))

    if (dir.exists(model_path)) {
      #===== Folder Contains Source Code =====#
      source_files <- file.path(
        model_path,         
        c("model.cpp", "model.h"))

      # check that header and implementation files exist
      for (file_name in source_files) {
        if (!file.exists(file_name)) {
          stop(cli::format_error(c(
            paste0("A file '", file_name, "' was not found."),
            i = "Check your model directory"
          )))
        }
      }

      # other optional files
      opt_src_files <- file.path(
        model_path, 
        c("template_params.R", "template_info.Rd", "template_augm.R"))
    
    } else if (model_path %in% list.files(file.path(housepi_dir, "models"))) {
      #===== Template Model Name =====#
      template_dir <- file.path(housepi_dir, "models", model_path)
      source_files <- file.path(template_dir,
                                c("model.cpp", "model.h", "template_params.R", "template_info.Rd"))
      opt_src_files <- file.path(template_dir, c("template_augm.R"))
      
    } else {
      #===== ERROR =====#
      stop(cli::format_error(paste0("'", model_path, "' is not a folder nor a template name.")))
    }

    #----- Add Optional files if they exist -----#
    for (file_name in opt_src_files) {
      if (file.exists(file_name)) {
        source_files <- c(source_files, file_name)
      }
    }

    #----- Save option has the priority -----#
    if ((!is.na(save)) && (is.character(save))){
      exec_dir <- save
    }
  
    #===== Create the Model Source Code =====#
      
    # copy the package to the new directory
    dir.create(exec_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(exec_dir)) {
      stop(cli::format_error(
        c("Impossibile to create a temporary directory for the model.",
          i = "Check your user permissions.")
      ))
    }

    # copy single package files
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(housepi_dir, "template", f)
      to_file <- file.path(exec_dir, f)
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.", i = "Check your user permissions.")
        ))
      }
    }

    # copy folders from template
    for (fldr in c("src", "R")) {
      from_path <- file.path(housepi_dir, "template", fldr)
      if (!file.copy(from_path, exec_dir, recursive = TRUE, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.", i = "Check your user permissions.")
        ))
      }
    }

    #----- Touch installation files -----#
    touch_time <- Sys.time()
    for (f in list.files(file.path(exec_dir,"src"), full.names = TRUE)) {
      Sys.setFileTime(f, touch_time)
    }

    #----- Copy source files in the new folder ------#
    for (from_file in source_files) {
      to_file <- file.path(exec_dir, "src", basename(from_file))
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup implementation files.", i = "Check your user permissions.")
        ))
      }
    }

    #===== COMPILE PACKAGE =====#
    # TODO : how to enhance pre-compilation?
    pkgbuild::compile_dll(exec_dir, quiet = quiet, debug = debug)
  }

  return(list(dll = executable(exec_dir), path = exec_dir))
}