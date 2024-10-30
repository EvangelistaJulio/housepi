#
# model_template.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title household model template
#' @description
#' Copy a default model to a given directory.
#' @param which (`character(1)`)\cr
#' Name of template model. If `NULL`, list of available templates.
#' @param path (`character(1)`)\cr
#' Where to copy the template scripts. If `NULL`, show template model description in RStudio Help window.
#' @export
household_model_template <- function(which = NULL, path = NULL) {
  # system path to housepi
  models_dir <- file.path(system.file(package = "housepi"), "models")
  
  if (is.null(which)) {
    cli::cli_alert_info("Template model names:")
    ulid <- cli::cli_ul()
    for (name in list.files(models_dir)) {
      cli::cli_li(name)
    }
    cli::cli_end(ulid)
    
  } else if (is.character(which)) {
    
    if (which %in% list.files(models_dir)) {
      if (is.null(path)) {
        #----- show model description -----#
        templ_info_file <- file.path(models_dir, which, "template_info.Rd")
        if (!file.exists((templ_info_file))) {
          stop(cli::format_error("Impossible to show model description. No template info were found."))
        }
        rstudioapi::previewRd(templ_info_file)
      } else {
        
        #----- copy files -----#
        dir.create(path, recursive = TRUE, showWarnings = TRUE)
        if (dir.exists(path)) {
          # FIXME : maybe don't copy the copyright section
          for (f in c("model.cpp", "model.h", "template_params.R", "template_info.Rd")) {
            from_file <- file.path(models_dir, which, f)
            to_file <- file.path(path, f)
            if (!file.copy(from_file, to_file, overwrite = TRUE)) {
              stop(cli::format_error(
                c(
                  "Failed to copy model files.",
                  i = "Check your user permissions."
                )
              ))
            }
          }
          
        } else {
          stop(cli::format_error(
            c(
              paste0("Impossibile to create directory `", path, "`."),
              i = "Check your user permissions."
            )
          ))
        }
      }
      
    } else {
      stop(cli::format_error(
        c(
          "Provided model is not a template.",
          i = "To see the available template models, run: template_model()"
        )
      ))
    }

  } else {    
    stop(cli::format_error("Template name should be a character."))
  }
}


household_model_generate <- function(model,
                                     save = NA,
                                     quiet = TRUE,
                                     force = FALSE,
                                     debug = FALSE) {
  # system path to housepi
  housepi_dir <- system.file(package = "housepi")

  # Returns the path/to/DLL file supposedly generated in this session
  executable <- function(exec_dir) {
    return(file.path(exec_dir, "src", paste0("housepienvironment", .Platform$dynlib.ext)))
  }
  # model environment name
  env_name <- paste0("housepi", gsub("_", "", model))
  # temporary session directory where to store model environment
  exec_dir <- file.path(tempdir(), env_name)
  
  if (dir.exists(model)) {
    #----- check if folder contains executable -----#
    src_fldr <- file.path(model)
    if (file.exists((executable(model)))) {
      exec_dir <- model
      src_fldr <- file.path(src_fldr, "src")
    }

    #----- custom user model -----#
    source_files <- file.path(src_fldr, c("model.cpp", "model.h"))
    opt_src_files <- file.path(src_fldr, c("template_params.R", "template_info.Rd", "template_augm.R"))
    # check that header and implementation files exist
    for (file_name in source_files) {
      if (!file.exists(file_name)) {
        stop(cli::format_error(c(
          paste0("A file '", file_name, "' was not found."),
          i = paste0("Check your model folder.")
        )))
      }
    }
    
  } else if (model %in% list.files(file.path(housepi_dir, "models"))) {
    #----- template model -----#
    source_files <- file.path(housepi_dir,
                              "models",
                              model,
                              c("model.cpp", "model.h", "template_params.R", "template_info.Rd"))
    opt_src_files <- file.path(housepi_dir, "models", model, c("template_augm.R"))
    
  } else {
    #----- error -----#
    stop(cli::format_error(paste0(
      "'", model, "' is not a folder nor a template name."
    )))
  }

  # add optional files if they exist
  for (file_name in opt_src_files) {
    if (file.exists(file_name)) {
      source_files <- c(source_files, file_name)
    }
  }

  # save option has the priority
  if ((!is.na(save)) && (is.character(save))){
    exec_dir <- save
  }
  
  # ----- Create the model package ----- #
  
  # option to force the removal of the executable
  if ((file.exists(executable(exec_dir)) && (force == TRUE))) {
    unlink(executable(exec_dir))
  }
  
  # deal with existing version
  if (file.exists((executable(exec_dir)))) {
    # use cached version # TODO: to be tuned with phgbuild compile_dll
    cli::cli_alert_info(
      "Using cached {.pkg {pkgload::pkg_name(exec_dir)}}"
    )

  } else {

    # copy the package to the new directory
    dir.create(exec_dir,
                recursive = TRUE,
                showWarnings = FALSE)
    if (!dir.exists(exec_dir)) {
      stop(cli::format_error(
        c(
          "Impossibile to create a temporary directory for the model.",
          i = "Check your user permissions."
        )
      ))
    }

    # copy single package files
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(housepi_dir, "template", f)
      to_file <- file.path(exec_dir, f)
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c(
            "Failed to setup internal files.",
            i = "Check your user permissions."
          )
        ))
      }
    }

    # copy folders from template
    for (fldr in c("src", "R")) {
      from_path <- file.path(housepi_dir, "template", fldr)
      if (!file.copy(from_path,
                      exec_dir,
                      recursive = TRUE,
                      overwrite = TRUE)) {
        stop(cli::format_error(
          c(
            "Failed to setup internal files.",
            i = "Check your user permissions."
          )
        ))
      }
    }

    # copy source files in the new folder
    for (from_file in source_files) {
      to_file <- file.path(exec_dir, "src", basename(from_file))
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c(
            "Failed to setup implementation files.",
            i = "Check your user permissions."
          )
        ))
      }
    }

    # ----- compile package ----- #
    # TODO : how to enhance pre-compilation?
    pkgbuild::compile_dll(exec_dir, quiet = quiet, debug = debug)
  }

  return(list(dll = executable(exec_dir), path = exec_dir))
}