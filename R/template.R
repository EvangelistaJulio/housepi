#      __                                 
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ / 
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /  
#  /_/ /_/\____/\____/____/\___/ ____/_/   
#                             /_/       
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# template.R
#

#' @description
#' Show list of all available template nmodel names.
#'
#' @noRd
show_model_names <- function() {

  cli::cli_alert_info("Available template models:")
  ulid <- cli::cli_ul()

  for (name in list.files(template_models_dir())) {
    cli::cli_li(name)
  }
  cli::cli_end(ulid)

  invisible()
}

#' @title
#' Template household model
#' @description
#' Copy a default model to a given directory.
#'
#' @param name `character(1)`\cr
#' Name of template model. If `NA`, list of available templates.
#' @param path `character(1)`\cr
#' Where to copy the template scripts.
#' If `NA`, show template model description in RStudio Help window.
#' @param replace `logical(1)`\cr
#' If `TRUE`, replaces previous content of the path.
#' Default is `FALSE`.
#'
#' @export
template_household_model <- function(name = NA,
                                     path = NA,
                                     replace = FALSE) {

  if (is.na(name)) {
    #===== Asking for model names =====#
    show_model_names()

  } else if (is.character(name)) {

    if (name %in% list.files(template_models_dir())) {
      if (is.na(path)) {
        #===== Asking for model description =====#
        info_file <- file.path(template_models_dir(),
                               name,
                               "template_info.Rd")
        if (!file.exists((info_file))) {
          stop_h(
            "Impossible to show model description",
            "No template info were foun for the chosen model."
          )
        }
        rstudioapi::previewRd(info_file)

      } else if (any(file.exists(source_files(path))) & isFALSE(replace)) {
        #===== Asking for copying template to folder with model =====#
        warning_h(
          "A model already exists in the specified folder.",
          "Remove its content before creating a new model."
        )

      } else {
        #===== Asking for copying template to empty folder =====#
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (dir.exists(path)) {
          templ_files <- c("model.cpp",
                           "model.h",
                           "template_params.R",
                           "template_info.Rd")
          for (f in templ_files) {
            from_file <- file.path(template_models_dir(), name, f)
            overwrite(from_file, file.path(path, f))
          }
          for (f in c("template_augm.R")) {
            from_file <- file.path(template_models_dir(),
                                   name, f)
            if (file.exists(from_file)) {
              overwrite(from_file, file.path(path, f))
            }
          }

        } else {
          stop_h(
            sprintf("Failed to create directory \"%s\".", path),
            "Check your user permissions."
          )
        }
      }

    } else {
      stop_h(
        sprintf("Provided model name \"%s\" is unknown.", name),
        "Run `template_household_model` to see available template names."
      )
    }

  } else {
    stop_h("Model name is not a character.")
  }
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @description
#' Find the vector of definition names required by the model according.
#'
#' @param refer_file `character(1)`\cr
#' path/to/.refer containing the model source.
#' @param regular_expr `character()`\cr
#' Vector containing regular expressions.
#'
#' @return
#' `character()`)
#'
#' @noRd
template_find_requirement <- function(refer_file,
                                      regular_expr) {

  # check if reference file exists
  if (file.exists(refer_file)) {

    # read the script and extract parameters
    script <- readChar(refer_file, file.info(refer_file)$size)

    # extract all calls in .refer source file
    required_set <- c()
    for (reg in regular_expr) {
      matches <- stringr::str_match_all(script, reg)
      required_set <- c(required_set, unlist(matches[[1]][, 2]))
    }

    # re-order alphabetically
    if (!is.null(required_set)) {
      required_set <- unique(required_set[order(required_set)])
    }

  } else {
    warning_h(
      sprintf("Impossible to find reference file \"%s\".", refer_file),
      "Try to generate a new model object."
    )

    required_set <- NULL
  }

  return(required_set)
}

#' @description
#' Compared required definitions to current definitions and find missing ones.
#'
#' @param required_set `character()`\cr
#' Vector with necessary definitions.
#' @param defined_set `character()`\cr
#' Vector with current definitions.
#'
#' @return
#' `character()`)
#'
#' @noRd
template_find_missing <- function(required_set,
                                  defined_set) {
  if (length(required_set) == 0) {
    missing_set <- NULL

  } else if (length(defined_set) == 0) {
    missing_set <- unlist(required_set)

  } else {
    # unique and ordered set of parameters
    required_set <- unique(required_set[order(required_set)])
    defined_set <- unique(defined_set[order(defined_set)])

    missing_set <- c()
    for (name in required_set) {
      if (!(name %in% defined_set)) {
        missing_set <- c(missing_set, name)
      }
    }
  }

  return(missing_set)
}

#' @description
#' Compared required definitions to current definitions and find unused ones.
#'
#' @param required_set `character()`\cr
#' Vector with necessary definitions.
#' @param regular_expr `character()`\cr
#' Vector with current definitions.
#'
#' @return
#' `character()`)
#'
#' @noRd
template_find_unused <- function(required_set,
                                 defined_set) {
  if (length(required_set) == 0) {
    unused_set <- unlist(defined_set)

  } else if (length(defined_set) == 0) {
    unused_set <- NULL

  } else {
    # unique and ordered set of parameters
    required_set <- unique(required_set[order(required_set)])
    defined_set <- unique(defined_set[order(defined_set)])
    
    unused_set <- c()
    for (name in defined_set) {
      if (!(name %in% required_set)) {
        unused_set <- c(unused_set, name)
      }
    }
  }

  return(unused_set)
}


