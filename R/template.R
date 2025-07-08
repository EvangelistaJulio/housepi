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

#' @title
#' Template models directory
#' @description
#' Returns the directory to housepi template models.
#'
#' @return `character(1)`
#' @keywords internal
template_models_dir <- function() {
  return(file.path(system.file(package = "housepi"), "models"))
}

#' @title
#' Show model names
#' @description
#' Show list of all available template nmodel names.
#'
#' @keywords internal
show_model_names <- function() {
  cli::cli_alert_info("Available template models:")
  ulid <- cli::cli_ul()
  for (name in list.files(template_models_dir())) cli::cli_li(name)
  cli::cli_end(ulid)
  invisible()
}

#' @title
#' Copy Template to Path
#' @description
#' Copy a default model to a given directory.
#'
#' @param name `character(1)`\cr
#' Name of template model. If `NA`, list of available templates.
#' @param path `character(1)`\cr
#' Where to copy the template scripts.
#' If `NA`, show template model description in RStudio Help window.
#'
#' @keywords internal
copy_template_to_path <- function(name, path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(path)) {
    templ_files <- c("model.cpp",
                     "model.h",
                     "template_param.R",
                     "template_info.Rd")
    for (f in templ_files) {
      from_file <- file.path(template_models_dir(), name, f)
      overwrite(from_file, file.path(path, f))
    }
    for (f in c("template_augm.R")) {
      from_file <- file.path(template_models_dir(), name, f)
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

#' @title
#' Template Household Model
#' @description
#' Navigate through housepi templates of household models.
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
template_household_model <- function(name = NA, path = NA, replace = FALSE) {

  if (is.na(name)) {
    #===== Asking for model names =====#
    show_model_names()

  } else if (is.character(name)) {

    if (name %in% list.files(template_models_dir())) {
      if (is.na(path)) {
        #===== Asking for model description =====#
        info_file <- file.path(template_models_dir(),
                               name, "template_info.Rd")
        if (!file.exists(info_file)) stop_h(
          "Impossible to show model description",
          "No template info were foun for the chosen model."
        )
        rstudioapi::previewRd(info_file)

      } else if (any(file.exists(source_files(path))) && isFALSE(replace)) {
        #===== Asking for copying template to folder with model =====#
        warning_h(
          "A model already exists in the specified folder.",
          "Remove its content before creating a new model."
        )

      } else {
        #===== Asking for copying template to empty folder =====#
        copy_template_to_path(name, path)
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

#' @title
#' Template assert expression
#' @description
#' Assert if a given string expression is present in source file.
#'
#' @param refer_file `character(1)`\cr
#' path/to/.refer containing the model source.
#' @param expr `character(1)`\cr
#' Expression.
#'
#' @return
#' `logical(1)`
#'
#' @keywords internal
template_assert_expression <- function(refer_file, expression) {
  # check if reference file exists
  if (file.exists(refer_file)) {
    # read the script and extract parameters
    script <- readChar(refer_file, file.info(refer_file)$size)

    # extract all calls in .refer source file
    matches <- stringr::str_match_all(script, expression)
    return(length(unlist(matches)) > 0)
  }
  return(FALSE)
}

#' @title
#' Template find requirement
#' @description
#' Find the vector of definition names required by the model according.
#'
#' @param refer_file `character(1)`\cr
#' path/to/.refer containing the model source.
#' @param regular_expr `character()`\cr
#' Vector containing regular expressions.
#'
#' @return
#' `character()`
#'
#' @keywords internal
template_find_requirement <- function(refer_file, regular_expr) {
  required <- c()

  # check if reference file exists
  if (file.exists(refer_file)) {
    # read the script and extract parameters
    script <- readChar(refer_file, file.info(refer_file)$size)

    # extract all calls in .refer source file
    for (reg in regular_expr) {
      matches <- stringr::str_match_all(script, reg)
      required <- c(required, unlist(matches[[1]][, 2]))
    }

    # re-order alphabetically
    if (length(required) > 0) required <- unique(required[order(required)])

  } else {
    warning_h(
      sprintf("Impossible to find reference file \"%s\".", refer_file),
      "Try to generate a new model object."
    )
  }

  return(required)
}

#' @title
#' Template find missing
#' @description
#' Compared required definitions to current definitions and find missing ones.
#'
#' @param required `character()`\cr
#' Vector with required definitions.
#' @param defined `character()`\cr
#' Vector with current definitions.
#'
#' @return
#' `character()`)
#'
#' @keywords internal
template_find_missing <- function(required, defined) {
  if (length(required) == 0) {
    miss <- NULL

  } else if (length(defined) == 0) {
    miss <- unlist(required)

  } else {
    # unique and ordered set of parameters
    required <- unique(required[order(required)])
    defined <- unique(defined[order(defined)])

    miss <- c()
    for (name in required) {
      if (!(name %in% defined)) miss <- c(miss, name)
    }
  }

  return(miss)
}

#' @title
#' Template Find Unused
#' @description
#' Compared required definitions to current definitions and find not used ones.
#'
#' @param required `character()`\cr
#' Vector with required definitions.
#' @param defined `character()`\cr
#' Vector with current definitions.
#'
#' @return
#' `character()`
#'
#' @keywords internal
template_find_unused <- function(required, defined) {
  if (length(required) == 0) {
    notused <- unlist(defined)

  } else if (length(defined) == 0) {
    notused <- NULL

  } else {
    # unique and ordered set of parameters
    required <- unique(required[order(required)])
    defined <- unique(defined[order(defined)])

    notused <- c()
    for (name in defined) {
      if (!(name %in% required)) notused <- c(notused, name)
    }
  }

  return(notused)
}

#' @title
#' Template Blank Dataset
#' @description
#' Template dataset associated to 'blank' model.
#'
#' @return `data.frame`
#'
#' @keywords internal
template_blank_dataset_ <- function() {
  data.frame(id = c("hh_01_1", "hh_01_2"), 
             house_id = "hh_01", 
             status = c(1, 0))
}