#
# template.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @description
#' Show list of all available template nmodel names.
#' @noRd
show_model_names <- function() {

  cli::cli_alert_info("Template model names:")
  ulid <- cli::cli_ul()

  for (name in list.files(template_models_dir())) {
    cli::cli_li(name)
  }
  cli::cli_end(ulid)

  invisible()
}

#' @title household model template
#' @description
#' Copy a default model to a given directory.
#' @param template_name (`character(1)`)\cr
#' Name of template model. If `NULL`, list of available templates.
#' @param path (`character(1)`)\cr
#' Where to copy the template scripts.
#' If `NULL`, show template model description in RStudio Help window.
#' @export
household_model_template <- function(template_name = NULL, path = NULL) {

  if (is.null(template_name)) {
    #===== Asking for model names =====#
    show_model_names()

  } else if (is.character(template_name)) {

    if (template_name %in% list.files(template_models_dir())) {
      if (is.null(path)) {
        #===== Asking for model description =====#
        info_file <- file.path(template_models_dir(),
                               template_name,
                               "template_info.Rd")
        if (!file.exists((info_file))) {
          stop(cli::format_error(
            "Impossible to show model description: no template info were found."
          ))
        }

        # FIXME: this doesn't work outside of Rstudio
        rstudioapi::previewRd(info_file)

      } else if (any(file.exists(source_files(path)))) {
        #===== Asking for copying template to folder with model =====#
        cli::cli_alert_warning(cli::format_warning(c(
          sprintf("A model already exists in folder \"%s\".", path),
          i = "To create a model in the same folder, delete its content first."
        )))

      } else {
        #===== Asking for copying template to empty folder =====#
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (dir.exists(path)) {
          # UPDATE : maybe don't copy the copyright section

          for (f in c("model.cpp", "model.h", "template_params.R", "template_info.Rd")) {
            from_file <- file.path(template_models_dir(),
                                   template_name, f)
            if (!overwrite(from_file, file.path(path, f))) {
              internal_setup_error()
            }
          }
          for (f in c("template_augm.R")) {
            from_file <- file.path(template_models_dir(),
                                   template_name, f)
            overwrite(from_file, file.path(path, f))
          }

        } else {
          internal_setup_error()
        }
      }

    } else {
      stop(cli::format_error(c(
        sprintf("Provided model name \"%s\" is unknown.", template_name),
        i = "Run {.code household_model_template()} to see available template names."
      )))
    }

  } else {
    stop(cli::format_error("Model name is not a character."))
  }
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @description
#' Assert that all required definition names are present.
#' @param refer_file (`character(1)`)\cr
#' path/to/.refer containing the model source.
#' @param regular_expr (`character()`)\cr
#' Vector containing regular expressions.
#' @return (`character()`)\cr
#' Vector with necessary definitions.
#' @noRd
template_requirement <- function(refer_file,
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
    required_set <- unique(required_set[order(required_set)])

  } else {
    cli::cli_alert_warning(cli::format_warning(c(
      sprintf("Impossible to find reference file \"%s\".", refer_file),
      i = "Try to generate a new model object with {.code force = TRUE}."
    )))

    required_set <- NULL
  } 
  
  return(required_set)
}

#' @description
#' Check for missing definitions.
#' @param required_set (`character()`)\cr
#' Vector with necessary definitions.
#' @param regular_expr (`character()`)\cr
#' Vector with current definitions.
#' @return (`character()`)\cr
#' Vector with missing definitions.
#' @noRd
template_missing <- function(required_set,
                            def_set) {
  if (is.null(required_set)) {
    missing_set <- NULL

  } else {
    missing_set <- c()
    for (name in required_set) {
      if (!(name %in% def_set)) {
        missing_set <- c(missing_set, name)
      }
    }
  }

  return(missing_set)
}

#' @description
#' Check for unused definitions.
#' @param required_set (`character()`)\cr
#' Vector with necessary definitions.
#' @param regular_expr (`character()`)\cr
#' Vector with current definitions.
#' @return (`character()`)\cr
#' Vector with unused definitions.
#' @noRd
template_unused <- function(required_set,
                            def_set) {
  if (is.null(required_set)) {
    unused_set <- NULL

  } else {
    unused_set <- c() 
    for (name in def_set) {
      if (!(name %in% required_set)) {
        unused_set <- c(unused_set, name)
      }
    }
  }

  return(unused_set)
}

character_to_string_list <- function(mylist) {
  return(paste(lapply(mylist, function(x) { paste0("\"", x ,"\"") } ),
               collapse=", "))
}
