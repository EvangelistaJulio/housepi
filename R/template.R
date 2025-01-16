#
# template.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# @description
# Show list of all available template nmodel names.
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
