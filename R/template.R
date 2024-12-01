#
# template.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

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
  # system path to housepi
  models_dir <- file.path(system.file(package = "housepi"), "models")

  if (is.null(template_name)) {
    cli::cli_alert_info("Template model names:")
    ulid <- cli::cli_ul()
    for (name in list.files(models_dir)) {
      cli::cli_li(name)
    }
    cli::cli_end(ulid)

  } else if (is.character(template_name)) {

    if (template_name %in% list.files(models_dir)) {
      if (is.null(path)) {
        #----- show model description -----#
        templ_info_file <- file.path(models_dir, template_name, "template_info.Rd")
        if (!file.exists((templ_info_file))) {
          stop(cli::format_error(
            "Impossible to show model description: no template info were found."
          ))
        }

        # FIXME: this doesn't work outside of Rstudio
        rstudioapi::previewRd(templ_info_file)

      } else if (any(file.exists(source_files(path)))) {
        cli::cli_alert_warning(cli::format_warning(c(
          sprintf("A model already exists in folder '%s'.", path),
          i = "To create a model in the same folder, delete its content first."
        )))
      } else {
        #----- copy files -----#
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (dir.exists(path)) {
          # UPDATE : maybe don't copy the copyright section
          
          for (f in c("model.cpp", "model.h", "template_params.R", "template_info.Rd")) {
            from_file <- file.path(models_dir, template_name, f)
            to_file <- file.path(path, f)
            if (!overwrite(from_file, to_file)) {
              stop(cli::format_error(c(
                x = "Failed to copy model files.",
                i = "Check your user permissions."
              )))
            }
          }
          for (f in c("template_augm.R")) {
            from_file <- file.path(models_dir, template_name, f)
            to_file <- file.path(path, f)
            overwrite(from_file, to_file)
          }

        } else {
          stop(cli::format_error(c(
            x = sprintf("Impossibile to create directory '%s'.", path),
            i = "Check your user permissions."
          )))
        }
      }

    } else {
      stop(cli::format_error(c(
        x = sprintf("Provided model name '%s' is unknown.", template_name),
        i = "To see the available template models, run: {.code template_model()}"
      )))
    }

  } else {
    stop(cli::format_error("Model name is not a character."))
  }
}
