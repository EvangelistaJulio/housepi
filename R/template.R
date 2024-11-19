#
# template.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title household model template
#' @description
#' Copy a default model to a given directory.
#' @param which (`character(1)`)\cr
#' Name of template model. If `NULL`, list of available templates.
#' @param path (`character(1)`)\cr
#' Where to copy the template scripts.
#' If `NULL`, show template model description in RStudio Help window.
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
