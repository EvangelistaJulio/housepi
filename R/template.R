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

spaced_template_choice <- function(x) {
  tolower(gsub("\\s+", " ", gsub("-|_| ", " ", x)))
}

# Returns the directory to all housepi template models.
template_models_dir <- function() {
  system.file("models", package = "housepi")
}

#' @title
#' Get template model directory
#' @description
#' Return the directory of chosen housepi template model.
#'
#' @param x `character(1)`\cr
#' Choice of the template.
#'
#' @return `character(1)`
#' @keywords internal
#' @noRd
get_template_dir <- function(x) {
  flat_choice <- spaced_template_choice(x)
  path_to_options <- list.files(template_models_dir(), full.names = TRUE)
  flat_options <- spaced_template_choice(basename(path_to_options))
  # extract path if present (unique since folder name)
  mask <- flat_options == flat_choice
  if (any(mask)) return(path_to_options[mask])
  NA
}

#' @title
#' Check if a template with similar name exists
#' @description
#' Return the directory of chosen housepi template model.
#'
#' @param x `character(1)`\cr
#' Choice of the template.
#'
#' @return `character(1)`
#' @keywords internal
#' @noRd
find_closest_template_name <- function(x) {
  flat_choice <- spaced_template_choice(x)
  options <- list.files(template_models_dir(), full.names = FALSE)
  flat_options <- spaced_template_choice(options)
  sms <- find_closest_match(flat_choice, flat_options, 5)
  if (stringr::str_ends(sms, "\\?")) return(sms)
  c("Show available templates calling 'template_household_model'.")
}

#' @title
#' Show model names
#' @description
#' Show list of all available template nmodel names.
#'
#' @keywords internal
#' @noRd
show_model_names <- function() {
  cli::cli_alert_info("Available template models:")
  ulid <- cli::cli_ul()
  for (name in list.files(template_models_dir())) {
    cli::cli_li(tolower(gsub("-|_", " ", name)))
  }
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
#' @noRd
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
#' Open Model Description
#' @description
#' Preview the model description in RStudio Help window.
#'
#' @param choice `character(1)`\cr
#' Path/to/model folder.
#'
#' @noRd
open_model_helper_ <- function(path) {
  # check for info file
  info_file <- file.path(path, "template_info.Rd")
  if (!file.exists(info_file)) stop_h(
    "Impossible to show model description",
    "No template info available for the given model."
  )

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    # print in RStudio Help window
    rstudioapi::previewRd(info_file)
  } else if (interactive()) {
    # print as plain text in the console
    cat(capture.output(tools::Rd2txt(info_file)), sep = "\n")
  } else {
    warning_h("Help method not supported in non-interactive mode.")
  }
  invisible()
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
    # asking for model names
    show_model_names()

  } else if (is.character(name)) {
    # asking for description or for copying template

    # if name is a template...
    if (name %in% list.files(template_models_dir())) {
      if (is.na(path)) {
        # ... ask for template model description
        open_model_helper_(file.path(template_models_dir(), name))
      } else if (any(file.exists(source_files(path))) && isFALSE(replace)) {
        # ... ask for copying template to folder with model
        warning_h(
          "A model already exists in the specified folder.",
          "Remove its content before creating a new model."
        )
      } else {
        # ... ask for copying template to empty folder
        copy_template_to_path(name, path)
      }

    } else {
      stop_h(
        sprintf("Provided model template \"%s\" is unknown.", name),
        "Show available templates calling 'template_household_model'."
      )
    }

  } else {
    stop_h("Model name is not a character.")
  }
}

############################################################################
############################################################################
############################################################################

# Associate same name to all elements of a vector
set_source_name_ <- function(x, name = NA) {
  if (!is.na(name) && is.character(name) && nchar(name) > 0) {
    if (!is.null(x)) x <- setNames(x, rep(name, length(x)))
  }
  return(x)
}

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
#' @noRd
template_assert_expression <- function(refer_file, expression) {
  # check if reference file exists
  if (file.exists(refer_file)) {
    # read the script and extract parameters
    script <- readChar(refer_file, file.info(refer_file)$size)

    # extract all calls in .refer source file
    matches <- stringr::str_match_all(script, expression)
    return(length(unlist(matches)) > 0)
  }
  FALSE
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
#' @noRd
template_required_find <- function(refer_file, regular_expr) {
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
    required <- sort(unique(required))

  } else {
    warning_h(
      sprintf("Impossible to find reference file \"%s\".", refer_file),
      "Try to generate a new model object."
    )
  }

  required
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
#' @noRd
template_missing_find <- function(required, defined) {
  if (length(required) == 0) NULL
  else if (length(defined) == 0) sort(unlist(required))
  else sort(setdiff(unlist(required), unlist(defined)))
}

#' @title
#' Template Missing Stop
#' @description
#' Stop if missing definitions are found
#'
#' @param required `character()`\cr
#' Vector with required definitions.
#' @param defined `character()`\cr
#' Vector with current definitions.
#' @param family `character(1)`\cr
#' Name of the definition family.
#'
#' @return
#' `character()`
#'
#' @keywords internal
#' @noRd
template_missing_stop <- function(required, defined, family) {
  miss <- template_missing_find(required, defined)
  len_miss <- length(miss)
  # find sources of missing definitions
  sources <- unique(names(required[required %in% miss]))
  # show error message
  if (len_miss > 0) stop_h(
    sprintf("There %s %s missing %s definition%s: %s",
            ifelse(len_miss > 1, "are", "is"), len_miss,
            family,
            ifelse(len_miss > 1, "s", ""),
            character_to_string_list(miss)),
    sprintf("The %s definition%s required by %s",
            family,
            ifelse(len_miss > 1, "s are", " is"),
            character_to_phrase(sources))
  )
  invisible()
}

############################################################################
############################################################################
############################################################################

#' @title
#' Template Find Not Used
#' @description
#' Compared required definitions to current definitions and find unused ones.
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
#' @noRd
template_not_used_find <- function(required, defined) {
  if (length(required) == 0) sort(unlist(defined))
  else if (length(defined) == 0) NULL
  else sort(setdiff(unlist(defined), unlist(required)))
}

#' @title
#' Template Not Used Warning
#' @description
#' Raise a warning if unused definitions are found
#'
#' @param required `character()`\cr
#' Vector with required definitions.
#' @param defined `character()`\cr
#' Vector with current definitions.
#' @param family `character(1)`\cr
#' Name of the definition family.
#'
#' @return
#' `character()`
#'
#' @keywords internal
#' @noRd
template_not_used_warn <- function(required, defined, family) {
  notused <- template_not_used_find(required, defined)

  len_not_used <- length(notused)
  if (len_not_used > 0) warning_h(
    sprintf("There %s %s unused %s definition%s: %s",
            ifelse(len_not_used > 1, "are", "is"),
            len_not_used, family,
            ifelse(len_not_used > 1, "s", ""),
            character_to_string_list(notused))
  )
  invisible()
}

############################################################################
############################################################################
############################################################################

#' @title
#' Template Blank Dataset
#' @description
#' Template data_wide associated to 'blank' model.
#'
#' @return `data.frame`
#'
#' @keywords internal
#' @noRd
template_blank_data_wide_ <- function() {
  data.frame(id = c("hh_01_1", "hh_01_2"),
             house_id = "hh_01",
             status = c(1, 0))
}
