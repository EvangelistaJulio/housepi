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

# Assert that given definition name is valid.
assert_valid_definition_name <- function(x, family) {
  if (is.na(x) || !is.character(x) || length(x) != 1 || nchar(x) == 0) stop_h(
    sprintf("Definition name of %s must be a non-empty character(1).", family)
  )
  else if (is.na(x)) stop_h(
    sprintf("Incomplete definition of %s: missing name.", family)
  )
  else if (!grepl("^[A-Za-z0-9_.]+$", x)) stop_h(
    sprintf("Definition name '%s' of %s contains invalid characters.",
            x, family),
    "Only letters, numbers, underscores (_) and dots (.) are allowed."
  )
}

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
#' @param x
#' A `character(1)` with the choice of the template.
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
  NA_character_
}

#' @title
#' Check if a template with similar name exists
#' @description
#' Return the directory of chosen housepi template model.
#'
#' @param x
#' A `character(1)` with the choice of the template.
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
  model_names <- list.files(template_models_dir())
  formatted_names <- tolower(gsub("-|_", " ", model_names))
  message(sprintf(
    "Available template models:\n%s\n",
    paste(paste0(cli::symbol$bullet, " ", formatted_names), collapse = "\n")
  ))
  invisible()
}

#' @title
#' Copy Template to Path
#' @description
#' Copy a default model to a given directory.
#'
#' @param name
#' A `character(1)` with the name of template model.
#' If `NA` returns list of available templates.
#' @param path
#' A `character(1)` with the path to where to copy the template scripts.
#' If `NA` shows template model description in RStudio Help window.
#'
#' @keywords internal
#' @noRd
copy_template_to_path <- function(name, path) {
  # avoid current directory
  if (path == ".") {
    path <- name
    warning_h(
      "Using current directory is not allowed for copying a template.",
      sprintf("Creating folder '%s' instead.", path)
    )
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(path)) {
    templ_files <- c("model.cpp",
                     "model.h",
                     "template_param.R",
                     "template_info.Rd")
    for (f in templ_files) {
      from_file <- file.path(get_template_dir(name), f)
      overwrite(from_file, file.path(path, f))
    }
    for (f in c("template_augm.R")) {
      from_file <- file.path(get_template_dir(name), f)
      if (file.exists(from_file)) {
        overwrite(from_file, file.path(path, f))
      }
    }

  } else {
    stop_h(
      sprintf("Failed to create directory '%s'.", path),
      "Check your user permissions."
    )
  }
}

#' @title
#' Open Model Description
#' @description
#' Preview the model description in RStudio Help window.
#'
#' @param choice
#' A `character(1)` with the path/to/model folder.
#'
#' @noRd
open_model_helper_ <- function(path) {
  # check for info file
  info_file <- file.path(path, "template_info.Rd")
  if (!file.exists(info_file)) stop_h(
    "Impossible to show model description.",
    "No template info available for the given model."
  )

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    # print in RStudio Help window
    rstudioapi::previewRd(info_file)
  } else if (interactive()) {
    # print as plain text in the console
    cat(utils::capture.output(tools::Rd2txt(info_file)), sep = "\n")
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
#' @param name
#' A `character(1)` with the name of template model.
#' If `NA` returns list of available templates.
#' @param path
#' A `character(1)` with the path to where to copy the template scripts.
#' If `NA` and name is not, shows template description in RStudio Help window.
#' @param replace
#' `logical(1)` to replace previous content of the path.
#'
#' @export
template_household_model <- function(name = NA_character_,
                                     path = NA_character_,
                                     replace = FALSE) {
  if (is.na(name)) show_model_names()
  else if (!is.character(name)) stop_h("Model name is not a character.")
  else if (is.na(get_template_dir(name))) stop_h(
    sprintf("Model '%s' was not a known template.", name),
    "Show available templates calling 'template_household_model'."
  )
  # asking for template model description
  else if (is.na(path)) open_model_helper_(get_template_dir(name))
  # asking for copying template to folder with model
  else if (all(file.exists(source_files(path)), !replace)) warning_h(
    "A model already exists in the specified folder.",
    "If you want to replace it, run again with option 'replace = TRUE'."
  )
  # askinf for copying template to empty folder
  else copy_template_to_path(name, path)
}

############################################################################
############################################################################
############################################################################

# Associate same name to all elements of a vector
set_source_name_ <- function(x, name = NA_character_) {
  if (!is.na(name) &&
        is.character(name) &&
        nchar(name) > 0 &&
        !is.null(x)) names(x) <- rep(name, length(x))
  x
}

#' @title
#' Template assert expression
#' @description
#' Assert if a given string expression is present in source file.
#'
#' @param refer_file
#' A `character(1)` with the path/to/.refer containing the model source.
#' @param expr
#' A `character(1)` with the expression.
#'
#' @return `logical(1)`
#'
#' @keywords internal
#' @noRd
template_assert_expression <- function(refer_file, expression) {
  # check if reference file exists
  if (!file.exists(refer_file)) return(FALSE)
  # read the script and extract parameters
  script <- readChar(refer_file, file.info(refer_file)$size)

  # extract all calls in .refer source file
  matches <- stringr::str_match_all(script, expression)
  length(unlist(matches)) > 0
}

#' @title
#' Template find requirement
#' @description
#' Find the vector of definition names required by the model according.
#'
#' @param refer_file
#' A `character(1)` with the path/to/.refer containing the model source.
#' @param regular_expr
#' A `character()` of regular expressions.
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
      sprintf("Impossible to find reference file '%s'.", refer_file),
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
#' @param required
#' A `character()` with required definitions.
#' @param defined
#' A `character()` with current definitions.
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
#' @param required `character()` with required definitions.
#' @param defined `character()` with current definitions.
#' @param family `character(1)` with the name of the definition family.
#' @param max_print `integer(1)` with the maximum number of elements to print.
#'
#' @return
#' `character()`
#'
#' @keywords internal
#' @noRd
template_missing_stop <- function(required, defined, family, max_print = 4) {
  miss <- template_missing_find(required, defined)
  len_miss <- length(miss)
  # find sources of missing definitions
  sources <- unique(names(required[required %in% miss]))
  # show error message
  output <- c()
  if (len_miss > 0) output <- c(
    sprintf("There %s %s missing %s definition%s: %s",
            ifelse(len_miss > 1, "were", "was"),
            ifelse(len_miss > max_print,
                   paste0("more than ", max_print), len_miss),
            family,
            ifelse(len_miss > 1, "s", ""),
            character_to_string_list(miss, max_print)),
    sprintf("The %s definition%s required by %s",
            family,
            ifelse(len_miss > 1, "s were", " was"),
            character_to_phrase(sources))
  )
  output
}

############################################################################
############################################################################
############################################################################

#' @title
#' Template Find Not Used
#' @description
#' Compared required definitions to current definitions and find unused ones.
#'
#' @param required
#' A `character()` with required definitions.
#' @param defined
#' A `character()` with current definitions.
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
#' @param required A `character()` with required definitions.
#' @param defined A `character()` with current definitions.
#' @param max_print `integer(1)` with the maximum number of elements to print.
#'
#' @return
#' `character()`
#'
#' @keywords internal
#' @noRd
template_not_used_warn <- function(required, defined, family, max_print = 5) {
  notused <- template_not_used_find(required, defined)

  len_not_used <- length(notused)
  if (len_not_used > 0) warning_h(
    sprintf("There %s %s unused %s definition%s%s",
            ifelse(len_not_used > 1, "were", "was"),
            ifelse(len_not_used > max_print,
                   paste0("more than ", max_print), len_not_used),
            family,
            ifelse(len_not_used > 1, "s:\n", ": "),
            character_to_string_list(notused, max_print))
  )
  invisible()
}

#' @title
#' Cli format assertion Not Used Warning
#' @description
#' Format the output of an assert message with CLI.
#'
#' @keywords internal
#' @noRd
cli_format_assertion <- function(sms, section) {
  if (length(sms) > 0) {
    sms <- c(paste0(cli::col_yellow(section), ":"), sms)
    sms <- gsub(sms, pattern = " were ", replacement = " are ")
    sms <- paste0(sms, collapse = "\n")
    cli::cli_alert_danger(sms)
  } else {
    cli::cli_alert_success(paste0(section, " definitions are complete."))
  }
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
