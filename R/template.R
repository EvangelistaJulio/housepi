## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom Rcpp Module

# Associate dynamic library extension to name
dll_file_ext <- function(name) {
  paste0(name, .Platform$dynlib.ext)
}

# Sidecar file recording the compile flags baked into a dll.
dll_flags_file <- function(dll) {
  paste0(dll, ".housepi.flags")
}

# Templated name convention (remove punctuation).
template_name <- function(path) {
  gsub("[[:punct:] ]", "_", basename(path))
}

# Add random string to templated name
random_name <- function(path) {
  # generate a random string
  str_len <- 8
  str_name <- paste0(sample(c(letters, LETTERS, 0:9), str_len, replace = TRUE),
                     collapse = "")
  paste0(template_name(path), str_name)
}

# DLL file with templated name convention
template_dll_file <- function(path) {
  file.path(path, dll_file_ext(template_name(path)))
}

# Return template name choice with spaces.
spaced_template_choice <- function(x) {
  trimws(tolower(gsub("\\s+", " ", gsub("-|_| ", " ", x))))
}

# Returns the directory to all housepi template models.
template_models_dir <- function() {
  system.file("models", package = "housepi")
}

# Get template name from DLL file
get_template_from_dll <- function(dll) {
  if (is.null(dll) || is.na(dll)) NA_character_
  else tools::file_path_sans_ext(basename(dll))
}

# Check if the same DLL is already loaded.
is_dll_loaded <- function(dll) {
  any(vapply(getLoadedDLLs(), function(x) x[["path"]] == dll, FALSE))
}

# Load model DLL if not already loaded.
# Auto-registration in the DLL's static init populates the ModelRegistry.
load_model_dll_ <- function(dll) {
  if (!is.null(dll) && !is.na(dll) && !is_dll_loaded(dll)) {
    dyn.load(dll)
  }
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
#' Suggest the closest matching template name
#' @description
#' Return a "did you mean" hint for a mistyped template name.
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
  # NOTE: avoid to attempt a suggestion for silly mistakes
  sms <- find_closest_match(flat_choice, flat_options)
  if (stringr::str_ends(sms, "\\?") && nchar(flat_choice) > 3) return(sms)
  c("Show available templates calling 'browse_housepi_templates'.")
}

#' @title
#' Show available template model names
#' @description
#' Show list of all available template model names.
#'
#' @keywords internal
#' @noRd
show_model_names <- function() {
  model_names <- list.files(template_models_dir())
  formatted_names <- tolower(gsub("-| |_", "_", model_names))
  sms <- c(i = "Available template models:")
  for (m in formatted_names) {
    link <- sprintf(
      "{.run [%s](housepi::browse_housepi_templates(\"%s\"))}", m, m
    )
    sms <- c(sms, paste0(cli::symbol$bullet, " ", link))
  }
  cli::cli_inform(sms)
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
    warn_h(
      "Using current directory is not allowed for copying a template.",
      sprintf("Creating folder '%s' instead.", path)
    )
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(path)) {
    # model.cpp + info.Rd are always present; param.R / augm.R are optional
    required <- c("model.cpp", "info.Rd")
    optional <- c("param.R", "augm.R")
    for (f in c(required, optional)) {
      from_file <- file.path(get_template_dir(name), f)
      if (file.exists(from_file)) {
        overwrite(from_file, file.path(path, f))
      } else if (f %in% required) {
        stop_h(sprintf("Template '%s' is missing required file '%s'.",
                       name, f))
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
  info_file <- file.path(path, "info.Rd")
  if (!file.exists(info_file)) stop_h(
    "Impossible to show model description.",
    "No template info available for the given model."
  )

  # load system + Rdpack macros so \doi and \insertRef{key}{housepi}
  # resolve against inst/REFERENCES.bib at render time
  sys_macros <- tools::loadRdMacros(
    file.path(R.home("share"), "Rd", "macros", "system.Rd")
  )
  macros <- tools::loadRdMacros(
    system.file("help/macros/refmacros.Rd", package = "Rdpack"),
    sys_macros
  )
  rd <- tools::parse_Rd(info_file, macros = macros)
  stages <- c("build", "install", "render")

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
    tmp <- tempfile(fileext = ".html")
    tools::Rd2HTML(rd, out = tmp, stages = stages)
    rstudioapi::viewer(tmp)
  } else if (interactive()) {
    cat(utils::capture.output(tools::Rd2txt(rd, stages = stages)), sep = "\n")
  } else {
    warn_h("Help method not supported in non-interactive mode.")
  }
  invisible()
}

#' @title Navigate through household model templates
#' @description
#' Print the list of available templates,
#' get their descriptions or copy source code to custom path.
#'
#' @param name
#' A `character(1)` with the name of template model.
#' If `NA` print the list of available templates (default).
#' @param path
#' A `character(1)` with the path to where to copy the template scripts.
#' If `NA` and name is not, print template description in RStudio Help window.
#' @param replace
#' `logical(1)` to replace previous content of the path.
#'
#' @examples
#' # list all available templates
#' browse_housepi_templates()
#' @export
browse_housepi_templates <- function(name = NA_character_,
                                     path = NA_character_,
                                     replace = FALSE) {
  if (is.na(name)) show_model_names()
  else if (!is.character(name)) stop_h("Model name is not a character.")
  else if (is.na(get_template_dir(name))) stop_h(
    sprintf("Model '%s' was not a known template.", name),
    find_closest_template_name(name)
  )
  # asking for template model description
  else if (is.na(path)) open_model_helper_(get_template_dir(name))
  # asking for copying template to folder with model
  else if (file.exists(file.path(path, "model.cpp")) && !replace) warn_h(
    "A model already exists in the specified folder.",
    "If you want to replace it, run again with option 'replace = TRUE'."
  )
  # asking for copying template to empty folder
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

#' Extract arguments following a key family in REGISTER_HOUSEPI_KEYS.
#' @param refer_file `character(1)`. Path to reference file.
#' @param key_family `character(1)`. Key family name to search for.
#' @return `character()`
#' @keywords internal
#' @noRd
template_extract_keys <- function(refer_file, key_family) {
  script <- paste(readLines(refer_file), collapse = " ")
  reg <- sprintf(
    "HOUSEPI_DECLARE_REGISTER_KEYS\\(%s,\\s*([^)]+)\\)",
    key_family
  )
  args <- stringr::str_match(script, reg)[1, 2]
  if (is.na(args)) return(character(0))
  args <- stringr::str_trim(stringr::str_split(args, ",")[[1]])
  args[nzchar(args)]
}

#' Compare required definitions to current definitions and find missing ones.
#' @param required `character()`. Required definitions.
#' @param defined `character()`. Current definitions.
#' @return `character()`
#' @keywords internal
#' @noRd
template_missing_find <- function(required, defined) {
  if (length(required) == 0) NULL
  else if (length(defined) == 0) sort(unlist(required))
  else sort(setdiff(unlist(required), unlist(defined)))
}

#' Stop if missing definitions are found
#' @param required `character()`. Required definitions.
#' @param defined `character()`. Current definitions.
#' @param family `character(1)`. Name of the definition family.
#' @param max_print `integer(1)`. Maximum number of elements to print.
#' @return `character()`
#' @keywords internal
#' @noRd
template_missing_stop <- function(required, defined, family, max_print = 4) {
  miss <- template_missing_find(required, defined)
  # find sources of missing definitions
  sources <- unique(names(required[required %in% miss]))
  # show error message
  output <- c()
  len_miss <- length(miss)
  if (len_miss > 0) output <- c(
    sprintf("%s missing %s %s%s",
            format_there_were_was(len_miss, max_print = max_print),
            family,
            as_plural("definition", len_miss, "s:"),
            character_to_string_list(miss, max_print = max_print)),
    sprintf("%s required by %s definitions.",
            as_plural("", len_miss, "They were", "It was"),
            character_to_phrase(sources))
  )
  output
}

#' Warn about overwritten definitions if missing definitions are found.
#' @param all_rmv `character()`. Removed definitions.
#' @param family `character(1)`. Name of the module family.
#' @param max_print `integer(1)`. Maximum number of elements to print.
#' @keywords internal
#' @noRd
template_overwrite_warn <- function(all_rmv, family, max_print = 4) {
  len_rmv <- length(all_rmv)
  if (len_rmv > 0) warn_h(sprintf(
    "Overwriting %s %s %s%s",
    ifelse(len_rmv > max_print,
           paste0("more than ", max_print), len_rmv),
    family,
    as_plural("definition", len_rmv, "s for:", " for"),
    character_to_string_list(all_rmv, max_print)
  ))
}

############################################################################
############################################################################
############################################################################

#' Compare required definitions to current definitions and find unused ones.
#' @param required `character()`. Required definitions.
#' @param defined `character()`. Current definitions.
#' @return `character()`
#' @keywords internal
#' @noRd
template_not_used_find <- function(required, defined) {
  if (length(required) == 0) sort(unlist(defined))
  else if (length(defined) == 0) NULL
  else sort(setdiff(unlist(defined), unlist(required)))
}

#' Raise a warning if unused definitions are found
#' @param required `character()`. Required definitions.
#' @param defined `character()`. Current definitions.
#' @param max_print `integer(1)`. Maximum number of elements to print.
#' @return `character()`
#' @keywords internal
#' @noRd
template_not_used_warn <- function(required, defined, family, max_print = 4) {
  notused <- template_not_used_find(required, defined)
  len_unused <- length(notused)
  if (len_unused > 0) warn_h(
    sprintf("%s unused %s %s%s",
            format_there_were_was(len_unused, max_print),
            family,
            as_plural("definition", len_unused, "s:"),
            character_to_string_list(notused, max_print))
  )
  invisible()
}

#' Format the output of an assert usage definition message with cli.
#' @param section `character(1)`. Section name.
#' @param sms `character()`. Messages to format.
#' @param minor `logical(1)`. Whether to format as minor message.
#' @keywords internal
#' @noRd
fmt_assert_cli <- function(section, sms, minor = FALSE) {
  sms <- gsub(sms, pattern = " were ", replacement = " are ")
  sms <- gsub(sms, pattern = " was ", replacement = " is ")
  lower_first <- function(col_sec, x) {
    paste0(col_sec, ": ", tolower(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  }

  if (length(sms) == 0) {
    c(v = paste0(cli::col_green(section), ": definitions complete."))
  } else if (minor) {
    sms <- gsub(sms, pattern = " (are|is) required by ",
                replacement = " can be specified for ")
    c(`!` = lower_first(cli::col_yellow(section), sms[1]), sms[-1])
  } else {
    c(x = lower_first(cli::col_red(section), sms[1]), sms[-1])
  }
}

#' Format a vector of messages for `cli::cli_alert_info()`.
#' @param sms `character()`. Messages, optionally named with `!`, `x`, `v`.
#' @keywords internal
#' @noRd
fmt_cli_info <- function(sms) {
  if (length(sms) == 0) return("")
  nms <- names(sms)
  if (is.null(nms)) nms <- rep("", length(sms))
  glyphs <- list(`!` = cli::col_yellow("!"),
                 `x` = cli::col_red(cli::symbol$cross),
                 `v` = cli::col_green(cli::symbol$tick))
  if (identical(nms[1], "i")) nms[1] <- ""
  out <- vapply(seq_along(sms), function(k) {
    g <- glyphs[[nms[k]]]
    if (is.null(g)) sms[[k]] else paste0(g, " ", sms[[k]])
  }, character(1))
  cli::cli_alert_info(paste(out, collapse = "\n"))
}
