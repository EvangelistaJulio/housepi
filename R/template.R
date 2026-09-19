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

# Templated name convention (non-identifier characters become underscores).
template_name <- function(path) {
  gsub("[^A-Za-z0-9_]", "_", basename(path))
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
# NOTE: auto-registration in the DLL's static init populates the model registry.
load_model_dll_ <- function(dll) {
  if (!is.null(dll) && !is.na(dll) && !is_dll_loaded(dll)) dyn.load(dll)
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
  c("Show available templates calling {.fn browse_templates}.")
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
      "{.run [%s](housepi::browse_templates(\"%s\"))}", m, m
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
    warn_h("Using current directory is not allowed for copying a template.",
           `!` = "Creating folder {.file {path}} instead.")
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) raise_permissions_error()

  # source cpp and info are always present
  required <- c("model.cpp", "info.Rd")
  # definition templates are optional
  optional <- unname(template_def_files)
  for (f in c(required, optional)) {
    from_file <- file.path(get_template_dir(name), f)
    if (file.exists(from_file)) overwrite(from_file, file.path(path, f))
    else if (f %in% required) stop_h(
      "Template {.val {name}} is missing required file {.file {f}}."
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
    `!` = "No template info available for the given model."
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
    # nocov start
    tmp <- tempfile(fileext = ".html")
    tools::Rd2HTML(rd, out = tmp, stages = stages)
    rstudioapi::viewer(tmp)
    # nocov end
  } else if (interactive()) {
    cat(utils::capture.output(tools::Rd2txt(rd, stages = stages)), sep = "\n")
  } else {
    warn_h("Help method not supported in non-interactive mode.")
  }
  invisible()
}

#' @title Navigate through household model templates
#' @description Print the list of available templates,
#' get their descriptions or copy source code to custom path.
#' @param name
#' A `character(1)` with the name of template model.
#' If `NA` print the list of available templates (default).
#' @param path
#' A `character(1)` with the path to where to copy the template scripts.
#' If `NA` and name is not, print template description in RStudio Help window.
#' @param replace
#' `logical(1)` to replace previous content of the path.
#' @return Invisibly returns `NULL`; called for printing the templates,
#' showing a description or copying the template files.
#' @examples
#' # list all available templates
#' browse_templates()
#' @export
browse_templates <- function(name = NA_character_,
                             path = NA_character_,
                             replace = FALSE) {
  if (is.na(name)) show_model_names()
  else if (!is.character(name)) stop_h("Model name is not a character.")
  else if (is.na(get_template_dir(name))) stop_h(
    "Model {.val {name}} was not a known template.",
    find_closest_template_name(name)
  )
  # asking for template model description
  else if (is.na(path)) open_model_helper_(get_template_dir(name))
  # asking for copying template to folder with model
  else if (file.exists(file.path(path, "model.cpp")) && !replace) warn_h(
    "A model already exists in the specified folder.",
    "To replace it, run again with option {.code replace = TRUE}."
  )
  # asking for copying template to empty folder
  else copy_template_to_path(name, path)
  invisible()
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

# Collapse definition sources into a phrase.
format_sources_phrase_ <- function(sources) {
  if (length(sources) == 0) return(NULL)
  character_to_phrase(ifelse(sources == "model", "transmission model",
                             paste(sources, "definitions")))
}

#' Compare required definitions to current definitions and find missing ones.
#' @param required `character()` Required definitions.
#' @param defined `character()` Current definitions.
#' @return `character()`
#' @keywords internal
#' @noRd
template_missing_find <- function(required, defined) {
  if (length(required) == 0) NULL
  else if (length(defined) == 0) sort(unlist(required))
  else sort(setdiff(unlist(required), unlist(defined)))
}

#' Stop if missing definitions are found
#' @param required `character()` Required definitions.
#' @param defined `character()` Current definitions.
#' @param holder `character(1)` Name of the definition holder.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @param noun `character(1)` Word naming the missing entries.
#' @param style `character(1)` cli markup class wrapping each entry.
#' @return `character()`
#' @keywords internal
#' @noRd
template_missing_stop <- function(required, defined, holder, max_print = 4,
                                  noun = "definition", style = "val") {
  miss <- template_missing_find(required, defined)
  # find sources of missing definitions
  sources <- unique(names(required[required %in% miss]))
  # show error message
  output <- c()
  len_miss <- length(miss)
  if (len_miss > 0) output <- c(
    sprintf("%s missing %s %s %s",
            format_there_were_was(len_miss, max_print = max_print),
            holder,
            as_plural(noun, len_miss, "s:"),
            character_to_string_list(miss, max_print = max_print,
                                     style = style)),
    sprintf("Required by %s.", format_sources_phrase_(sources))
  )
  output
}

#' Warn about overwritten definitions if missing definitions are found.
#' @param all_rmv `character()` Removed definitions.
#' @param holder `character(1)` Name of the definition holder.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @keywords internal
#' @noRd
template_overwrite_warn <- function(all_rmv, holder, max_print = 4) {
  len_rmv <- length(all_rmv)
  if (len_rmv > 0) warn_h(sprintf(
    "Overwriting %s %s %s %s",
    ifelse(len_rmv > max_print,
           paste0("more than ", max_print), len_rmv),
    holder,
    as_plural("definition", len_rmv, "s for:", " for"),
    character_to_string_list(all_rmv, max_print)
  ))
}

#' Warn about overwritten definitions now, or cache them on the holder.
#' @param self Definition holder with `family` and `store_warn_over`.
#' @param rmv_name `character()` Removed definitions.
#' @param warn `logical(1)` Warn now instead of caching.
#' @keywords internal
#' @noRd
warn_or_store_overwrite_ <- function(self, rmv_name, warn) {
  if (warn) template_overwrite_warn(rmv_name, self$family)
  else self$store_warn_over <- c(self$store_warn_over, rmv_name)
}

############################################################################
############################################################################
############################################################################

#' Compare required definitions to current definitions and find unused ones.
#' @param required `character()` Required definitions.
#' @param defined `character()` Current definitions.
#' @return `character()`
#' @keywords internal
#' @noRd
template_not_used_find <- function(required, defined) {
  if (length(required) == 0) sort(unlist(defined))
  else if (length(defined) == 0) NULL
  else sort(setdiff(unlist(defined), unlist(required)))
}

#' Raise a warning if unused definitions are found.
#' @param required `character()` Required definitions.
#' @param defined `character()` Current definitions.
#' @param holder `character(1)` Name of the definition holder.
#' @param max_print `integer(1)` Maximum number of elements to print.
#' @return `character()`
#' @keywords internal
#' @noRd
template_not_used_warn <- function(required, defined, holder, max_print = 4) {
  notused <- template_not_used_find(required, defined)
  len_unused <- length(notused)
  if (len_unused > 0) warn_h(
    sprintf("%s unused %s %s %s",
            format_there_were_was(len_unused, max_print),
            holder,
            as_plural("definition", len_unused, "s:"),
            character_to_string_list(notused, max_print))
  )
  invisible()
}

#' Format an assert usage definition message as a `cli::cli_bullets()` block.
#'
#' Returns a named character vector following the `cli::cli_bullets()` glyph
#' convention: the headline is named after its severity (`x` error, `!` minor,
#' `v` satisfied) and continuation lines are named `" "` to indent under it.
#' @param section `character(1)` Section name.
#' @param sms `character()` Messages to format.
#' @param minor `logical(1)` Whether to format as minor message.
#' @keywords internal
#' @noRd
fmt_assert_cli <- function(section, sms, minor = FALSE) {
  sms <- gsub(sms, pattern = " were ", replacement = " are ")
  sms <- gsub(sms, pattern = " was ", replacement = " is ")
  # colour the section name to match its bullet
  if (length(sms) == 0) return(stats::setNames(cli::col_green(section), "v"))
  col <- if (minor) cli::col_yellow else cli::col_red
  if (minor) sms <- sub(sms, pattern = "^Required by ",
                        replacement = "Can be specified for ")
  head <- paste0(col(section), ": ", tolower(substr(sms[1], 1, 1)),
                 substr(sms[1], 2, nchar(sms[1])))
  stats::setNames(c(head, sms[-1]),
                  c(if (minor) "!" else "x", rep(" ", length(sms) - 1)))
}
