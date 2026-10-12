## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' @importFrom Rcpp Module

# Name with the dynamic library extension.
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

# Append a random string to a templated name.
random_name <- function(path) {
  # generate a random string
  str_len <- 8
  str_name <- paste0(sample(c(letters, LETTERS, 0:9), str_len, replace = TRUE),
                     collapse = "")
  paste0(template_name(path), str_name)
}

# DLL file name, by the templated name convention.
template_dll_file <- function(path) {
  file.path(path, dll_file_ext(template_name(path)))
}

# Directory of all housepi template models.
template_models_dir <- function() {
  system.file("models", package = "housepi")
}

# Template name of a DLL file.
get_template_from_dll <- function(dll) {
  if (is.null(dll) || is.na(dll)) NA_character_
  else tools::file_path_sans_ext(basename(dll))
}

# Whether the DLL is already loaded.
is_dll_loaded <- function(dll) {
  any(vapply(getLoadedDLLs(), \(x) x[["path"]] == dll, FALSE))
}

# Load the model DLL unless already loaded.
# Its static init registers the model in the model registry.
load_model_dll_ <- function(dll) {
  if (!is.null(dll) && !is.na(dll) && !is_dll_loaded(dll)) dyn.load(dll)
}

#' Directory of template model `x`.
#' @param x `character(1)` Template name.
#' @noRd
get_template_dir <- function(x) {
  flat_choice <- spaced_choice(x)
  path_to_options <- list.files(template_models_dir(), full.names = TRUE)
  flat_options <- spaced_choice(basename(path_to_options))
  # extract path if present (unique since folder name)
  mask <- flat_options == flat_choice
  if (any(mask)) return(path_to_options[mask])
  NA_character_
}

#' Suggest the closest template name to a mistyped `x`, as "did you mean".
#' @noRd
find_closest_template_name <- function(x) {
  flat_choice <- spaced_choice(x)
  options <- list.files(template_models_dir(), full.names = FALSE)
  flat_options <- spaced_choice(options)
  # no suggestion for trivial mistakes
  sms <- find_closest_match(flat_choice, flat_options)
  if (endsWith(sms, "?") && nchar(flat_choice) > 3) return(sms)
  c("Show available templates calling {.fn browse_templates}.")
}

#' Show the names of the available template models.
#' @noRd
show_model_names <- function() {
  model_names <- list.files(template_models_dir())
  links <- vapply(
    tolower(gsub("-| |_", "_", model_names)),
    \(m) inline_cli_("{.run [{m}](housepi::browse_templates(\"{m}\"))}"),
    character(1), USE.NAMES = FALSE
  )
  cli::cli_inform(c(i = "Available template models:",
                    paste(cli::symbol$bullet, links)))
}

#' Copy template model `name` to directory `path`.
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

  # source cpp is always present
  required <- "model.cpp"
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

#' Open the help page of template model `name`.
#' @noRd
open_model_helper_ <- function(name) {
  print(utils::help(paste0("template_", name), package = "housepi"))
  invisible()
}

#' @title Navigate through household model templates
#' @description Print the list of available templates,
#' get their descriptions or copy their source code to a custom path.
#' @param name `character(1)` Template name.
#' If `NA`, print the list of available templates.
#' @param path `character(1)` Path to copy the template scripts to.
#' If `NA` and `name` is not, open the template help page.
#' @param replace `logical(1)` Whether to replace the content of `path`.
#' @return `NULL`, invisibly.
#' The function prints the templates, opens a description or copies the files.
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
  else if (is.na(path)) open_model_helper_(basename(get_template_dir(name)))
  # asking for copying template to folder with model
  else if (file.exists(file.path(path, "model.cpp")) && !replace) warn_h(
    "A model already exists in the specified folder.",
    "To replace it, run again with option {.code replace = TRUE}."
  )
  # asking for copying template to empty folder
  else copy_template_to_path(name, path)
  invisible()
}

########  Missing definitions  ########

# Name every element of a vector alike.
set_source_name_ <- function(x, name = NA_character_) {
  if (isTRUE(nzchar(name, keepNA = TRUE)) && length(x))
    names(x) <- rep(name, length(x))
  x
}

#' Required definitions missing from `defined`.
#' @noRd
template_missing_find <- function(required, defined) {
  if (!length(required)) NULL
  else if (!length(defined)) sort(unlist(required))
  else sort(setdiff(unlist(required), unlist(defined)))
}

#' Messages on the required definitions missing from `defined`.
#' @param holder `character(1)` Name of the definition holder.
#' @param noun `character(1)` Word naming the missing entries.
#' @param style `character(1)` cli markup class wrapping each entry.
#' @return `character()` | `NULL` Messages, `NULL` when none is missing.
#' @noRd
template_missing_stop <- function(required, defined, holder,
                                  noun = "definition", style = "val") {
  miss <- template_missing_find(required, defined)
  if (!length(miss)) return(NULL)
  # find sources of missing definitions
  sources <- unique(names(required[required %in% miss]))
  sources <- ifelse(sources == "model", "transmission model",
                    paste(sources, "definitions"))
  c(inline_cli_("{cli::qty(miss)}There {?was/were} {length(miss)} missing ",
                "{holder} {noun}{cli::qty(miss)}{?s}: {.", style,
                " {miss}}."),
    inline_cli_("Required by {sources}."))
}

#' Warn on the definitions `all_rmv` overwrote.
#' @param all_rmv `character()` Removed definitions.
#' @param holder `character(1)` Name of the definition holder.
#' @noRd
template_overwrite_warn <- function(all_rmv, holder) {
  len_rmv <- length(all_rmv)
  if (len_rmv > 0) warn_h(paste(
    "Overwriting {len_rmv} {holder} definition{cli::qty(len_rmv)}{?s}",
    "for: {.val {all_rmv}}."
  ))
}

#' Warn about overwritten definitions now, or cache them on the holder.
#' @param self Definition holder with `family` and `store_warn_over`.
#' @param rmv_name `character()` Removed definitions.
#' @param warn `logical(1)` Warn now instead of caching.
#' @noRd
warn_or_store_overwrite_ <- function(self, rmv_name, warn) {
  if (warn) template_overwrite_warn(rmv_name, self$family)
  else self$store_warn_over <- c(self$store_warn_over, rmv_name)
}

########  Unused definitions  ########

#' Definitions in `defined` that are not required.
#' @noRd
template_not_used_find <- function(required, defined) {
  if (!length(required)) sort(unlist(defined))
  else if (!length(defined)) NULL
  else sort(setdiff(unlist(defined), unlist(required)))
}

#' Warn on the definitions in `defined` that are not required.
#' @param holder `character(1)` Name of the definition holder.
#' @noRd
template_not_used_warn <- function(required, defined, holder) {
  notused <- template_not_used_find(required, defined)
  len_unused <- length(notused)
  if (len_unused > 0) warn_h(paste(
    "{cli::qty(len_unused)}There {?was/were} {len_unused} unused {holder}",
    "definition{cli::qty(len_unused)}{?s}: {.val {notused}}."
  ))
  invisible()
}

#' Format an assertion message as a `cli::cli_bullets()` block.
#' * headline -> named after its severity: `x` error, `!` minor, `v` satisfied;
#' * continuation line -> named `" "`, indented under it.
#' @param section `character(1)` Section name.
#' @param sms `character()` Messages to format.
#' @param minor `logical(1)` Whether to format as minor message.
#' @noRd
fmt_assert_cli <- function(section, sms, minor = FALSE) {
  sms <- gsub(sms, pattern = " were ", replacement = " are ")
  sms <- gsub(sms, pattern = " was ", replacement = " is ")
  # colour the section name to match its bullet
  if (!length(sms)) return(stats::setNames(cli::col_green(section), "v"))
  col <- if (minor) cli::col_yellow else cli::col_red
  if (minor) sms <- sub(sms, pattern = "^Required by ",
                        replacement = "Can be specified for ")
  head <- paste0(col(section), ": ", tolower(substr(sms[1], 1, 1)),
                 substr(sms[1], 2, nchar(sms[1])))
  stats::setNames(c(head, sms[-1]),
                  c(if (minor) "!" else "x", rep(" ", length(sms) - 1)))
}
