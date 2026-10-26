## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' @importFrom R6 R6Class
#' @importFrom pkgbuild compile_dll
NULL

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

# Whether the build cache would overlap the model source.
cache_overlaps_source <- function(cache, source) {
  norm <- \(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  cache <- norm(cache)
  source <- norm(source)
  cache == source || startsWith(source, paste0(cache, "/"))
}

# Whether model `name` is the test model.
is_test_model_ <- function(name) {
  identical(name, "test")
}

# Assert valid household data.
# * nothing provided -> dummy_household.
assert_household_data <- function(data, model) {
  if (is.null(data)) {
    # data may be supplied later with add_data_wide(): held until the wall
    if (!is_test_model_(model)) return(NULL)
    env <- new.env(parent = emptyenv())
    utils::data("dummy_household", package = "housepi", envir = env)
    return(env$dummy_household)
  }
  if (!is.data.frame(data)) {
    sms <- c("Data should be an object of class {.cls data.frame}.",
             `!` = format_type_error(data))
    if (is.character(data) && length(data) == 1) sms <- c(sms,
      paste0("Did you mean {.run [new_housepi_model(model = \"{data}\")]",
             "(housepi::new_housepi_model(model = \"{data}\"))}?")
    )
    stop_h(sms)
  }
  if (!nrow(data)) stop_h("Provided {.arg data} is empty.")
  assert_valid_data_columns(colnames(data))
  data
}

########  Compilation  ########

# Interpret the model choice.
# * nothing given -> `test`, with a warning.
model_choice_ <- function(model) {
  if (all(is.na(model) | is.null(model))) {
    warn_h(
      "No model specified: falling back to {.val test} model.",
      browse_templates_hint_
    )
    return("test")
  }
  if (length(model) != 1 || !is.character(model)) stop_h(
    "Model choice must be a non-empty {.cls character}.",
    format_type_error(model)
  )
  if (model %in% list(".")) stop_h(
    "Invalid choice {.val {model}} for housepi model.",
    "Move your implementation in a different folder."
  )
  model
}

# Process user options for compilation.
compile_options_ <- function(quiet, ...) {
  assert_advanced_settings_(c("discr", "flags", "libs", "paths"), ...)
  options <- list(...)
  paths <- as.character(options[["paths"]])
  flags <- as.character(options[["flags"]])
  libs <- as.character(options[["libs"]])

  # process time discrete option
  if (isTRUE(options[["discr"]])) flags <- c(flags, "_HOUSEPI_DISCRETE_TIME")

  list(quiet = isTRUE(quiet), paths = paths, flags = flags, libs = libs)
}

# Check whether source model contains a literal macro/token.
model_reads_flag_ <- function(source_dir, flag) {
  model_cpp <- file.path(source_dir, "model.cpp")
  lines <- readLines(model_cpp, warn = FALSE)
  patt <- paste0("\\b", gsub("([\\W])", "\\\\\\1", flag), "\\b")
  any(grepl(patt, lines, perl = TRUE))
}

# Assert the model source reads every requested flag.
# @param label `character(1)` Model name, as printed.
assert_model_flags_ <- function(source_dir, flags, label) {
  for (f in flags) {
    if (model_reads_flag_(source_dir, f)) next
    if (f == "_HOUSEPI_DISCRETE_TIME") stop_h(
      "No visible discrete time option for model {.val {label}}."
    )
    stop_h("No visible option associated to flag {.val {f}}.")
  }
  invisible()
}

# Whether the flags match the existing dll.
same_flags_as_dll_ <- function(dll, flags) {
  # custom dll carries its flags in a sidecar for comparison
  if (!is.na(dll) && file.exists(dll)) {
    flags_file <- dll_flags_file(dll)
    # unknown flags -> rebuild
    if (!file.exists(flags_file)) return(FALSE)
    return(setequal(flags, readLines(flags_file)))
  }
  # built-in templates have no extra flags
  !length(flags)
}

# Assert source files exist.
assert_source_files_ <- function(source_dir) {
  file_name <- file.path(source_dir, "model.cpp")
  if (!file.exists(file_name)) stop_h(
    "File {.file {file_name}} was not found.",
    "Check your model directory."
  )
  invisible()
}

# Cached dll still matching the source and the flags.
# * nothing cached, or several dlls -> `NA`.
find_cached_dll_ <- function(cache_dir, source_dir, flags) {
  # the extension anchored so the flags sidecar next to the dll
  # is not picked up as a second match
  dlls <- list.files(cache_dir,
                     pattern = paste0("\\", .Platform$dynlib.ext, "$"),
                     full.names = TRUE)
  if (length(dlls) == 1 &&
        !untracked_changes(source_dir, dlls) &&
        same_flags_as_dll_(dlls, flags))
    return(dlls)
  NA_character_
}

# Create the cache folder and copy the source code there.
copy_source_code_ <- function(cache_dir, source_dir) {
  # clean all previous files in the cache
  unlink(list.files(cache_dir, full.names = TRUE), recursive = TRUE)

  # 'housepi' template package directory
  housepi_dir <- system.file("template", package = "housepi")

  # copy R folder from template
  overwrite(file.path(housepi_dir, "R"), cache_dir, recursive = TRUE)

  # copy single package files
  for (f in file.path(housepi_dir, c("NAMESPACE", "DESCRIPTION"))) {
    overwrite(f, file.path(cache_dir, basename(f)))
  }

  # copy model specific source files in the new folder
  dir.create(file.path(cache_dir, "src"),
             showWarnings = FALSE, recursive = TRUE)
  overwrite(file.path(source_dir, "model.cpp"),
            file.path(cache_dir, "src", "model.cpp"))
  invisible()
}

# Prepare the cache folder of a model to rebuild.
# * `as_template` -> dll stored as a local template, replaced unless loaded;
# * `force` -> compiled even when a cached model exists;
# * otherwise -> a cached dll matching source and flags is reused.
# @param compiled `list` As compile_model_() returns it.
# @return `list` The `dll`, its `cache_dir`, and whether to `compile` it.
stage_build_ <- function(compiled, as_template, force) {
  name <- compiled$name
  source_dir <- compiled$source_dir
  dll <- compiled$dll
  cache_dir <- file.path(tempdir(), template_name(name))

  # nest one level if cache overlaps source so the latter is never deleted
  if (cache_overlaps_source(cache_dir, source_dir))
    cache_dir <- file.path(cache_dir, ".build")

  if (!dir.exists(cache_dir))
    if (!dir.create(cache_dir, recursive = TRUE,
                    showWarnings = FALSE)) raise_permissions_error()

  compile <- TRUE
  if (as_template) {
    if (file.exists(dll)) {
      # check if dll is loaded
      dll_info <- base::getLoadedDLLs()[[template_name(name)]]

      if (is.null(dll_info) || isFALSE(dll_info[["dynamicLookup"]])) {
        # safely remove previous template dll
        unlink(dll)

      } else {
        compile <- FALSE
        warn_h(
          "Impossible to replace existing template in current session.",
          `!` = paste("Warning raised because model changed",
                      "or {.code force = TRUE}."),
          "Start a new R session to re-compile your model as a template."
        )
      }
    }

  } else if (!force) {
    # cached dll matching source and flags -> no compilation
    cached <- find_cached_dll_(cache_dir, source_dir, compiled$flags)
    if (!is.na(cached)) {
      dll <- cached
      compile <- FALSE
    }
  }

  if (compile) {
    # change dll name if not as template
    if (!as_template)
      dll <- file.path(cache_dir, dll_file_ext(random_name(source_dir)))
    assert_source_files_(source_dir)
    copy_source_code_(cache_dir, source_dir)
  }
  list(dll = dll, cache_dir = cache_dir, compile = compile)
}

# Generate flags used for environment compilation.
compile_envvars_ <- function(dll, paths = NULL, flags = NULL, libs = NULL) {
  flags <- c(
    # flag to pass the model name for model factory auto-registration
    paste0("_HOUSEPI_MODULE=", get_template_from_dll(dll)),
    flags
  )
  paths <- c(system.file("include", package = "housepi"), paths)
  housepi_lib <- getLoadedDLLs()[["housepi"]][["path"]]

  # collect RcppParallel flags
  rcpp_parallel_libs <- paste(utils::capture.output(
    RcppParallel::RcppParallelLibs()
  ), collapse = " ")

  # quote each token: make hands the recipe to the shell
  c(
    PKG_CPPFLAGS = paste(
      shQuote(c(paste0("-D", flags), paste0("-I", paths))), collapse = " "
    ),
    PKG_LIBS = paste(
      c(rcpp_parallel_libs,
        shQuote(c(if (length(libs)) paste0("-L", libs), housepi_lib))
      ), collapse = " "
    )
  )
}

# Relink a custom DLL to the housepi image .engine uses, macOS only.
# Under devtools::load_all(), the custom DLL links housepi.so by basename.
# A reload in the same session maps two images under that basename:
# * dyld -> binds the custom DLL to the first one;
# * .engine -> uses the most recently loaded.
# C++ registries are singletons per image, so the dependency is rewritten
# to the image .engine uses.
bind_to_active_housepi_ <- function(custom_dll) {
  if (Sys.info()[["sysname"]] != "Darwin") return(invisible()) # nocov
  # R's PACKAGE="housepi" lookup picks the last loaded entry, so the
  # last one is what .engine resolves against.
  hp_entries <- Filter(\(d) d[["name"]] == "housepi", getLoadedDLLs())
  housepi_path <- hp_entries[[length(hp_entries)]][["path"]]
  system2("install_name_tool",
          c("-change", "housepi.so",
            shQuote(housepi_path), shQuote(custom_dll)),
          stdout = FALSE, stderr = FALSE)
  invisible()
}

# Compile and load the source files inside the cache directory.
# @param compiled `list` As compile_model_() returns it.
# @param options `list` As compile_options_() returns it.
compile_dll_ <- function(cache_dir, compiled, options) {
  # nolint next: object_usage_linter. read by the cli glue below
  label <- spaced_choice(compiled$name)
  if (!pkgbuild::has_build_tools()) stop_h(
    "No C++ compiler found to compile model {.val {label}}.",
    `!` = "Custom and discrete-time models compile at run time.",
    "Diagnose with {.run pkgbuild::check_build_tools(debug = TRUE)}."
  )
  # compile the dynamic library
  withr::with_envvar(
    compile_envvars_(compiled$dll, paths = options$paths,
                     flags = options$flags, libs = options$libs),
    pkgbuild::compile_dll(cache_dir,
                          quiet = options$quiet,
                          compile_attributes = FALSE,
                          debug = FALSE)
  )
  # copy DLL outside of the source directory
  custom_dll <- file.path(cache_dir, "src", dll_file_ext("housepi.custom"))
  bind_to_active_housepi_(custom_dll)
  overwrite(custom_dll, compiled$dll)
  # record the flags this dll was compiled with, for later flag checks
  writeLines(compiled$flags, dll_flags_file(compiled$dll))
  # drop other files
  unlink(file.path(cache_dir, c("NAMESPACE", "DESCRIPTION", "src", "R")),
         recursive = TRUE)
  invisible()
}

# Compile the transmission model, or reuse its dll.
# Arguments as in new_housepi_model().
# @return `list` The compiled model:
# * `name` -> as chosen, unformatted;
# * `source_dir` -> folder of its `model.cpp`;
# * `dll` -> dynamic library, `NA` for a template the package ships built;
# * `flags` -> compiler flags.
compile_model_ <- function(model = NA_character_, as_template = FALSE,
                           force = FALSE, quiet = FALSE, ...) {
  choice <- model_choice_(model)
  name <- template_name(choice)
  label <- spaced_choice(name)
  # fallback to 'test' -> no spurious template
  if (all(is.na(model) | is.null(model))) as_template <- FALSE
  options <- compile_options_(quiet, ...)

  if (dir.exists(choice)) {
    source_dir <- choice
    dll <- template_dll_file(choice)
  } else if (!is.na(get_template_dir(choice))) {
    if (isTRUE(as_template)) warn_h(
      "Cannot store template {.val {label}} as a new template.",
      "Copy it to a writable folder first, and then load it:",
      ` ` = "#> {.code browse_templates(\"{label}\", \"your/dir\")}",
      ` ` = "#> {.code new_housepi_model(data = ..., model = \"your/dir\")}"
    )
    as_template <- FALSE
    source_dir <- get_template_dir(choice)
    dll <- NA_character_
  } else {
    stop_h("{.val {choice}} is not a folder nor a template name.",
           find_closest_template_name(choice))
  }
  assert_model_flags_(source_dir, options$flags, label)

  compiled <- list(name = name, source_dir = source_dir, dll = dll,
                   flags = options$flags)
  # custom model -> rebuilt when its dll is missing or outdated
  stale <- !is.na(dll) &&
    (!file.exists(dll) || untracked_changes(source_dir, dll))
  if (force || stale || !same_flags_as_dll_(dll, options$flags)) {
    build <- stage_build_(compiled, as_template, force)
    compiled$dll <- build$dll
    if (build$compile) compile_dll_(build$cache_dir, compiled, options)
  }
  compiled
}

# Create the MCMC runner of a compiled model, loading its dll.
# @param compiled `list` As compile_model_() returns it.
create_mcmc_runner_ <- function(compiled) {
  load_model_dll_(compiled$dll)
  module_name <- get_template_from_dll(compiled$dll)
  if (is.na(module_name)) module_name <- compiled$name
  .engine$create_mcmc_runner(module_name)
}
