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

# Templated name convention (non-identifier characters become underscores).
template_name <- function(path) {
  gsub("[^A-Za-z0-9_]", "_", basename(path))
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

# Path of the housepi image .engine uses:
# R's PACKAGE="housepi" lookup picks the last loaded entry.
active_housepi_lib_ <- function() {
  hp_entries <- Filter(\(d) d[["name"]] == "housepi", getLoadedDLLs())
  hp_entries[[length(hp_entries)]][["path"]]
}

# Module name of a model build: `name` and a hash of what its dll is built from,
# the source, the compile options, housepi headers and image.
module_name_ <- function(name, source_dir, options) {
  headers <- list.files(system.file("include", package = "housepi"),
                        recursive = TRUE, full.names = TRUE)
  files <- c(file.path(source_dir, "model.cpp"), headers, active_housepi_lib_())
  key <- rlang::hash(list(rlang::hash_file(files), sort(options$flags),
                          options$paths, options$libs))
  paste0(name, "_", substr(key, 1, 12))
}

# Whether model `name` is the test model.
is_test_model_ <- function(name) {
  identical(name, "test")
}

# Assert valid household data.
# * nothing provided -> dummy_household.
assert_household_data <- function(data, model) {
  if (is.null(data)) {
    # data may be supplied later with add_dataset(): held until the wall
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

# Assert source files exist.
assert_source_files_ <- function(source_dir) {
  file_name <- file.path(source_dir, "model.cpp")
  if (!file.exists(file_name)) stop_h(
    "File {.file {file_name}} was not found.",
    "Check your model directory."
  )
  invisible()
}

# Create the build folder and copy the source code there.
copy_source_code_ <- function(build_dir, source_dir) {
  # leftovers of an interrupted build
  unlink(build_dir, recursive = TRUE)
  dir.create(file.path(build_dir, "src"),
             showWarnings = FALSE, recursive = TRUE)

  # 'housepi' template package directory
  housepi_dir <- system.file("template", package = "housepi")

  # copy R folder from template
  overwrite(file.path(housepi_dir, "R"), build_dir, recursive = TRUE)

  # copy single package files
  for (f in file.path(housepi_dir, c("NAMESPACE", "DESCRIPTION"))) {
    overwrite(f, file.path(build_dir, basename(f)))
  }

  # copy model specific source files in the new folder
  overwrite(file.path(source_dir, "model.cpp"),
            file.path(build_dir, "src", "model.cpp"))
  invisible()
}

# Generate flags used for environment compilation.
compile_envvars_ <- function(dll, paths = NULL, flags = NULL, libs = NULL) {
  flags <- c(
    # flag to pass the model name for model factory auto-registration
    paste0("_HOUSEPI_MODULE=", get_template_from_dll(dll)),
    flags
  )
  paths <- c(system.file("include", package = "housepi"), paths)
  housepi_lib <- active_housepi_lib_()

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
  out <- suppressWarnings(system2(
    "install_name_tool",
    c("-change", "housepi.so",
      shQuote(active_housepi_lib_()), shQuote(custom_dll)),
    stdout = TRUE, stderr = TRUE
  ))
  if (!is.null(attr(out, "status"))) stop_h(
    "Failed to link the model to the loaded {.pkg housepi} library.",
    `!` = "{.code install_name_tool} reported: {out}"
  )
  invisible()
}

# Compile the source files in a build folder into the dll of `compiled`.
# @param compiled `list` As compile_model_() returns it.
# @param options `list` As compile_options_() returns it.
compile_dll_ <- function(compiled, options) {
  # nolint next: object_usage_linter. read by the cli glue below
  label <- spaced_choice(compiled$name)
  if (!pkgbuild::has_build_tools()) stop_h(
    "No C++ compiler found to compile model {.val {label}}.",
    `!` = "Custom and discrete-time models compile at run time.",
    "Diagnose with {.run pkgbuild::check_build_tools(debug = TRUE)}."
  )
  # build folder in the session cache, dropped once the dll is copied out
  build_dir <- file.path(tempdir(), "housepi_models",
                         get_template_from_dll(compiled$dll))
  on.exit(unlink(build_dir, recursive = TRUE), add = TRUE)
  copy_source_code_(build_dir, compiled$source_dir)
  withr::with_envvar(
    compile_envvars_(compiled$dll, paths = options$paths,
                     flags = options$flags, libs = options$libs),
    pkgbuild::compile_dll(build_dir,
                          quiet = options$quiet,
                          compile_attributes = FALSE,
                          debug = FALSE)
  )
  custom_dll <- file.path(build_dir, "src", dll_file_ext("housepi.custom"))
  bind_to_active_housepi_(custom_dll)
  overwrite(custom_dll, compiled$dll)
  invisible()
}

# Compile the transmission model, or reuse its dll.
# Arguments as in new_housepi_model().
# @return `list` The compiled model:
# * `name` -> as chosen, unformatted;
# * `source_dir` -> folder of its `model.cpp`;
# * `dll` -> dynamic library, `NA` for a template the package ships built;
# * `flags` -> compiler flags;
# * `module` -> engine identity, `name` and the hash of its build.
compile_model_ <- function(model = NA_character_, as_template = FALSE,
                           force = FALSE, quiet = FALSE, ...) {
  choice <- model_choice_(model)
  name <- template_name(choice)
  label <- spaced_choice(name)
  # fallback to 'test' -> no spurious template
  if (all(is.na(model) | is.null(model))) as_template <- FALSE
  options <- compile_options_(quiet, ...)

  custom <- dir.exists(choice)
  if (custom) {
    source_dir <- choice
  } else if (!is.na(get_template_dir(choice))) {
    if (isTRUE(as_template)) warn_h(
      "Cannot store template {.val {label}} as a new template.",
      "Copy it to a writable folder first, and then load it:",
      ` ` = "#> {.code browse_templates(\"{label}\", \"your/dir\")}",
      ` ` = "#> {.code new_housepi_model(data = ..., model = \"your/dir\")}"
    )
    as_template <- FALSE
    source_dir <- get_template_dir(choice)
  } else {
    stop_h("{.val {choice}} is not a folder nor a template name.",
           find_closest_template_name(choice))
  }
  assert_source_files_(source_dir)
  assert_model_flags_(source_dir, options$flags, label)

  compiled <- list(name = name, source_dir = source_dir, dll = NA_character_,
                   flags = options$flags, module = name)
  # built-in template without flags -> the dll the package ships
  if (!custom && !length(options$flags) && !force) return(compiled)

  # dll named after its build: an edited source or option gets a new one
  # * a template in the model folder, then the session cache;
  # * `force` -> rebuilt in the last, under the next suffix free in session.
  module <- compiled$module <- module_name_(name, source_dir, options)
  dirs <- c(if (custom) source_dir,
            if (!as_template) file.path(tempdir(), "housepi_models"))
  dlls <- file.path(dirs, dll_file_ext(module))
  built <- dlls[file.exists(dlls)]
  compiled$dll <- if (force || !length(built)) dlls[length(dlls)] else built[1]
  if (!force && length(built)) return(compiled)
  n <- 0
  while (is_dll_loaded(compiled$dll)) {
    n <- n + 1
    compiled$dll <- file.path(dirs[length(dirs)],
                              dll_file_ext(paste0(module, "_", n)))
  }
  compile_dll_(compiled, options)

  # earlier builds of a template, unless loaded in session
  if (as_template) {
    pattern <- paste0("^", name, "_[0-9a-f]{12}(_[0-9]+)?\\",
                      .Platform$dynlib.ext, "$")
    old <- setdiff(list.files(source_dir, pattern, full.names = TRUE),
                   compiled$dll)
    unlink(Filter(Negate(is_dll_loaded), old))
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
