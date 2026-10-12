## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' @importFrom R6 R6Class
#' @importFrom pkgbuild compile_dll
NULL

# Whether the build cache would overlap the model source.
cache_overlaps_source <- function(cache, source) {
  norm <- \(p) normalizePath(p, winslash = "/", mustWork = FALSE)
  cache <- norm(cache)
  source <- norm(source)
  cache == source || startsWith(source, paste0(cache, "/"))
}

# Assert valid household data.
# * nothing provided -> dummy_household.
assert_household_data <- function(data, model) {
  is_test_model <- is.character(model) && length(model) == 1 &&
    isTRUE(model == "test")

  if (is.null(data)) {
    # data may be supplied later with add_data_wide(): held until the wall
    if (!is_test_model) return(NULL)
    env <- new.env(parent = emptyenv())
    utils::data("dummy_household", package = "housepi", envir = env)
    return(env$dummy_household)
  }
  if (!is.data.frame(data)) {
    sms <- c("Data should be an object of class {.cls data.frame}.",
             `!` = format_type_error(data))
    if (is.character(data) && length(data) == 1) sms <- c(sms,
      "Did you mean {.code new_housepi_model(model = \"{data}\")}?"
    )
    stop_h(sms)
  }
  if (!nrow(data)) stop_h("Provided {.arg data} is empty.")
  assert_valid_data_columns(colnames(data))
  data
}

########  HousepiGenerator  ########

#' Compile and instantiate the Rcpp model object.
#' @noRd
HousepiGenerator <- R6::R6Class(
  classname = "HousepiGenerator",
  cloneable = FALSE,
  lock_class = TRUE,

  ####  Private  ####

  private = list(
    name_ = NA_character_,
    dll_ = NA_character_,
    cache_dir_ = NA_character_,
    source_dir_ = NA_character_,
    needs_compile_ = FALSE,
    mcmc_ptr_ = NULL,
    last_flags_ = character(),

    # `integer()` Infected statuses in the data.
    infec_status_set_ = integer(),

    # `character()` Columns of the last imported data_wide.
    cols_data_wide_ = character(),

    # `list` Factor levels of data_wide.
    info_data_wide_ = NULL,

    # `logical()` Whether each data_wide column is discrete, by name.
    discr_data_wide_ = NULL,

    # `character()` Columns of the last imported data_long.
    cols_data_long_ = character(),

    # `list` Factor levels of data_long.
    info_data_long_ = NULL,

    # `character()` Ordered infected-status labels, to refactor the status.
    status_levels_ = character(),

    # `character(1)` Label of status code 0, outside declared infected statuses.
    escaped_ = "escaped",

    # `character()` Infected-status labels of character or factor data,
    # in code order, the source of truth for the dictionary.
    # * numeric status -> empty.
    status_data_labels_ = character(),

    ####  Options  ####

    # Process user options for compilation.
    pre_process_user_options_ = function(quiet, ...) {
      options <- list(...)
      paths <- character()
      flags <- character()
      libs <- character()

      if ("paths" %in% names(options)) paths <- as.character(options$paths)
      if ("flags" %in% names(options)) flags <- as.character(options$flags)
      if ("libs" %in% names(options)) libs <- as.character(options$libs)

      # process time discrete option
      if (isTRUE(options$discr)) flags <- c(flags, "_HOUSEPI_DISCRETE_TIME")

      list(quiet = isTRUE(quiet), paths = paths, flags = flags, libs = libs)
    },

    # Check whether source model contains a literal macro/token.
    assert_allows_flag_ = function(source_dir, flag) {
      model_cpp <- file.path(source_dir, "model.cpp")
      lines <- readLines(model_cpp, warn = FALSE)
      patt <- paste0("\\b", gsub("([\\W])", "\\\\\\1", flag), "\\b")
      any(grepl(patt, lines, perl = TRUE))
    },

    # Assign the source directory and keep the requested flags it reads.
    # * flag absent from model.cpp -> warning.
    store_source_and_filter_flags_ = function(source_dir) {
      # check the model flags make sense
      keep <- character()
      for (f in private$last_flags_) {
        if (private$assert_allows_flag_(source_dir, f)) keep <- c(keep, f)
        else if (f == "_HOUSEPI_DISCRETE_TIME") stop_h(
          "No visible discrete time option for model {.val {self$name}}."
        )
        else stop_h(
          "No visible option associated to flag {.val {f}}."
        )
      }
      # store just effective flags
      private$last_flags_ <- keep

      # store source directory for later checks
      private$source_dir_ <- source_dir

      invisible(self)
    },

    # Whether the user flags match the existing dll.
    has_same_flags_as_dll_ = function() {
      # custom dll carries its flags in a sidecar for comparison
      if (!is.null(private$dll_) &&
            !is.na(private$dll_) &&
            file.exists(private$dll_)) {
        flags_file <- dll_flags_file(private$dll_)
        # unknown flags -> rebuild
        if (!file.exists(flags_file)) return(FALSE)
        return(setequal(private$last_flags_, readLines(flags_file)))
      }
      # built-in templates have no extra flags
      !length(private$last_flags_)
    },

    ####  Paths  ####

    # Point a custom model to its source code.
    init_custom_model_ = function(choice, as_template, force) {
      private$store_source_and_filter_flags_(choice)
      private$dll_ <- template_dll_file(choice)
      if (force ||
            !file.exists(private$dll_) ||
            untracked_changes(private$source_dir_, private$dll_) ||
            !private$has_same_flags_as_dll_())
        private$build_as_custom_model_(as_template, force)

      invisible(self)
    },

    # Point a template model to its template.
    init_template_model_ = function(choice, force) {
      private$store_source_and_filter_flags_(get_template_dir(choice))
      if (force || !private$has_same_flags_as_dll_())
        private$build_as_custom_model_(FALSE, force)

      invisible(self)
    },

    ####  Compilation  ####

    # Whether the model choice is interpretable.
    helper_init_model_choice_ = function(model) {
      choice <- "test"
      if (all(is.na(model) | is.null(model))) warn_h(
        "No model specified: falling back to {.val test} model.",
        "To show available templates, run {.run browse_templates()}."
      )
      else if (length(model) != 1 || !is.character(model)) stop_h(
        "Model choice must be a non-empty {.cls character}.",
        format_type_error(model)
      )
      else if (model %in% list(".")) stop_h(
        "Invalid choice {.val {model}} for housepi model.",
        "Move your implementation in a different folder."
      )
      else choice <- model

      private$name_ <- template_name(choice)
      choice
    },

    # Assert source files exist.
    assert_source_files_ = function() {
      file_name <- file.path(private$source_dir_, "model.cpp")
      if (!file.exists(file_name)) stop_h(
        "File {.file {file_name}} was not found.",
        "Check your model directory."
      )
      invisible(self)
    },

    # Prepare the temporary folder for the installation.
    build_temporary_cache_ = function(as_template, force) {
      private$cache_dir_ <- file.path(tempdir(), template_name(private$name_))

      # nest one level if cache overlaps source so the latter is never deleted
      if (cache_overlaps_source(private$cache_dir_, private$source_dir_))
        private$cache_dir_ <- file.path(private$cache_dir_, ".build")

      if (!dir.exists(private$cache_dir_))
        if (!dir.create(private$cache_dir_, recursive = TRUE,
                        showWarnings = FALSE)) raise_permissions_error()

      if (as_template) {
        if (file.exists(private$dll_)) {
          # check if dll is loaded
          namespace <- template_name(private$name_)
          dll_info <- base::getLoadedDLLs()[[namespace]]

          if (is.null(dll_info) || isFALSE(dll_info[["dynamicLookup"]])) {
            # safely remove previous template dll
            unlink(private$dll_)

          } else {
            private$needs_compile_ <- FALSE
            warn_h(
              "Impossible to replace existing template in current session.",
              `!` = paste("Warning raised because model changed",
                          "or {.code force = TRUE}."),
              "Start a new R session to re-compile your model as a template."
            )
          }
        }

      } else if (force) {
        private$needs_compile_ <- TRUE

      } else {
        # search cached dll and check if to be recompiled
        private$needs_compile_ <- private$search_for_cached_dll_()
      }
      invisible(self)
    },

    # Build a new custom model.
    # * `as_template` -> dynamic library stored as a local template,
    #   else recompiled in every new session;
    # * `force` -> compiled even when a cached model exists.
    build_as_custom_model_ = function(as_template, force) {
      private$needs_compile_ <- TRUE
      # temporary directory where to store source code
      # already cached -> needs_compile_ FALSE
      private$build_temporary_cache_(as_template, force)

      # check whether it still needs compiling
      if (!private$needs_compile_) return(invisible(self))

      # change dll name if not as template
      if (!as_template) private$dll_ <- file.path(
        private$cache_dir_, dll_file_ext(random_name(private$source_dir_))
      )
      # source directory must contain source files
      private$assert_source_files_()
      # organize cache directory to perform compilation
      private$copy_source_code_()

      invisible(self)
    },

    # Whether the model needs recompiling: nothing cached, or changes untracked.
    search_for_cached_dll_ = function() {
      # look for the dynamic library, the extension anchored so the flags
      # sidecar next to the dll is not picked up as a second match
      cache_dll_files <- list.files(
        private$cache_dir_,
        pattern = paste0("\\", .Platform$dynlib.ext, "$"),
        full.names = TRUE
      )
      # multiple dlls -> treated as an error
      if (length(cache_dll_files) == 1) {
        private$dll_ <- cache_dll_files[1]
        # cached dll is the one if changes are tracked and flags still match
        if (!untracked_changes(private$source_dir_, private$dll_) &&
              private$has_same_flags_as_dll_())
          return(FALSE)
      }
      # no dll cached or multiple dlls require recompilation
      TRUE
    },

    # Create temporary folder and copy source code there.
    copy_source_code_ = function() {
      # clean all previous files in the cache
      unlink(list.files(private$cache_dir_, full.names = TRUE),
             recursive = TRUE)

      # 'housepi' template package directory
      housepi_dir <- system.file("template", package = "housepi")

      # copy R folder from template
      overwrite(file.path(housepi_dir, "R"), private$cache_dir_,
                recursive = TRUE)

      # copy single package files
      for (f in file.path(housepi_dir, c("NAMESPACE", "DESCRIPTION"))) {
        overwrite(f, file.path(private$cache_dir_, basename(f)))
      }

      # copy model specific source files in the new folder
      dir.create(file.path(private$cache_dir_, "src"),
                 showWarnings = FALSE, recursive = TRUE)
      overwrite(file.path(private$source_dir_, "model.cpp"),
                file.path(private$cache_dir_, "src", "model.cpp"))

      invisible(self)
    },

    # Generate flags used for environment compilation.
    build_compile_flags_ = function(paths = NULL, flags = NULL, libs = NULL) {
      flags <- c(
        # flag to pass the model name for model factory auto-registration
        paste0("_HOUSEPI_MODULE=", get_template_from_dll(private$dll_)),
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
    },

    # Relink a custom DLL to the housepi image .engine uses, macOS only.
    # Under devtools::load_all(), the custom DLL links housepi.so by basename.
    # A reload in the same session maps two images under that basename:
    # * dyld -> binds the custom DLL to the first one;
    # * .engine -> uses the most recently loaded.
    # C++ registries are singletons per image, so the dependency is rewritten
    # to the image .engine uses.
    bind_to_active_housepi_ = function(custom_dll) {
      if (Sys.info()[["sysname"]] != "Darwin") return(invisible())
      # R's PACKAGE="housepi" lookup picks the last loaded entry, so the
      # last one is what .engine resolves against.
      hp_entries <- Filter(\(d) d[["name"]] == "housepi", getLoadedDLLs())
      housepi_path <- hp_entries[[length(hp_entries)]][["path"]]
      system2("install_name_tool",
              c("-change", "housepi.so",
                shQuote(housepi_path), shQuote(custom_dll)),
              stdout = FALSE, stderr = FALSE)
      invisible()
    },

    # Compile and load the source files inside the cache directory.
    # @param quiet `logical(1)` Whether to hide compilation commands.
    # @param user_paths `character()` Additional include paths.
    # @param user_flags `character()` Additional compiler flags.
    # @param user_libs `character()` Additional linker flags.
    compile_ = function(quiet = TRUE,
                        paths = character(),
                        flags = character(),
                        libs = character()) {
      if (!pkgbuild::has_build_tools()) stop_h(
        "No C++ compiler found to compile model {.val {self$name}}.",
        `!` = "Custom and discrete-time models compile at run time.",
        "Diagnose with {.code pkgbuild::check_build_tools(debug = TRUE)}."
      )
      # compile the dynamic library
      withr::with_envvar(
        private$build_compile_flags_(paths = paths,
                                     flags = flags,
                                     libs = libs),
        pkgbuild::compile_dll(private$cache_dir_,
                              quiet = quiet,
                              compile_attributes = FALSE,
                              debug = FALSE)
      )
      # copy DLL outside of the source directory
      custom_dll <- file.path(private$cache_dir_, "src",
                              dll_file_ext("housepi.custom"))
      private$bind_to_active_housepi_(custom_dll)
      overwrite(custom_dll, private$dll_)
      # record the flags this dll was compiled with, for later flag checks
      writeLines(private$last_flags_, dll_flags_file(private$dll_))
      # drop other files
      unlink(file.path(private$cache_dir_,
                       c("NAMESPACE", "DESCRIPTION", "src", "R")),
             recursive = TRUE)

      invisible(self)
    },

    # Create new Rcpp pointers to Model and MCMC.
    create_new_pointers_ = function() {
      if (!is.null(private$mcmc_ptr_)) return(invisible())
      load_model_dll_(private$dll_)
      module_name <- get_template_from_dll(private$dll_)
      if (is.na(module_name)) module_name <- private$name_
      private$mcmc_ptr_ <- .engine$create_mcmc_runner(module_name)
    }
  ),

  ####  Active  ####

  active = list(
    # `character(1)` Model name.
    name = function(value) {
      spaced_choice(private$name_)
    },

    # `character(1)` Model name as given, unformatted.
    raw_name = function(value) {
      private$name_
    },

    # Rcpp pointer to the MCMC runner.
    mcmc_ptr = function(value) {
      private$mcmc_ptr_
    },

    # `character(1)` Model source directory.
    source_dir = function(value) {
      private$source_dir_
    },

    # `character()` Data required by the model.
    data_required = function(value) {
      set_source_name_(.engine$key_Data(private$mcmc_ptr_), "model")
    },

    # `character()` Parameters required by the model.
    param_required = function(value) {
      set_source_name_(.engine$key_Param(private$mcmc_ptr_), "model")
    },

    # `character()` Group risks required by the model.
    group_risks_required = function(value) {
      set_source_name_(.engine$key_GroupRisk(private$mcmc_ptr_), "model")
    },

    # `character()` Latent variables required by the model.
    latent_required = function(value) {
      set_source_name_(.engine$key_Latent(private$mcmc_ptr_), "model")
    },

    # `character()` Parameters the model requires discrete.
    param_discr_required = function(value) {
      .engine$key_Param_discrete(private$mcmc_ptr_)
    },

    # `character()` Latent variables the model requires discrete.
    latent_discr_required = function(value) {
      .engine$key_Latent_discrete(private$mcmc_ptr_)
    },

    # `character()` Data columns the model requires discrete.
    data_wide_discr_required = function(value) {
      .engine$key_Data_discrete(private$mcmc_ptr_)
    },

    # `logical(1)` Whether the model allows contact risk definitions.
    has_contact_risk = function(value) {
      .engine$has_contact_risk(private$mcmc_ptr_)
    },

    # `list` Baseline parameters and the risks scaling them, keyed by parameter.
    # Each entry lists the `risks` group names,
    # and whether `contacts` scale the same term.
    risk_anchors = function(value) {
      .engine$risk_anchors(private$mcmc_ptr_)
    },

    # `character()` Group risks the model reads on infected members only.
    risk_status_stratified = function(value) {
      .engine$risk_status_stratified(private$mcmc_ptr_)
    },

    # `logical(1)` Whether the model allows infectivity profile definitions.
    has_infec_profile = function(value) {
      .engine$has_infec_profile(private$mcmc_ptr_)
    },

    # `logical(1)` Whether the model steps one whole day at once.
    is_discrete_time = function(value) {
      "_HOUSEPI_DISCRETE_TIME" %in% private$last_flags_
    },

    # `character()` Columns of the imported wide data.
    cols_data_wide = function(value) {
      private$cols_data_wide_
    },

    # `list` Levels of factors of the wide data.
    info_data_wide = function(value) {
      private$info_data_wide_
    },

    # `logical()` Whether each imported data_wide column is discrete, by name.
    discr_data_wide = function(value) {
      private$discr_data_wide_
    },

    # `integer()` Infected statuses in the wide data.
    infec_status_set = function(value) {
      private$infec_status_set_
    },

    # `character()` Infected labels seen in data.
    status_data_labels = function(value) {
      private$status_data_labels_
    },

    # `character(1)` Label for status code 0.
    escaped = function(value) {
      private$escaped_
    }
  ),

  ####  Public  ####

  public = list(
    # Create the generator, naming the dynamic library.
    # Arguments as in new_housepi_model().
    initialize = function(model = NA_character_, as_template = FALSE,
                          force = FALSE, quiet = FALSE,
                          escaped = "escaped", ...) {
      # helper construction option for copying
      if (isTRUE(.housepi_env$hollow)) return(invisible(self))

      if (length(escaped) != 1 || !is.character(escaped)) stop_h(
        "{.arg escaped} must be a non-empty {.cls character}.",
        format_type_error(escaped)
      )
      private$escaped_ <- escaped

      # assign name to the model
      choice <- private$helper_init_model_choice_(model)
      # fallback to 'test' -> no spurious template
      if (all(is.na(model) | is.null(model))) as_template <- FALSE

      # assign quiet field and load user options
      user_options <- private$pre_process_user_options_(quiet, ...)
      # store raw flags for later checks
      private$last_flags_ <- user_options$flags

      if (dir.exists(choice)) {
        private$init_custom_model_(choice, as_template, force)
      } else if (!is.na(get_template_dir(choice))) {
        if (isTRUE(as_template)) warn_h(
          "Cannot store template {.val {self$name}} as a new template.",
          "Copy it to a writable folder first, and then load it:",
          ` ` = "#> {.code browse_templates(\"{self$name}\", \"your/dir\")}",
          ` ` = "#> {.code new_housepi_model(data = ..., model = \"your/dir\")}"
        )
        private$init_template_model_(choice, force)
      } else {
        stop_h("{.val {choice}} is not a folder nor a template name.",
               find_closest_template_name(choice))
      }

      # if necessary, recompile only the model script
      if (private$needs_compile_) {
        # drop previous object
        private$mcmc_ptr_  <- NULL
        do.call(private$compile_, user_options)
        private$needs_compile_ <- FALSE
      }
      # generate the Rcpp model object if missing
      private$create_new_pointers_()

      invisible(self)
    },

    # Whether the model is the test model.
    is_test_model = function() {
      private$name_ == "test"
    },

    # Record the ordered infected-status labels for refactoring.
    # @param levels `character()` Status labels in code order (1..K).
    set_status_levels = function(levels) {
      private$status_levels_ <- as.character(levels)
      invisible(self)
    },

    # Refactor a wide-format data.frame back to user levels.
    # @todo drop the refactoring of variant data
    refactor_data_wide = function(x) {
      cpp_data_refactor(x, private$info_data_wide_, private$cols_data_wide_,
                        private$status_levels_, private$escaped_)
    },

    # Refactor a long-format data.frame back to user levels.
    # @todo drop the refactoring of variant data
    refactor_data_long = function(x) {
      cpp_data_refactor(x, private$info_data_long_, private$cols_data_long_,
                        private$status_levels_, private$escaped_)
    },

    # Allocate fresh Rcpp pointers (used after a copy).
    init_pointers = function() {
      private$create_new_pointers_()
      invisible(self)
    },

    # Validate and factorize wide-format data, populating `cols_data_wide`
    # for the assertions before the push to the engine.
    # @return `data.frame` Factorized data, ready for push_data_wide().
    factorize_data_wide = function(data) {
      embed <- cpp_data_factorize(data, private$escaped_)
      private$info_data_wide_ <- embed$info
      private$cols_data_wide_ <- embed$cols
      private$discr_data_wide_ <- embed$discr
      private$infec_status_set_ <- embed$status
      private$status_data_labels_ <- embed$status_labels
      embed$data
    },

    # Validate, factorize, and push long-format data into the engine.
    import_data_long = function(data_long) {
      if (!is.data.frame(data_long)) stop_h(
        "Data-long should be an object of class {.cls data.frame}.",
        format_type_error(data_long)
      )
      assert_valid_data_columns(colnames(data_long), "id")

      embed <- cpp_data_factorize(data_long, private$escaped_)
      private$info_data_long_ <- embed$info
      private$cols_data_long_ <- embed$cols
      .engine$import_data_long(private$mcmc_ptr_, embed$data)
      invisible(self)
    }
  )
)
