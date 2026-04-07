## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom R6 R6Class
#' @importFrom pkgbuild compile_dll
NULL

############################################################################
############################################################################
############################# Housepi Generator ############################
############################################################################
############################################################################

#' @title Housepi Generator
#' @description
#' Help with generating Rcpp model object.
#' @keywords internal
#' @noRd
HousepiGenerator <- R6::R6Class(
  classname = "HousepiGenerator",
  cloneable = FALSE,
  lock_class = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    name_ = NA_character_,
    dll_ = NA_character_,
    cache_dir_ = NA_character_,
    source_dir_ = NA_character_,
    needs_compile_ = FALSE,
    mcmc_ptr_ = NULL,
    last_flags_ = character(),

    # @field infec_status_set_ `integer()`\cr
    # All positive status extracted from data.
    infec_status_set_ = integer(),

    # @field cols_data_wide_ `character()`\cr
    # Headers of last imported data wide.
    cols_data_wide_ = character(),

    # @field info_data_wide_ `list`\cr
    # List of levels of the factors in data wide.
    info_data_wide_ = NULL,

    # @field cols_data_long_ `character()`\cr
    # Headers of the last data long imported.
    cols_data_long_ = character(),

    # @field info_data_long_ `list`\cr
    # List of levels of the factors in the data long.
    info_data_long_ = NULL,

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
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

    # @description
    # Check whether source model contains a literal macro/token.
    assert_allows_flag_ = function(source_dir, flag) {
      model_cpp <- file.path(source_dir, "model.cpp")
      lines <- readLines(model_cpp, warn = FALSE)
      patt <- paste0("\\b", gsub("([\\W])", "\\\\\\1", flag), "\\b")
      any(grepl(patt, lines, perl = TRUE))
    },

    # @description
    # Assign source directory and filter requested flags according to source.
    # If a flag is not present in model.cpp, raise a warning.
    store_source_and_filter_flags_ = function(source_dir) {
      # asses if model flags actually make sense
      keep <- character()
      for (f in private$last_flags_) {
        if (private$assert_allows_flag_(source_dir, f)) keep <- c(keep, f)
        else if (f == "_HOUSEPI_DISCRETE_TIME") stop_h(
          sprintf("No visible discrete time option for model '%s'.", self$name)
        )
        else stop_h(
          sprintf("No visible option associated to flag '%s'.", f),
          "Only flags used in model.cpp are applicable at compile time."
        )
      }
      # store just effective flags
      private$last_flags_ <- keep

      # store source directory for later checks
      private$source_dir_ <- source_dir

      invisible(self)
    },

    # @description
    # Wheter user flags are the same as existing dll.
    has_same_flags_as_dll_ = function() {
      # FIXME: this is simply wrong!
      if (is.null(private$mcmc_ptr_)) {
        if (length(private$last_flags_) == 0) return(TRUE)
        else return(FALSE)
      }
      load_model_dll_(private$dll_)
      module_name <- get_template_from_dll(private$dll_)
      if (is.na(module_name)) module_name <- private$name_
      flags_dll <- .engine$compile_flags(.engine$new_model(module_name))
      setequal(private$last_flags_, flags_dll)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Initialize custom model loading path to model source code.
    init_custom_model_ = function(choice, as_template, force) {
      private$store_source_and_filter_flags_(choice)
      private$dll_ <- template_dll_file(choice)
      if (force ||
            !file.exists(private$dll_) ||
            untracked_changes(private$source_dir_, private$dll_) ||
            !private$has_same_flags_as_dll_()) {
        private$build_as_custom_model_(as_template, force)
      }
      invisible(self)
    },

    # @description
    # Initialize template model loading template model name.
    init_template_model_ = function(choice) {
      private$store_source_and_filter_flags_(get_template_dir(choice))
      if (!private$has_same_flags_as_dll_()) {
        private$build_as_custom_model_(FALSE, FALSE)
      }
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Wheter model choice is interpretable.
    helper_init_choice_ = function(choice) {
      choice_ <- "blank"
      if (is.na(choice)) warning_h(
        "No model specified. Falling back to 'blank' model.",
        "To show available templates, run 'template_household_model'."
      )
      else if (as.character(choice) %in% list(".")) stop_h(
        sprintf("Invalid choice '%s' for housepi model.", choice),
        "Move your implementation in a different folder."
      )
      else choice_ <- as.character(choice)
      private$name_ <- template_name(choice_)
      choice_
    },

    # @description
    # Assert source files exist.
    assert_source_files_ = function() {
      file_name <- file.path(private$source_dir_, "model.cpp")
      if (!file.exists(file_name)) stop_h(
        sprintf("File '%s' was not found.", file_name),
        "Check your model directory."
      )
      invisible(self)
    },

    # @description
    # Predispose temporare fodler for the installation.
    build_temporary_cache_ = function(as_template, force) {
      private$cache_dir_ <- file.path(tempdir(), template_name(private$name_))

      if (!dir.exists(private$cache_dir_)) {
        if (!dir.create(private$cache_dir_, recursive = TRUE,
                        showWarnings = FALSE)) stop_h(
          "Impossibile to create a cache directory for compilation.",
          "Check your user permissions."
        )
      }

      if (as_template) {
        if (file.exists(private$dll_)) {
          # check if dll is loaded
          namespace <- template_name(private$name_)
          dll_info <- base::getLoadedDLLs()[[namespace]]

          if (is.null(dll_info) | isFALSE(dll_info[["dynamicLookup"]])) {
            # safely remove previous template dll
            unlink(private$dll_)

          } else {
            private$needs_compile_ <- FALSE
            warning_h(
              "Impossible to replace existing template in current session.",
              "Start a new R session to re-compile your model as a template.",
              "Warning raised because model changed or 'force = TRUE'."
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

    # @description Build a new custom model.
    # @param as_template `logical(1)`.
    # If `TRUE` store the dynamic library as a local template.
    # If `FALSE` recompile the model at every new session.
    # @param force `logical(1)`.
    # If `TRUE` compile even if a cached model already exists.
    build_as_custom_model_ = function(as_template, force) {
      private$needs_compile_ <- TRUE
      # temporary directory where to store source code
      # note: if it's been already cached, mark needs_compile_ false
      private$build_temporary_cache_(as_template, force)

      # asses if still to be compiled after further checks
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

    # @description
    # Check for previous cached dll files and tell if to be recompiled.
    # @return `logical(1)`.
    # Wheter nothing is cached or changes are untracked.
    search_for_cached_dll_ = function() {
      # look for dynamic library name
      cache_dll_files <- list.files(private$cache_dir_,
                                    pattern = .Platform$dynlib.ext,
                                    full.names = TRUE)
      # NOTE: multiple dlls are treated as an error
      if (length(cache_dll_files) == 1) {
        # if all changes are tracked, dll is the one
        if (!untracked_changes(private$source_dir_, cache_dll_files[1])) {
          private$dll_ <- cache_dll_files[1]
          return(FALSE)
        }
      }
      # no dll cached or multiple dlls require recompilation
      TRUE
    },

    # @description
    # Create temporary folder and copy source code there.
    copy_source_code_ = function() {
      # clean all previous file in the cache
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

      # copy helper files in the output directory
      for (f in helper_files(private$source_dir_)) {
        if (file.exists(f)) {
          overwrite(f, file.path(dirname(private$dll_), basename(f)))
        }
      }

      invisible(self)
    },

    # @descrtiption
    # Generate flags used for environment compilation.
    build_compile_flags_ = function(paths = NULL, flags = NULL, libs = NULL) {
      flags <- c(
        # flag to pass the model name for ModelRegistry auto-registration
        paste0("_HOUSEPI_MODULE=", get_template_from_dll(private$dll_)),
        flags
      )
      paths <- c(system.file("include", package = "housepi"), paths)
      housepi_lib <- getLoadedDLLs()[["housepi"]][["path"]]

      c(PKG_CXXFLAGS = paste0(c("$(SHLIB_OPENMP_CXXFLAGS)",
                                paste0("-D", flags),
                                paste0("-I", paths)),
                              collapse = " "),
        PKG_LIBS = paste0(c("$(SHLIB_OPENMP_CXXFLAGS)",
                            if (length(libs)) paste0("-L", libs),
                            housepi_lib),
                          collapse = " "))
    },

    # @description
    # Compile and load the source files inside cached directory.
    # @param quiet `r rd_mdl_quiet_`
    # @param user_paths `character()` Additional include paths.
    # @param user_flags `character()` Additional compiler flags.
    # @param user_libs `character()` Additional linker flags.
    compile_ = function(quiet = TRUE,
                        paths = character(),
                        flags = character(),
                        libs = character()) {
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
      overwrite(custom_dll, private$dll_)
      # drop other files
      unlink(file.path(private$cache_dir_,
                       c("NAMESPACE", "DESCRIPTION", "src", "R")),
             recursive = TRUE)

      invisible(self)
    },

    # @description Initialize new Rcpp pointers to Model and MCMC.
    create_new_pointers_ = function() {
      if (!is.null(private$mcmc_ptr_)) return(invisible())
      load_model_dll_(private$dll_)
      module_name <- get_template_from_dll(private$dll_)
      if (is.na(module_name)) module_name <- private$name_
      private$mcmc_ptr_ <- .engine$create_mcmc_runner(module_name)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description Wheter model is blank.
    # @return `logical(1)`
    is_blank_ = function() {
      private$name_ == "blank"
    },

    # @description Wheter model is template.
    # @return `logical(1)`
    is_template_ = function() {
      !is.na(get_template_dir(private$name_))
    },

    # @description Reference script.
    # @return `character(1)`
    get_reference_ = function() {
      output <- file.path(private$source_dir_, "model.cpp")
      # check if reference file exists
      if (!file.exists(output)) stop_h(
        sprintf("Impossible to find reference file: '%s'.", output),
        "Try to re-generate the model with 'force = TRUE'."
      )
      output
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # UPDATE: I shall drop necessity for refactoring with variant data.
    refactor_data_wide_ = function(x) {
      cpp_data_refactor(x, private$info_data_wide_, private$cols_data_wide_)
    },

    # UPDATE: I shall drop necessity for refactoring with variant data.
    refactor_data_long_ = function(x) {
      cpp_data_refactor(x, private$info_data_long_, private$cols_data_long_)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field name `character(1)`. Model name.
    #' Assign a new `character(1)` to change model name.
    name = function(value) {
      if (!missing(value)) stop_read_only("name")
      spaced_template_choice(private$name_)
    },

    #' @field setup `list`. Last setup used by the build process.
    setup = function(value) {
      if (!missing(value)) stop_read_only("setup")
      .engine$export_setup(private$mcmc_ptr_)
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Initialize model generator defining name of dynamic library.
    #' @param choice `r rd_mdl_choice_`
    #' @param as_template `r rd_mdl_astmpl_`
    #' @param force `r rd_mdl_force_`
    #' @param quiet `r rd_mdl_quiet_`
    #' @param ... `r rd_mdl_flags_`
    initialize = function(choice = NA_character_, as_template = FALSE,
                          force = FALSE, quiet = FALSE, ...) {
      # assign name to the model
      choice <- private$helper_init_choice_(choice)

      # assign quiet field and load user options
      user_options <- private$pre_process_user_options_(quiet, ...)
      # store raw flags for later checks
      private$last_flags_ <- user_options$flags

      if (dir.exists(choice)) {
        private$init_custom_model_(choice, as_template, force)
      } else if (!is.na(get_template_dir(choice))) {
        private$init_template_model_(choice)
      } else {
        stop_h(sprintf("'%s' is not a folder nor a template name.", choice),
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
    }
  )
)
