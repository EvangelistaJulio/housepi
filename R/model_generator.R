## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom R6 R6Class
#' @importFrom pkgbuild compile_dll
NULL

# @description
# Assert source files exist.
generator_assert_source_files <- function(source_dir) {
  for (file_name in source_files(source_dir)) {
    if (!file.exists(file_name)) stop_h(
      sprintf("File '%s' was not found.", file_name),
      "Check your model directory."
    )
  }
}

# @description
# Check whether source model contains a literal macro/token.
generator_assert_allows_flag <- function(source_dir, flag) {
  candidates <- c(file.path(source_dir, "src", "model.cpp"),
                  file.path(source_dir, "model.cpp"))
  model_cpp <- candidates[file.exists(candidates)][1]
  if (is.na(model_cpp)) return(FALSE)
  lines <- readLines(model_cpp, warn = FALSE)
  patt <- paste0("\\b", gsub("([\\W])", "\\\\\\1", flag), "\\b")
  any(grepl(patt, lines, perl = TRUE))
}

# @description
# Process user options for compilation.
generator_process_user_options <- function(...) {
  user_paths <- character()
  user_flags <- character()
  user_linkers <- character()
  user_options <- list(...)

  # process time discrete option
  if ("discr" %in% names(user_options) && isTRUE(user_options$discr)) {
    user_flags <- c(user_flags, "_HOUSEPI_DISCRETE_TIME")
  }

  list("user_paths" = user_paths,
       "user_flags" = user_flags,
       "user_linkers" = user_linkers)
}

# @description
# Wheter model choice is interpretable.
generator_init_helper_choice <- function(choice) {
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

  choice_
}

# @description Reference script.
# @param seed `integer(1)`
# @return `integer(1)`
generator_adjust_seed <- function(seed) {
  max_r_int <- .Machine$integer.max
  if (is.na(seed)) seed <- sample.int(max_r_int, 1)
  else seed <- abs(as.integer(seed))

  # control that seed is 32bit (for custom seed)
  if (seed > max_r_int) stop_h(
    "Supplied seed is not supported.",
    sprintf("Maximum value for an integer in R is %.0f.", max_r_int)
  )
  seed
}

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
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    # @description
    # When x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
    # each field, with the name and value.
    # UPDATE: avoid re-importing dataset if not necessary using shared pointers
    deep_clone = function(name, value) {
      # make a new object of rcpp model class
      if (name == ".rcpp") {
        output <- NULL
        if (!is.null(value)) {
          output <- private$.new_rcpp_model()
          # note: no factorization is needed here
          output$import_setup(char_setup_input_(value$export_setup()))
          output$import_data_wide(char_data_input_(value$export_data_wide()))
          output$import_data_long(char_data_input_(value$export_data_long()))
        }
        return(output)
      }
      # for all other fields, use generalized helper for deep clone
      helper_clone(value)
    },

    .name = NA_character_,
    .dll = NA_character_,
    .cache_dir = NA_character_,
    .source_dir = NA_character_,
    .quiet = FALSE,
    .need_compile = FALSE,
    .rcpp = NULL,
    .last_user_flags = character(),

    # @field .infec_status_set `integer()`\cr
    # All positive status extracted from data.
    .infec_status_set = integer(),

    # @field .cols_data_wide `character()`\cr
    # Headers of last imported data wide.
    .cols_data_wide = NULL,

    # @field .info_data_wide `list`\cr
    # List of levels of the factors in data wide.
    .info_data_wide = NULL,

    # @field .cols_data_long `character()`\cr
    # Headers of the last data long imported.
    .cols_data_long = NULL,

    # @field .info_data_long `list`\cr
    # List of levels of the factors in the data long.
    .info_data_long = NULL,

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assign source directory and filter requested flags according to source.
    # If a flag is not present in model.cpp, raise a warning.
    .store_source_and_filter_flags = function(source_dir) {
      # asses if model flags actually make sense
      keep <- character()
      for (f in private$.last_user_flags) {
        if (generator_assert_allows_flag(source_dir, f)) keep <- c(keep, f)
        else if (f == "_HOUSEPI_DISCRETE_TIME") warning_h(
          sprintf("No visible discrete time option for model '%s'.", self$name)
        )
        else warning_h(
          sprintf("Ignoring unused flag '%s' not found in model source.", f),
          "Only flags present in model.cpp are applied at compile time."
        )
      }
      # store just effective flags
      private$.last_user_flags <- keep

      # store source directory for later checks
      private$.source_dir <- source_dir

      invisible(self)
    },

    # @description
    # Wheter user flags are the same as existing dll.
    .wheter_same_flags_as_dll = function() {
      if (private$.is_blank()) {
        flags_dll <- c()
      } else {
        m <- get_rcpp_module_from_dll(private$.dll)
        y <- m$housepi_compile_flags()
        flags_dll <- strsplit(y, ":")[[1]]
      }
      setequal(private$.last_user_flags, flags_dll)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Initialize to "blank" model.
    .initialize_blank_model = function() {
      private$.last_user_flags <- character()
      private$.source_dir <- system.file("template/src", package = "housepi")
      # NOTE: blank model dll is never used for loading
      private$.dll <- file.path(private$.source_dir,
                                dll_file_ext("housepienvironment"))
      invisible(self)
    },

    # @description
    # Initialize custom model loading path to model source code.
    .initialize_custom_model = function(choice, as_template, force) {
      # UPDATE: storing should happen after assert source!
      # otherwise idiot warning...
      private$.store_source_and_filter_flags(choice)
      private$.dll <- template_dll_file(choice)
      private$.new_custom_model(as_template, force)
      invisible(self)
    },

    # @description
    # Initialize template model loading template model name.
    .initialize_template_model = function(choice, force) {
      private$.store_source_and_filter_flags(get_template_dir(choice))
      private$.dll <- template_dll_file(private$.source_dir)
      private$.new_custom_model(FALSE, force)
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Predispose temporare fodler for the installation.
    .new_temporary_cache = function(as_template, force) {
      private$.cache_dir <- file.path(tempdir(), template_name(private$.name))

      if (!dir.exists(private$.cache_dir)) {
        if (!dir.create(private$.cache_dir, recursive = TRUE,
                        showWarnings = FALSE)) stop_h(
          "Impossibile to create a cache directory for compilation.",
          "Check your user permissions."
        )
      }

      if (as_template) {
        if (file.exists(private$.dll)) {
          # check if dll is loaded
          namespace <- template_name(private$.name)
          dll_info <- base::getLoadedDLLs()[[namespace]]

          if (is.null(dll_info) | isFALSE(dll_info[["dynamicLookup"]])) {
            # safely remove previous template dll
            unlink(private$.dll)

          } else {
            private$.need_compile <- FALSE
            warning_h(
              "Impossible to replace existing template in current session.",
              "Start a new R session to re-compile your model as a template.",
              "Warning raised because model changed or 'force = TRUE'."
            )
          }
        }

      } else if (force) {
        private$.need_compile <- TRUE

      } else {
        # search cached dll and check if to be recompiled
        private$.need_compile <- private$.search_for_cached_dll()
      }
      invisible(self)
    },

    # @description Initialize a new custom model.
    # @param as_template `logical(1)`.
    # If `TRUE` store the dynamic library as a local template.
    # If `FALSE` recompile the model at every new session.
    # @param force `logical(1)`.
    # If `TRUE` compile even if a cached model already exists.
    .new_custom_model = function(as_template, force) {
      # force re-compilation or check for changes
      private$.need_compile <- TRUE
      if (!force &&
            file.exists(private$.dll) &&
            !untracked_changes(private$.source_dir, private$.dll) &&
            private$.wheter_same_flags_as_dll()) {
        private$.need_compile <- FALSE
      }

      # temporary directory where to store source code
      if (private$.need_compile) {
        private$.new_temporary_cache(as_template, force)
      }

      # asses if still to be compiled after further checks
      if (private$.need_compile) {
        # change dll name if not as template
        if (!as_template) private$.dll <- file.path(
          private$.cache_dir, dll_file_ext(random_name(private$.source_dir))
        )
        # source directory must contain source files
        generator_assert_source_files(private$.source_dir)
        # organize cache directory to perform compilation
        private$.copy_source_code()
      }

      invisible(self)
    },

    # @description
    # Check for previous cached dll files and tell if to be recompiled.
    # @return `logical(1)`.
    # Wheter nothing is cached or changes are untracked.
    .search_for_cached_dll = function() {
      # look for dynamic library name
      cache_dll_files <- list.files(private$.cache_dir,
                                    pattern = .Platform$dynlib.ext,
                                    full.names = TRUE)
      # NOTE: multiple dlls are treated as an error
      if (length(cache_dll_files) == 1) {
        # if all changes are tracked, dll is the one
        if (!untracked_changes(private$.source_dir, cache_dll_files[1])) {
          private$.dll <- cache_dll_files[1]
          return(FALSE)
        }
      }
      # no dll cached or multiple dlls require recompilation
      TRUE
    },

    # @description
    # Create temporary folder and copy source code there.
    .copy_source_code = function() {

      # clean all previous file in the cache
      unlink(list.files(private$.cache_dir, full.names = TRUE),
             recursive = TRUE)

      # 'housepi' template package directory
      housepi_dir <- system.file("template", package = "housepi")

      # copy folders from template
      for (fldr in file.path(housepi_dir, c("src", "R"))) {
        overwrite(fldr, private$.cache_dir, recursive = TRUE)
      }

      # copy single package files
      for (f in file.path(housepi_dir, c("NAMESPACE", "DESCRIPTION"))) {
        overwrite(f, file.path(private$.cache_dir, basename(f)))
      }

      # copy model specific source files in the new folder
      for (f in source_files(private$.source_dir)) {
        overwrite(f, file.path(private$.cache_dir, "src", basename(f)))
      }

      # copy helper files in the output directory
      for (f in helper_files(private$.source_dir)) {
        if (file.exists(f)) {
          overwrite(f, file.path(dirname(private$.dll), basename(f)))
        }
      }

      invisible(self)
    },

    # @description
    # Compile and load the source files inside cached directory.
    # Note: I am tricking pkgbuild::compile_dll and R CMD INSTALL.
    # To avoid recompiling immutable files I just touch them.
    # @param user_paths A `character()` vector of additional include paths.
    # @param user_flags A `character()` vector of additional compiler flags.
    # @param user_linkers A `character()` vector of additional linker flags.
    .compile = function(user_paths = character(),
                        user_flags = character(),
                        user_linkers = character()) {
      # drop previous Rcpp model object
      private$.rcpp <- NULL

      # source folder in cache directory
      cache_dir_src <- file.path(private$.cache_dir, "src")
      # remove object file
      unlink(file.path(cache_dir_src, "model.o"))

      # touch installation files
      touch_time <- Sys.time() - 2
      for (f in list.files(cache_dir_src, full.names = TRUE)) {
        Sys.setFileTime(f, touch_time)
      }
      # touch as source model file was more recent
      Sys.setFileTime(file.path(cache_dir_src, "model.cpp"), touch_time + 1)
      # same for registry about current flags if any
      # NOTE: default installation has empty flags
      if (length(private$.last_user_flags) > 0) {
        Sys.setFileTime(file.path(cache_dir_src, "flags.cpp"), touch_time + 1)
      }

      # compile the dynamic library
      withr::with_envvar(
        get_cpp_envvar_flags(paths = user_paths,
                             flags = user_flags,
                             linkers = user_linkers),
        pkgbuild::compile_dll(private$.cache_dir,
                              quiet = private$.quiet,
                              compile_attributes = FALSE,
                              debug = FALSE)
      )
      private$.need_compile <- FALSE

      # copy DLL outside of the source directory
      from_file <- file.path(cache_dir_src, dll_file_ext("housepienvironment"))
      overwrite(from_file, private$.dll)

      invisible(self)
    },

    # @description Initialize new object of Rcpp Model.
    .new_rcpp_model = function() {
      # expose internally the module
      if (private$.is_blank()) m <- get_rcpp_module_from_model(NULL)
      else m <- get_rcpp_module_from_dll(private$.dll)
      y <- m$HousepiRcppModel$new()
      y
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description Wheter model is blank.
    # @return `logical(1)`
    .is_blank = function() {
      private$.name == "blank"
    },

    # @description Wheter model is template.
    # @return `logical(1)`
    .is_template = function() {
      if (private$.is_blank()) return(TRUE)
      dirname(private$.path()) == template_models_dir()
    },

    # @description Environment path.
    # @return `character(1)`
    .path = function() {
      # reference script
      output <- dirname(private$.dll)
      # check if reference file exists
      if (!dir.exists(output)) stop_h(
        sprintf("Impossible to find environment path: '%s'.", output),
        "Try to re-generate the model with 'force = TRUE'."
      )
      output
    },

    # @description Reference script.
    # @return `character(1)`
    .reference = function() {
      output <- file.path(private$.source_dir, "model.cpp")
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

    .refactor_data_wide = function(x) {
      cpp_data_refactor(x, private$.info_data_wide, private$.cols_data_wide)
    },

    .refactor_data_long = function(x) {
      cpp_data_refactor(x, private$.info_data_long, private$.cols_data_long)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field name `character(1)`. Model name.
    #' Assign a new `character(1)` to change model name.
    name = function(value) {
      if (missing(value)) return(spaced_template_choice(private$.name))
      private$.name <- basename(as.character(value))
    },

    #' @field setup `list`. Last setup used by the build process.
    setup = function(value) {
      if (!missing(value)) stop_h("'setup' is read only.")
      private$.rcpp$export_setup()
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
    initialize = function(choice, as_template = FALSE,
                          force = FALSE, quiet = FALSE, ...) {
      withr::local_preserve_seed()
      # assign name to the model
      choice <- generator_init_helper_choice(choice)
      self$name <- choice

      # assign quiet field
      private$.quiet <- isTRUE(quiet)
      # load user options
      user_options <- generator_process_user_options(...)
      # store raw flags for later checks
      private$.last_user_flags <- user_options$user_flags

      if (choice == "blank") {
        private$.initialize_blank_model()
      } else if (dir.exists(choice)) {
        private$.initialize_custom_model(choice, as_template, force)
      } else if (!is.na(get_template_dir(choice))) {
        private$.initialize_template_model(choice, force)
      } else {
        stop_h(sprintf("'%s' is not a folder nor a template name.", choice),
               find_closest_template_name(choice))
      }

      # finally recompile only the model script, if necessary.
      if (private$.need_compile) do.call(private$.compile, user_options)
      # generate the Rcpp model object.
      if (is.null(private$.rcpp)) private$.rcpp <- private$.new_rcpp_model()

      invisible(self)
    }
  )
)
