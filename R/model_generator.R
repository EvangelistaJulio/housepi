#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# model_generator.R
#

############################################################################
############################################################################
############################# Housepi Generator ############################
############################################################################
############################################################################

#' @title Housepi Generator R6 Class
#' @description
#' Help with generating Rcpp model object.
#'
#' @keywords internal
HousepiGenerator <- R6::R6Class(
  classname = "HousepiGenerator",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .name = NA,
    .dll = NA,
    .cache_dir = NA,
    .source_dir = NA,
    .quiet = TRUE,
    .need_compile = FALSE,
    .rcpp = NULL,

    # @description
    # Assert source files exist.
    #
    # @param source_dir `character(1)`\cr
    # The path to source code
    .assert_source_files = function(source_dir) {
      for (file_name in source_files(source_dir)) {
        if (!file.exists(file_name)) {
          stop_h(
            sprintf("File \"%s\" was not found.", file_name),
            "Check your model directory."
          )
        }
      }
      invisible(self)
    },

    # @description
    # Assert extra allocations in source code match calls.
    #
    # @param source_dir `character(1)`\cr
    # The path to source code
    .assert_extra = function(source_dir) {
      # use source file as a reference
      refer_file <- file.path(source_dir, "model.cpp")

      # dictionary of allocate\call functions
      checks <- list(
        "allocate_extra_house_integer_vector" = ".extra_integer_vector",
        "allocate_extra_house_double_vector" = ".extra_double_vector",
        "allocate_extra_house_double_matrix" = ".extra_double_matrix"
      )

      missing_set <- c()
      unused_set <- c()
      for (i in seq_along(checks)) {
        # regular expression
        reg <- '%s\\s*\\(\\s*["\'](.*?)["\']'
        # .extra_<type>("string")
        reg_called <- sprintf(reg, checks[[i]])
        # allocate_extra_house_<type>("string")
        reg_allocated <- sprintf(reg, names(checks)[[i]])

        # set of extra variables
        requirements <- template_find_requirement(refer_file, reg_called)
        defined <- template_find_requirement(refer_file, reg_allocated)

        # point out extra variables
        unused_set <- c(unused_set,
                        template_find_unused(requirements, defined))
        # find missing extra variables
        missing_set <- c(missing_set,
                         template_find_missing(requirements, defined))
      }

      # warn for unused extra variables
      if (length(unused_set) > 0) {
        warning_h(
          sprintf("Allocated extra variable%s not used in the model: %s",
                  ifelse(length(unused_set) > 1, "s", ""),
                  character_to_string_list(unused_set)),
          find_closest_match(unused_set, requirements)
        )
      }
      # error for not allocated extra variables
      if (length(missing_set) > 0) {
        stop_h(
          sprintf("Missing allocation for extra variable%s in the model:",
                  ifelse(length(missing_set) > 1, "s", "")),
          character_to_string_list(missing_set)
        )
      }

      invisible(self)
    },

    # @description
    # Initialize a new custom model.
    #
    # @param source_dir `character(1)`\cr
    # The path to source code.
    # @param as_template `logical(1)`\cr
    # If `TRUE` store the dynamic library as a local template.
    # If `FALSE` recompile the model at every new session.
    # @param force `logical(1)`\cr
    # If `TRUE` compile even if a cached model already exists.
    .new_custom_model = function(source_dir,
                                 as_template,
                                 force) {
      # asses if model choice requires recompilation
      private$.need_compile <- TRUE

      # force re-compilation or check for changes
      if (file.exists(private$.dll)) {
        if (!(force || untracked_changes(source_dir, private$.dll))) {
          private$.need_compile <- FALSE
        }
      }

      # temporary directory where to store source code
      if (private$.need_compile) {
        private$.cache_dir <- file.path(tempdir(), template_name(private$.name))

        if (!dir.exists(private$.cache_dir)) {
          if (!dir.create(private$.cache_dir, recursive = TRUE,
                          showWarnings = FALSE)) {
            stop_h(
              "Impossibile to create a cache directory for compilation.",
              "Check your user permissions."
            )
          }
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
          private$.need_compile <- private$.search_for_cached_dll(source_dir)
        }
      }

      # asses if still to be compiled after further checks
      if (private$.need_compile) {
        # change dll name if not as template
        if (!as_template) {
          private$.dll <- file.path(private$.cache_dir,
                                    dll_file_ext(random_name(source_dir)))
        }
        # source directory must contain source files
        private$.assert_source_files(source_dir)
        # organize cache directory to perform compilation
        private$.assert_extra(source_dir)
        private$.copy_source_code(source_dir)
      }

      invisible(self)
    },

    # @description
    # Check for previous cached dll files and tell if to be recompiled.
    #
    # @param source_dir
    # The path to source code.
    #
    # @return `logical(1)`\cr
    # Wheter nothing is cached or changes are untracked.
    .search_for_cached_dll = function(source_dir) {
      # look for dynamic library name
      cache_dll_files <- list.files(private$.cache_dir,
                                    pattern = dll_file_ext(),
                                    full.names = TRUE)
      # NOTE: at the moment multiple dlls are treated as an error
      if (length(cache_dll_files) == 1) {
        # if all changes are tracked, dll is the one
        if (!untracked_changes(source_dir, cache_dll_files[1])) {
          private$.dll <- cache_dll_files[1]
          return(FALSE)
        }
      }

      # no dll cached or multiple dlls require recompilation
      return(TRUE)
    },

    # @description
    # Create temporary folder and copy source code there.
    #
    # @param source_dir
    # The path to source code.
    .copy_source_code = function(source_dir) {

      # clean all previous file in the cache
      unlink(list.files(private$.cache_dir, full.names = TRUE),
             recursive = TRUE)

      # 'housepi' template package directory
      housepi_dir <- file.path(system.file(package = "housepi"), "template")

      # copy folders from template
      for (fldr in file.path(housepi_dir, c("src", "R"))) {
        overwrite(fldr, private$.cache_dir, recursive = TRUE)
      }

      # copy single package files
      for (f in file.path(housepi_dir, c("NAMESPACE", "DESCRIPTION"))) {
        overwrite(f, file.path(private$.cache_dir, basename(f)))
      }

      # copy model specific source files in the new folder
      for (f in source_files(source_dir)) {
        overwrite(f, file.path(private$.cache_dir, "src", basename(f)))
      }

      # copy helper files in the output directory
      for (f in helper_files(source_dir)) {
        if (file.exists(f)) {
          overwrite(f, file.path(dirname(private$.dll), basename(f)))
        }
      }

      # Create a reference script
      self$save_reference()

      invisible(self)
    },

    # @description
    # Compile and load the source files inside cached directory.
    # Note: I am tricking pkgbuild::compile_dll and R CMD INSTALL.
    # To avoid recompiling immutable files I just touch them.
    .compile = function() {
      # source folder within cache directory
      cache_dir_src <- file.path(private$.cache_dir, "src")
      # remove object file
      unlink(file.path(cache_dir_src, "model.o"))

      # touch installation files
      touch_time <- Sys.time()
      for (f in list.files(file.path(cache_dir_src), full.names = TRUE)) {
        Sys.setFileTime(f, touch_time)
      }
      # touch as source model file was more recent
      Sys.setFileTime(file.path(cache_dir_src, "model.cpp"), touch_time + 1)

      # compile the dynamic library
      pkgbuild::compile_dll(private$.cache_dir,
                            quiet = self$quiet,
                            debug = FALSE)
      private$.need_compile <- FALSE

      # copy DLL outside of the source directory
      from_file <- file.path(cache_dir_src, dll_file_ext("housepienvironment"))
      overwrite(from_file, private$.dll)

      invisible(self)
    },

    # @description
    # Return new instance of Rcpp Model.
    .new_rcpp_model = function() {
      # load dynamic library
      dyn.load(private$.dll)
      # expose internally the module
      package <- tools::file_path_sans_ext(basename(private$.dll))
      m <- Rcpp::Module("housepi_model_core", PACKAGE = package)

      return(m$ModelHousepi$new())
    },

    # @description
    # When x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
    # each field, with the name and value.
    deep_clone = function(name, value) {
      # make a new instance of rcpp model object
      if (name == ".rcpp") {
        output <- private$.new_rcpp_model()
        output$import_setup(
          data.frame(lapply(value$export_setup(), as.character))
        )
        output$import_dataset(
          data.frame(lapply(value$export_dataset(), as.character))
        )
        return(output)
      }

      # for all other fields, just return the value
      return(value)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field name `character(1)`\cr
    #' Model name.
    name = function(value) {
      if (missing(value)) return(private$.name)
      stop_h("'name' is read only.")
    },

    #' @field cache_dir `character(1)`\cr
    #' Cache directory where the compilation is stored.
    cache_dir = function(value) {
      if (missing(value)) return(private$.cache_dir)
      stop_h("'cache_dir' is read only.")
    },

    #' @field is_blank `logical(1)`\cr
    #' Wheter model is blank.
    is_blank = function(value) {
      if (missing(value)) return(private$.name == "blank")
      stop_h("'is_blank' is read only.")
    },

    #' @field quiet `logical(1)`\cr
    #' Wheter quiet option is on.
    quiet = function(value) {
      if (missing(value)) return(private$.quiet)
      stop_h("'quiet' is read only.")
    },

    #' @field path `character(1)`\cr
    #' Envirnoment path.
    path = function(value) {
      if (missing(value)) {
        # reference script
        output <- dirname(private$.dll)
        # check if reference file exists
        if (!dir.exists(output)) {
          stop_h(
            sprintf("Impossible to find environment path: \"%s\".", output),
            "Try to re-generate the model with 'force = TRUE'."
          )
        }
        return(output)
      }
      stop_h("'path' is read only.")
    },

    #' @field source_dir `character(1)`\cr
    #' Directory containing the orginal source code.
    source_dir = function(value) {
      if (missing(value)) return(private$.source_dir)
      stop_h("'source_dir' is read only.")
    },

    #' @field reference `character(1)`\cr
    #' Reference script.
    reference = function(value) {
      if (missing(value)) {
        if (self$is_blank) {
          output <- file.path(self$source_dir, "model.cpp")
        } else {
          # reference script
          output <- file.path(self$path, ".refer")
          # check if reference file exists
          if (!file.exists(output)) {
            stop_h(
              sprintf("Impossible to find reference file: \"%s\".", output),
              "Try to re-generate the model with 'force = TRUE'."
            )
          }
        }
        return(output)
      }
      stop_h("'reference' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' Initialize model generator defining name of dynamic library.
    #'
    #' @param choice `character(1)`\cr
    #' Choose a model providing one between:
    #' - "blank" template name;
    #' - template name;
    #' - path/to/folder with implementation;
    #' - path/to/folder with previously compiled model.
    #' @param as_template `logical(1)`\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' If `FALSE` recompile the model at every new session.
    #' @param force `logical(1)`\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param name `chracter(1)`\cr
    #' Custom name of the model (default is inferred from `choice`).
    #' @param quiet `logical(1)`\cr
    #' If `TRUE` hide `pkgbuild` compilation commands.
    initialize = function(choice, as_template = FALSE, force = FALSE,
                          name = NA, quiet = TRUE) {
      # check model choice
      if (is.na(choice)) {
        choice <- "blank"
        warning_h(
          "No model specified. Falling back to \"blank\" model.",
          "To show available templates, run function 'template_household_model'."
        )

      } else if (choice %in% list(".")) {
        # prohibit directories
        stop_h(
          sprintf("Invalid choice \"%s\" for housepi model.", choice),
          "Move your implementation in a different folder."
        )
      }

      # assign name to the model
      private$.name <- ifelse(is.na(name), basename(choice), name)

      # assign quiet field
      private$.quiet <- isTRUE(quiet)

      if (choice == "blank") {
        # 1: Blank Model (Developer Mode) 
        private$.dll <- file.path(system.file(package = "housepi"),
                                  "template", "src",
                                  dll_file_ext("housepienvironment"))
        private$.source_dir <- dirname(private$.dll)

      } else if (dir.exists(choice)) {
        # 2: Model Path Contains Source Code
        private$.source_dir <- choice
        private$.dll <- template_dll_file(choice)
        private$.new_custom_model(choice, as_template, force)

      } else if (choice %in% list.files(template_models_dir())) {
        # 3: Plain Template Model Name
        private$.source_dir <- file.path(template_models_dir(), choice)
        private$.dll <- template_dll_file(private$.source_dir)

      } else {
        # 4: Error
        stop_h(
          sprintf("\"%s\" is not a folder nor a template name.", choice),
          "To show available templates, run function 'template_household_model'."
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Generate the Rcpp model object.
    #' If necessary, it recompiles only the model script.
    generate = function() {
      if (private$.need_compile) {
        # compile model
        cli::cli_alert_info("Compiling {.pkg {private$.name}}")
        private$.compile()
        private$.rcpp <- private$.new_rcpp_model()

      } else if (is.null(private$.rcpp)) {
        # load if not already generated
        cli::cli_alert_info("Loading {.pkg {private$.name}}")
        private$.rcpp <- private$.new_rcpp_model()
      }

      invisible(self)
    },

    #' @description
    #' Save the reference script.
    save_reference = function() {
      # reference script
      refer_file <- file.path(dirname(private$.dll), ".refer")
      overwrite(file.path(private$.cache_dir, "src", "model.cpp"),
                refer_file)

      invisible(self)
    }
  )
)
