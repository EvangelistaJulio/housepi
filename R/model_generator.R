#
# model_generator.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

# @title Housepi Generator R6 Class
# @description
# Help with generating Rcpp model object.
#
# @noRd
HousepiGenerator <- R6::R6Class(
  classname = "HousepiGenerator",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .name = NA,
    .dll = NA,
    .cache_dir = NA,
    .quiet_pkgbuild = TRUE,
    .need_compile = FALSE,
    .rcpp = NULL,

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
      if (file.exists(private$.dll)) {
        # force template re-compilation or check for changes
        if (force | untracked_changes(source_dir, private$.dll)) {
          private$.need_compile <- TRUE
        } else {
          private$.need_compile <- FALSE
        }
      } else {
        private$.need_compile <- TRUE
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

            if (is.null(dll_info) | (!dll_info[["dynamicLookup"]])) {
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

        } else {
          if (force) {
             private$.need_compile <- TRUE
          } else {
            # search cached dll and check if to be recompiled 
            private$.search_for_cached_dll(source_dir)
          }
        }
      }

      # asses if still to be compiled after further checks
      if (private$.need_compile) {
        # change dll name if not as template
        if (!as_template) {
          private$.dll <- file.path(private$.cache_dir,
                                    dll_file_ext(random_name(source_dir)))
        }
        # organize cache directory to perform compilation
        private$.copy_source_code(source_dir)
      }
    },

    # @description
    # Check for previous cached dll files and tell if to be recompiled.
    #
    # @param source_dir
    # The path to source code.
    .search_for_cached_dll = function(source_dir) {
      # look for dynamic library name
      cache_dll_files <- list.files(private$.cache_dir,
                                    pattern = dll_file_ext(),
                                    full.names = TRUE)

      if (length(cache_dll_files) == 1) {
        # check if it needs recompile
        # FIXME: #UPDATE: can look for the correct one in the list...
        private$.need_compile <- untracked_changes(
          source_dir, cache_dll_files[1]
        )
        if (!private$.need_compile) {
          private$.dll <- cache_dll_files[1]
        }
      } else {
        # no dll cached or multiple dlls
        private$.need_compile <- TRUE
      }

      invisible(self)
    },

    # @description
    # Create temporary folder and copy source code there.
    #
    # @param source_dir
    # The path to source code.
    .copy_source_code = function(source_dir) {

      # source directory must contain source files
      assert_source_files(source_dir)

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
      # touch as if more recent the source model file
      Sys.setFileTime(file.path(cache_dir_src, "model.cpp"), touch_time + 1)

      # compile the dynamic library
      pkgbuild::compile_dll(private$.cache_dir,
                            quiet = private$.quiet_pkgbuild,
                            debug = FALSE)
      private$.need_compile <- FALSE

      # copy DLL outside of the source directory
      from_file <- file.path(cache_dir_src, dll_file_ext("housepienvironment"))
      overwrite(from_file, private$.dll)

      invisible(self)
    },

    # @description
    # When x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
    # each field, with the name and value.
    deep_clone = function(name, value) {
      if (name == ".rcpp") {
        # make a new instance of rcpp model object
        output <- private$.new_rcpp_model()
        output$import_setup(
          data.frame(lapply(value$export_setup(), as.character))
        )
        output$import_dataset(
          data.frame(lapply(value$export_dataset(), as.character))
        )
        return(output)

      } else {
        # for all other fields, just return the value
        return(value)
      }
    }
  ),

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    # @field model_name `character(1)`\cr
    # Model name.
    model_name = function(value) {
      if (missing(value)) {
        return(private$.name)
      }
      stop_h("'model_name' is read only.")
    },
    
    # @field cache_dir `character(1)`\cr
    # Cache directory where the compilation is stored.
    cache_dir = function(value) {
      if (missing(value)) {
        return(private$.cache_dir)
      }
      stop_h("'cache_dir' is read only.")
    },

    # @field is_blank `logical(1)`\cr
    # Wheter model is blank.
    is_blank = function(value) {
      if (missing(value)) {
        return(private$.name == "blank")
      }
      stop_h("'is_blank' is read only.")
    },

    # @field path `character(1)`\cr
    # Envirnoment path.
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

    # @field reference `character(1)`\cr
    # Reference script.
    reference = function(value) {
      if (missing(value)) {
        # reference script
        output <- file.path(self$path, ".refer")
        # check if reference file exists
        if (!file.exists(output)) {
          stop_h(
            sprintf("Impossible to find reference file: \"%s\".", output),
            "Try to re-generate the model with 'force = TRUE'."
          )
        }
        return(output)
      }
      stop_h("'reference' is read only.")
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    # @description
    # Initialize model generator defining name of dynamic library.
    #
    # @param model_choice `character(1)`\cr
    # One between:
    # - "blank" template name;
    # - template name;
    # - path/to/folder with implementation;
    # - path/to/folder with previously compiled model.
    # @param as_template `logical(1)`\cr
    # If `TRUE` store the dynamic library as a local template.
    # If `FALSE` recompile the model at every new session.
    # @param force `logical(1)`\cr
    # If `TRUE` compile even if a cached model already exists.
    # @param quiet_pkgbuild `logical(1)`\cr
    # If `TRUE` hide `pkgbuild` compilation commands.
    initialize = function(model_choice,
                          as_template = FALSE,
                          force = FALSE,
                          quiet_pkgbuild = TRUE) {
      # assign quiet_pkgbuild field
      private$.quiet_pkgbuild <- isTRUE(quiet_pkgbuild)

      if (model_choice == "blank") {
        #===== Blank Model (Developer Mode) =====#
        private$.name <- model_choice
        private$.dll <- file.path(system.file(package = "housepi"),
                                  "template", "src",
                                  dll_file_ext("housepienvironment"))

      } else if (model_choice %in% list.files(template_models_dir())) {
        #===== Plain Template Model Name =====#
        private$.name <- model_choice
        source_dir <- file.path(template_models_dir(), model_choice)
        private$.dll <- template_dll_file(source_dir)

      } else if (model_choice %in% list(".")) {
        # prohibit current directory
        stop_h(
          sprintf("Invalid folder name \"%s\" for housepi model.",
                  model_choice),
          "Move your implementation in a different folder."
        )

      } else if (dir.exists(model_choice)) {
        #===== Model Path Contains Source Code =====#
        private$.name <- basename(model_choice)
        private$.dll <- template_dll_file(model_choice)
        private$.new_custom_model(model_choice, as_template, force)

      } else {
        #===== Error =====#
        stop_h(
          sprintf("\"%s\" is not a folder nor a template name.", model_choice),
          "To see available templates, run function 'template_household_model'."
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Name of model dynamic library. Extract name of the dynamic library.
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

    # @description
    # Save reference script.
    save_reference = function() {
      # reference script
      refer_file <- file.path(dirname(private$.dll), ".refer")
      overwrite(file.path(private$.cache_dir, "src", "model.cpp"),
                refer_file)

      invisible(self)
    }
  )
)
