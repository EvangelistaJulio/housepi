#
# generator.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title Housepi Generator Class
#' @description
#' Help with generating model Rcpp object.
#' 
#' @noRd
HousepiGeneratorClass <- R6::R6Class(
  classname = "HousepiGeneratorClass",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .name = NA,
    .dll = NA,
    .quiet = TRUE,
    .need_compile = FALSE,

    #' @description
    #' Create temporary folder and copy source code there.
    copy_source_code = function(source_dir) {

      # source directory must contain source files
      assert_source_files(source_dir)

      # clean all previous file in the cache
      unlink(list.files(self$cache_dir, full.names = TRUE),
             recursive = TRUE)

      # 'housepi' template package directory
      housepi_dir <- file.path(system.file(package = "housepi"), "template")

      # copy folders from template
      for (fldr in file.path(housepi_dir, c("src", "R"))) {
        overwrite(fldr, self$cache_dir, recursive = TRUE)
      }

      # copy single package files
      for (f in file.path(housepi_dir, c("NAMESPACE", "DESCRIPTION"))) {
        overwrite(f, file.path(self$cache_dir, basename(f)))
      }

      # copy model specific source files in the new folder
      for (f in source_files(source_dir)) {
        overwrite(f, file.path(self$cache_dir, "src", basename(f)))
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

    #' @description
    #' Remove all cache files.
    clear_cache = function() {
      cache_files <- file.path(self$cache_dir,
                               c("src", "R", "NAMESPACE", "DESCRIPTION"))
      unlink(cache_files, recursive = TRUE)

      invisible(self)
    },

    #' @description
    #' Load model dynamic library.
    load = function() {
      # load dynamic library
      dyn.load(private$.dll)

      # expose internally the module
      package <- tools::file_path_sans_ext(basename(private$.dll))
      m <- Rcpp::Module("model_housepi_rcpp", PACKAGE = package)

      # assign the model object
      self$rcpp <- m$ModelHousepiRcpp$new()

      invisible(self)
    },

    #' @description
    #' Compile and load the source files inside cached directory.
    #' Note: I am tricking pkgbuild::compile_dll and R CMD INSTALL.
    #' To avoid recompiling immutable files I just touch them.
    compile = function() {
      # source folder within cache directory
      cache_dir_src <- file.path(self$cache_dir, "src")

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
      pkgbuild::compile_dll(self$cache_dir,
                            quiet = private$.quiet, debug = FALSE)
      private$.need_compile <- FALSE

      # copy DLL outside of the source directory
      from_file <- file.path(cache_dir_src, dll_file_ext("housepienvironment"))
      overwrite(from_file, private$.dll)

      #load the compiled package
      private$load()

      invisible(self)
    }
  ),

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(

    #' @field path (`character(1)`)\cr
    #' Model name.
    model_name = function(value) {
      if (missing(value)) {
        return(private$.name)
      } else {
        stop.h("'model_name' is read only.")
      }
    },

    #' @field path (`logical(1)`)\cr
    #' Wheter model is blank.
    is_blank = function(value) {
      if (missing(value)) {
        return(private$.name == "blank")
      } else {
        stop.h("'is_blank' is read only.")
      }
    },

    #' @field path (`character(1)`)\cr
    #' Envirnoment path.
    path = function(value) {
      if (missing(value)) {
        # reference script
        output <- dirname(private$.dll)
        # check if reference file exists
        if (!dir.exists(output)) {
          stop.h(
            paste(
             sprintf("Impossible to find environment path: \"%s\".", output),
             "Try to re-generate the model with `force = TRUE`."
            )
          )
        }
        return(output)
      } else {
        stop.h("'path' is read only.")
      }
    },

    #' @field reference (`character(1)`)\cr
    #' Reference script.
    reference = function(value) {
      if (missing(value)) {
        # reference script
        output <- file.path(self$path, ".refer")
        # check if reference file exists
        if (!file.exists(output)) {
          stop.h(
            paste(
             sprintf("Impossible to find reference file: \"%s\".", output),
             "Try to re-generate the model with `force = TRUE`."
            )
          )
        }
        return(output)
      } else {
        stop.h("'reference' is read only.")
      }
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    rcpp = NULL,
    cache_dir = NA,

    #' @description
    #' Initialize model generator defining name of dynamic library.
    #' 
    #' @param model_choice (`character(1)`)\cr
    #' One between:
    #' - template name;
    #' - path/to/folder with implementation;
    #' - path/to/folder with previously compiled model;
    #' - "blank" template name.
    #' @param quiet (`logical(1)`)\cr
    #' If `TRUE` hide compilation commands.
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param as_template (`logical(1)`)\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' If `FALSE` recompile the model at every new session.
    initialize = function(model_choice,
                          quiet = TRUE,
                          force = FALSE,
                          as_template = FALSE) {

      if (model_choice == ".") {
        # prohibit current directory
        stop.h(
          paste(
            "Invalid folder name \".\" for housepi model.",
            "Move your implementation in a different folder."
          )
        )
      }

      # assign quiet field
      private$.quiet <- isTRUE(quiet)

      # interpret model choice
      if (isTRUE(model_choice == "blank")) {
        #===== Developer Mode =====#
        private$.name <- "blank"
        private$.dll <- file.path(system.file(package = "housepi"),
                                 "template", "src",
                                 dll_file_ext("housepienvironment"))

      } else if (model_choice %in% list.files(template_models_dir())) {
        #===== Plain Template Model Name =====#
        private$.name <- model_choice
        model_choice <- file.path(template_models_dir(), model_choice)
        private$.dll <- template_dll_file(model_choice)

      } else if (dir.exists(model_choice)) {
        #===== Model Path Containing Source Code =====#
        private$.name <- basename(model_choice)
        private$.dll <- template_dll_file(model_choice)

        # asses if to be compiled
        if ((!file.exists(private$.dll)) || isTRUE(force)) {
          # standard temporary directory where to store source code
          self$cache_dir <- file.path(tempdir(), template_name(private$.name))
          if (!dir.exists(self$cache_dir)) {
            if (!dir.create(self$cache_dir,
                            recursive = TRUE,
                            showWarnings = FALSE)) {
              stop.h(
                paste(
                  "Impossibile to create a cache directory for compilation.",
                  "Check your user permissions."
                )
              )
            }
          }

          # check for previous cached dll files
          cache_dll_files <- list.files(self$cache_dir,
                                        pattern = dll_file_ext())

          # define dynamic library name
          if ((length(cache_dll_files) == 1) && isFALSE(force)) {
            private$.dll <- cache_dll_files[1]

          } else {
            # change dll name if not as template
            if (isFALSE(as_template)) {
              private$.dll <- file.path(self$cache_dir,
                                       dll_file_ext(random_name(model_choice)))
            }
            private$.need_compile <- TRUE
            # organize cache directory to perform compilation
            private$copy_source_code(model_choice)
          }
        }

      } else {
        #===== ERROR =====#
        stop.h(
          paste(
            sprintf("\"%s\" is not a folder nor a template name.", model_choice),
            "Run `household_model_template` to see available template names."
          )
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Name of model dynamic library. Extract name of the dynamic library.
    generate = function() {
      if (isTRUE(private$.need_compile)) {
        cli::cli_alert_info("Compiling {.pkg {private$.name}}")
        private$compile()

      } else if (is.null(self$rcpp)) {
        # load if not already generated
        cli::cli_alert_info("Loading {.pkg {private$.name}}")
        private$load()
      }

      invisible(self)
    },

    #' @description
    #' Genera reference script.
    save_reference = function() {
      # reference script
      refer_file <- file.path(dirname(private$.dll), ".refer")
      overwrite(file.path(self$cache_dir, "src", "model.cpp"),
                refer_file)

      invisible(self)
    }
  )
)

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @description
#' Re-compiles all required objects of the housepi library.
#' Passed to .onLoad in file zzz.R
#' 
#' @param pkgname
#' Name of the package.
#' @param force
#' Force compilation even if it already exists.
#' 
#' @noRd
re_compile_library <- function(pkgname,
                               force = FALSE) {
  # compile the core
  housepi_models_dir <- file.path(system.file(package = pkgname), "template")
  pkgbuild::compile_dll(housepi_models_dir, debug = FALSE)

  # compile each model
  models_dir <- file.path(system.file(package = pkgname), "models")
  for (m in list.files(models_dir, full.names = TRUE)) {
    if ((!file.exists(template_dll_file(m))) || isTRUE(force)) {
      HousepiGeneratorClass$new(m,
                                force = force,
                                quiet = FALSE,
                                as_template = TRUE)$generate()
    }
  }
}
