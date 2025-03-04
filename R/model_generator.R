#
# model_generator.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

#' @title Housepi Generator R6 Class
#' @description
#' Help with generating Rcpp model object.
#'
#' @noRd
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
    .quiet_pkgbuild = TRUE,
    .need_compile = FALSE,
    .rcpp = NULL,

    # @description
    # Create temporary folder and copy source code there.
    .copy_source_code = function(source_dir) {

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

    #' @field path (`character(1)`)\cr
    #' Model name.
    model_name = function(value) {
      if (missing(value)) {
        return(private$.name)
      }
      stop_h("'model_name' is read only.")
    },

    #' @field path (`logical(1)`)\cr
    #' Wheter model is blank.
    is_blank = function(value) {
      if (missing(value)) {
        return(private$.name == "blank")
      }
      stop_h("'is_blank' is read only.")
    },

    #' @field path (`character(1)`)\cr
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

    #' @field reference (`character(1)`)\cr
    #' Reference script.
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
    },

    #' @field dataset (`data.frame`)\cr
    #' Dataset with individual level information of current dataset.
    #' If the case, it also return the associated augmented variables.
    dataset = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export dataset before model build.")
        }
        return(private$.rcpp$export_dataset())
      }
      stop_h("'dataset' is read only.")
    },

    #' @field setup (`list()`)\cr
    #' Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export setup before model build.")
        }
        return(private$.rcpp$export_setup())
      }
      stop_h("'setup' is read only.")
    },

    #' @field (`matrix()` | `list()` of `matrix()`)
    #' Export the traces of last inference call.
    traces = function(value) {
      if (missing(value)) {
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export traces before inference.")

        } else {
          # load traces from last inference
          df <- private$.rcpp$export_traces()

          if (dim(df)[1] == 0) {
            # note: this should never be possible...
            stop_h("Traces are empty.")

          } else if ("CHAIN__" %in% colnames(df)) {
            df_list <- lapply(
              split(seq_len(nrow(df)), df[, "CHAIN__"]),
              function(idx) df[idx, colnames(df) != "CHAIN__", drop = FALSE]
            )
            return(df_list)

          } else {
            return(df)
          }
        }
      }
      stop_h("'traces' is read only.")
    },

    #' @field (`matrix()` | `list()` of `matrix()`)
    #' Export the traces of last inference call without loglikelihood.
    traces_no_loglike = function(value) {
      if (missing(value)) {
        trs <- self$traces
        if (is.list(trs)) {
          return(lapply(trs, function(x) {
            x <- x[, colnames(x) != "loglikelihood", drop = FALSE]
          }))
        } else {
          return(trs[, colnames(trs) != "loglikelihood", drop = FALSE])
        }
      }
      stop_h("'traces_no_loglike' is read only.")
    },

    #' @field traces_as_coda (`coda::mcmc` | `coda::mcmc.list`)\cr
    #' Traces of last inference as a `coda` object.
    traces_as_coda = function(value) {
      if (missing(value)) return(coerce_traces_to_coda(self$traces))
      stop_h("'traces_as_coda' is read only.")
    },

    #' @field traces_results (`data.frame`)\cr
    #' Results extracted from the traces of last inference.
    traces_results = function(value) {
      if (missing(value)) {
        return(retrieve_traces_results(self$traces_no_loglike))
      }
      stop_h("'traces_results' is read only.")
    },

    #' @field reference (`matrix()`)\cr
    #' Household loglikelihood traces.
    household_loglike_matrix = function(value) {
      if (missing(value)) {
        # check model has been built
        if (is.null(private$.rcpp)) {
          stop_h("Impossible to export household traces before inference.")
        }
        # extract household loglikelihood traces
        df <- private$.rcpp$export_household_loglike_matrix()
        # check
        if (dim(df)[1] == 0) {
          stop_h("Household traces are empty.")
        }
        return(df)
      }
      stop_h("'household_loglike_traces' is read only.")
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
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
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param quiet_pkgbuild (`logical(1)`)\cr
    #' If `TRUE` hide `pkgbuild` compilation commands.
    #' @param as_template (`logical(1)`)\cr
    #' If `TRUE` store the dynamic library as a local template.
    #' If `FALSE` recompile the model at every new session.
    initialize = function(model_choice,
                          force = FALSE,
                          quiet_pkgbuild = TRUE,
                          as_template = FALSE) {

      if (model_choice == ".") {
        # prohibit current directory
        stop_h(
          "Invalid folder name \".\" for housepi model.",
          "Move your implementation in a different folder."
        )
      }

      # assign quiet_pkgbuild field
      private$.quiet_pkgbuild <- isTRUE(quiet_pkgbuild)

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
              stop_h(
                "Impossibile to create a cache directory for compilation.",
                "Check your user permissions."
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
            private$.copy_source_code(model_choice)
          }
        }

      } else {
        #===== ERROR =====#
        stop_h(
          sprintf("\"%s\" is not a folder nor a template name.", model_choice),
          "Run function 'template_household_model' to see available templates."
        )
      }

      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Name of model dynamic library. Extract name of the dynamic library.
    #'
    #' @param verbose (`logical(1)`)\cr
    #' If `TRUE` show message about generation.
    generate = function(verbose = TRUE) {
      if (isTRUE(private$.need_compile)) {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Compiling {.pkg {private$.name}}")
        }
        # compile model
        private$.compile()
        private$.rcpp <- private$.new_rcpp_model()

      } else if (is.null(private$.rcpp)) {
        if (isTRUE(verbose)) {
          cli::cli_alert_info("Loading {.pkg {private$.name}}")
        }
        # load if not already generated
        private$.rcpp <- private$.new_rcpp_model()
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
    },

    ############################################################################
    ############################################################################
    ############################################################################

    #' @description
    #' Wrapper of the Rcpp import_dataset method.
    #'
    #' @param dataset (`data.frame`)\cr
    #' Dataset with rowise individual information.
    import_dataset = function(dataset) {
      private$.rcpp$import_dataset(
        data.frame(lapply(dataset, as.character))
      )
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp import_setup method.
    #'
    #' @param setup_list (`list()`)\cr
    #' Setup definitions.
    import_setup = function(setup_list) {
      # change scientific notation locally
      withr::local_options(list(scipen = 500))
      private$.rcpp$import_setup(
        lapply(setup_list, as.character)
      )
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp import_parameter method.
    #'
    #' @param defs (`data.frame`)\cr
    #' Formatted parameter definitions.
    import_parameter = function(defs) {
      private$.rcpp$import_parameter(
        data.frame(lapply(defs, as.character))
      )
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp import_augmented method.
    #'
    #' @param defs (`data.frame`)\cr
    #' Formatted augmented variable definitions.
    import_augmented = function(defs) {
      private$.rcpp$import_augmented(
        data.frame(lapply(defs, as.character))
      )
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp import_infec_profile method.
    #'
    #' @param defs (`list()`)\cr
    #' Infectivity profile definitions.
    import_infectivity_profile = function(defs) {
      private$.rcpp$import_infectivity_profile(
        lapply(defs, as.character)
      )
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp build_all_rates method.
    build_all_rates = function() {
      private$.rcpp$build_all_rates()
      invisible(self)
    },

    #' @description
    #' Export instantaneous and cumulative values of the infectivity profile.
    #'
    #' @param opts (`double()`)\cr
    #' Additional values for the computation.
    #' @param xlim (`double(2)`)\cr
    #' Customize the left and/or right limits for the x axis.
    #' To keep the default value of one limit, set it to `NA`.
    #' @param n_bins (`integer(1)`)\cr
    #' Number of bins for continuous profile (ignored if discrete).
    #'
    #' @return
    #' (`data.frame`)\cr
    export_infectivity_profile = function(opts = list(),
                                          xlim = double(),
                                          n_bins = 2000) {
      # load options
      opts_ <- as.list(as.numeric(opts))

      # load x-limits
      xlim_ <- list()
      xlim <- unlist(xlim)
      if (length(xlim) == 2) {
        if (!all(is.na(xlim))) {
          # check
          if (!any(is.na(xlim))) {
            if (xlim[[1]] == xlim[[2]]) {
              stop_h("Left x limit larger than right limit.")
            } else if (xlim[[1]] > xlim[[2]]) {
              stop_h("Left x limit equal to right limit.")
            }
          }
          xlim_ <- lapply(xlim, as.character)
          names(xlim_) <- c("min", "max")
        }

      } else if (length(xlim) > 0) {
        stop_h("Provide exactly 2 limits for the x axis.")
      }

      return(private$.rcpp$export_infectivity_profile(opts_, xlim_, n_bins))
    },

    #' @description
    #' Wrapper of the Rcpp evaluate method.
    evaluate = function() {
      return(private$.rcpp$evaluate())
    },

    #' @description
    #' Wrapper of the Rcpp infer method.
    infer = function() {
      private$.rcpp$infer()
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp infer method.
    simulate = function() {
      private$.rcpp$simulate()
      invisible(self)
    },

    #' @description
    #' Wrapper of the Rcpp estimate_evidence method.
    #'
    #' @param scheme (`character(1)`)\cr
    #' If "adaptive" uses adaptive scheme for thermodynamic ladder (default).
    #' Other available scheme is "power-fractions".
    #' @param n_rungs (`integer`)\cr
    #' Number of rungs in the thermodynamic ladder.
    estimate_evidence = function(scheme = "adaptive",
                                 n_rungs = 40) {
      return(private$.rcpp$estimate_evidence(scheme, n_rungs))
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
      HousepiGenerator$new(m,
                           force = force,
                           quiet_pkgbuild = FALSE,
                           as_template = TRUE)$generate()
    }
  }
}
