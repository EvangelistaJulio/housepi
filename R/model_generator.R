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

    # @field .cols_dataset `character()`\cr
    # Headers of the last dataset imported.
    .cols_dataset = NULL,

    # @field .cols_datalist `character()`\cr
    # Headers of the last datalist imported.
    .cols_datalist = NULL,

    # @field .levels_dataset `list()`\cr
    # List of levels of the factors in the dataset.
    .levels_dataset = NULL,

    # @field .levels_datalist `list()`\cr
    # List of levels of the factors in the datalist.
    .levels_datalist = NULL,

    # @field .cri_level `numeric(1)`\cr
    # Width of ther credible interval.
    .cri_level = 0.95,

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # When x$clone(deep=TRUE) is called, the deep_clone gets invoked once for
    # each field, with the name and value.
    deep_clone = function(name, value) {
      # make a new instance of rcpp model object
      if (name == ".rcpp") {
        if (is.null(value)) {
          output <- NULL
        } else {
          output <- private$.new_rcpp_model()
          output$import_setup(fmt_setup_input_(value$export_setup()))
          output$import_dataset(fmt_data_input_(value$export_dataset()))
          output$import_datalist(fmt_data_input_(value$export_datalist()))
          return(output)
        }
      }
      # for all other fields, just return the value
      return(value)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description
    # Assert source files exist.
    #
    # @param source_dir `character(1)`\cr
    # The path to source code
    .assert_source_files = function(source_dir) {
      for (file_name in source_files(source_dir)) {
        if (!file.exists(file_name)) stop_h(
          sprintf("File \"%s\" was not found.", file_name),
          "Check your model directory."
        )
      }
      invisible(self)
    },

    # @description
    # Assert cache allocations in source code match calls.
    #
    # @param source_dir `character(1)`\cr
    # The path to source code
    .assert_cache = function(source_dir) {
      # use source file as a reference
      refer_file <- file.path(source_dir, "model.cpp")

      # dictionary of allocate\call functions
      checks <- list(
        "cache_house_integer_vector" = ".cache_integer_vector",
        "cache_house_double_vector" = ".cache_double_vector",
        "cache_house_double_matrix" = ".cache_double_matrix"
      )
    
      required <- c()
      defined <- c()
      for (i in seq_along(checks)) {
        # regular expression
        reg <- '%s\\s*\\(\\s*["\'](.*?)["\']'
        # .extra_<type>("string")
        reg_called <- sprintf(reg, checks[[i]])
        # allocate_extra_house_<type>("string")
        reg_allocated <- sprintf(reg, names(checks)[[i]])

        # set of extra variables
        required <- c(required, template_required_find(refer_file, reg_called))
        defined <- c(defined, template_required_find(refer_file, reg_allocated))
      }
      required <- set_source_name_(unique(required), "model")
      # point out unused cached variables
      template_not_used_warn(required, defined, "cached variable")
      # stop if cached variables are missing
      template_missing_stop(required, defined, "cached variable")

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
    .new_custom_model = function(source_dir, as_template, force) {
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
        private$.assert_cache(source_dir)
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
      private$.save_reference()

      invisible(self)
    },

    # @description
    # Compile and load the source files inside cached directory.
    # Note: I am tricking pkgbuild::compile_dll and R CMD INSTALL.
    # To avoid recompiling immutable files I just touch them.
    .compile = function() {
      # drop previous instance
      private$.rcpp <- NULL

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
                            quiet = private$.quiet,
                            debug = FALSE)
      private$.need_compile <- FALSE

      # copy DLL outside of the source directory
      from_file <- file.path(cache_dir_src, dll_file_ext("housepienvironment"))
      overwrite(from_file, private$.dll)

      invisible(self)
    },

    # @description
    # Initialize new instance of Rcpp Model.
    .new_rcpp_model = function() {
      # load dynamic library
      dyn.load(private$.dll)
      # expose internally the module
      m <- Rcpp::Module(
        "housepi_model_core",
        PACKAGE = get_template_from_dll(private$.dll)
      )$ModelHousepi$new()
      m$init_advanced_setup()
      return(m)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description Save the reference script.
    .save_reference = function() {
      # reference script
      refer_file <- file.path(dirname(private$.dll), ".refer")
      overwrite(file.path(private$.cache_dir, "src", "model.cpp"),
                refer_file)

      invisible(self)
    },

    # @description Wheter model is blank.
    # @return `logical(1)`
    .is_blank = function() {
      return(private$.name == "blank")
    },

    # @description Envirnoment path.
    # @return `character(1)`
    .path = function() {
      # reference script
      output <- dirname(private$.dll)
      # check if reference file exists
      if (!dir.exists(output)) stop_h(
        sprintf("Impossible to find environment path: \"%s\".", output),
        "Try to re-generate the model with 'force = TRUE'."
      )
      return(output)
    },

    # @description Reference script.
    # @return `character(1)`
    .reference = function() {
      if (private$.is_blank()) {
        output <- file.path(private$.source_dir, "model.cpp")
      } else {
        # reference script
        output <- file.path(private$.path(), ".refer")
        # check if reference file exists
        if (!file.exists(output)) stop_h(
          sprintf("Impossible to find reference file: \"%s\".", output),
          "Try to re-generate the model with 'force = TRUE'."
        )
      }
      return(output)
    },

    # @description Reference script.
    # @param seed `integer(1)`
    # @return `integer(1)`
    .adjust_seed = function(seed) {
      max_r_int <- .Machine$integer.max
      if (is.na(seed)) seed <- sample.int(max_r_int, 1)
      else seed <- abs(seed)

      # control that seed is 32bit (for custom seed)
      if (seed > max_r_int) stop_h(
        "Supplied seed is not supported.",
        sprintf("Maximum value for an integer in R is %.0f.", max_r_int)
      )
      return(seed)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description 
    # Dataset headers required by the model.
    # @return `character()`
    .dataset_required = function() {
      # vector of regular expression to use for assertion
      regular_expr <- c(
        # from : .data("string")
        '.data\\s*\\(\\s*["\'](.*?)["\']'
      )
      # retrieve required dataset headers from regular expressions
      refer_file <- private$.reference()
      output <- template_required_find(refer_file, regular_expr)
      return(set_source_name_(output, "model"))
    },

    # @description 
    # Assert if dataset headers match implementation file.
    # @param dataset_headers `character()`\cr
    # Vector with the dataset headers.
    .assert_dataset = function(dataset_headers) {
      if (length(dataset_headers) == 0) stop_h(
        "Impossible to proceed. No database was provided."
      )
      # find missing dataset headers
      required <- private$.dataset_required()
      template_missing_stop(required, dataset_headers, "dataset")
      invisible(self)
    },

    ############################################################################
    ############################################################################
    ############################################################################

    # @description 
    # Assert if dataset headers match implementation file.
    # @param dataset_headers `character()`\cr
    # Vector with the dataset headers.
    .refactor_dataset = function(x) {
      cpp_data_refactor(x, private$.levels_dataset, private$.cols_dataset)
    },

    # @description 
    # Assert if dataset headers match implementation file.
    # @param dataset_headers `character()`\cr
    # Vector with the dataset headers.
    .refactor_datalist = function(x) {
      cpp_data_refactor(x, private$.levels_datalist, private$.cols_datalist)
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

    #' @field cri_level `double(1)`\cr
    #' Level for the credible interval.
    #' Value between 0.5 and 1 (default is 0.95).
    cri_level = function(value) {
      if (missing(value)) return(private$.cri_level)
      else if (is.na(value)) stop_h("Credible interval cannot be NA.")
      else if ((value <= 0.5) || (value >= 1)) stop_h(
        "Credible interval level must be a value between 0.5 and 1."
      )
      else private$.cri_level <- value
    },

    #' @field path_to_dll `character(1)`\cr
    #' Path to dynamic library file of the model.
    path_to_dll = function(value) {
      if (missing(value)) return(private$.dll)
      stop_h("'path_to_dll' is read only.")
    },

    #' @field setup `list()`\cr
    #' Last setup used by the build process.
    setup = function(value) {
      if (missing(value)) return(private$.rcpp$export_setup())
      stop_h("'setup' is read only.")
    },

    #' @field summary `data.frame`\cr
    #' A summary of statuses in the dataset.
    summary = function(value) {
      if (missing(value)) return(private$.rcpp$export_dataset_report())
      stop_h("'summary' is read only.")
    },

    #' @field summary_house_size `data.frame`\cr
    #' A summary of house sizes in the dataset.
    summary_house_size = function(value) {
      if (missing(value)) return(private$.rcpp$export_house_size_report())
      stop_h("'summary_house_size' is read only.")
    },

    #' @field traces `posterior::draws_df`\cr
    #' Posterior traces of last inference.
    traces = function(value) {
      if (missing(value)) {
        trs <- private$.rcpp$export_traces()
        return(cpp_traces_to_posterior(trs))
      }
      stop_h("'traces' is read only.")
    },

    #' @field house_loglike_traces `posterior::draws_df`\cr
    #' Traces of pointwise loglikelihood for all households.
    house_loglike_traces = function(value) {
      if (missing(value)) {
        trs <- private$.rcpp$export_house_loglike_traces()
        return(cpp_traces_to_posterior(trs))
      }
      stop_h("'house_loglike_traces' is read only.")
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
        warning_h(
          "No model specified. Falling back to \"blank\" model.",
          "To show available templates, run 'template_household_model'."
        )
        choice <- "blank"
      }
      else if (choice %in% list(".")) stop_h(
        sprintf("Invalid choice \"%s\" for housepi model.", choice),
        "Move your implementation in a different folder."
      )

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

      } else stop_h(
        sprintf("\"%s\" is not a folder nor a template name.", choice),
        "To show available templates, run 'template_household_model'."
      )

      # finally recompile only the model script, if necessary.
      if (private$.need_compile) private$.compile()
      # generate the Rcpp model object.
      if (is.null(private$.rcpp)) private$.rcpp <- private$.new_rcpp_model()

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info("Household model {.pkg {self$name}}.")
      invisible(self)
    },

    #' @description
    #' Show model description in RStudio Help window.
    help = function() {
      open_model_helper_(private$.path())
      invisible(self)
    },

    #' @description
    #' A snapshot of all current internal state of a given household.
    #' 
    #' @param house_id `character(1)`\cr
    #' Id of the household to export.
    #' 
    #' @return
    #' `list()`
    export_house_snapshot = function(house_id) {
      house_id <- as.character(house_id)
      output <- private$.rcpp$export_house_snapshot(house_id)
      output$dataset <- private$.refactor_dataset(output$dataset)
      output$datalist <- private$.refactor_datalist(output$datalist)

      # first draft to reformat infectivity profile output
      lst <- output$infec_profile
      df <- data.frame(id = lst$id, discr = lst$discr, name = lst$name,
                       stringsAsFactors = FALSE)
      df$args <- lst$args
      df$range <- lst$range
      output$infec_profile <- df

      return(output)
    }
  )
)
