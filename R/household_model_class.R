#
# household_model_class.R
#
# Copyright (C) 2024 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# Define a custom get function (like in Python dictionaries)
get_val_ <- function(lst, key, default = NA) {
  if (key %in% names(lst)) {
    return(lst[[key]])
  } else {
    return(default)
  }
}

# Format NA types with defaults according to values
default_types_ <- function(typs, vals) {
  for (i in which(is.na(typs))) {
    y <- vals[[i]]
    if (is.na(y)) {
      typs[i] <- NA
    } else if (is.numeric(y)) {
      typs[i] <- "S"
    } else if (is.character(y)) {
      if (is.na(suppressWarnings(as.numeric(y)))) {
        typs[i] <- "D"
      } else {
        typs[i] <- "S"
      }
    } else {
      typs[i] <- "D"
    }
  }
  return(typs)
}

# Split a typed string or character vector
str_split_ <- function(lst) {
  vals <- sapply(
    strsplit(as.character(lst), "\\@"), function(x) x[1]
  )
  types <- sapply(
    strsplit(as.character(lst), "\\@"), function(x) x[2]
  )
  return(list(values = vals, types = default_types_(types, vals)))
}

# Join elements of a dictionary to pass to csv
str_join_ <- function(x) {
  return(paste(x, collapse = "/"))
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

#' @title household model template
#' @description
#' Copy a default model to a given directory.
#' @param which (`character(1)`)\cr
#' Name of template model. If NULL, list of available templates.
#' @param path (`character(1)`)\cr
#' Where to copy the template scripts.
#' @export
household_model_template <- function(which = NULL, path = NULL) {
  # system path to housepi
  models_dir <- file.path(system.file(package = "housepi"), "models")
  
  if (is.null(which)) {
    cli::cli_alert_info("Default model names:")
    ulid <- cli::cli_ul()
    for (name in list.files(models_dir)) {
      cli::cli_li(name)
    }
    cli::cli_end(ulid)
    
  } else if (is.character(which)) {
    if (which %in% list.files(models_dir)) {
      # create directory
      if (is.null(path)) {
        new_dir <- file.path(getwd(), paste0("user_", which))
      } else {
        new_dir <- path
      }
      dir.create(new_dir, recursive = TRUE, showWarnings = TRUE)
      if (!dir.exists(new_dir)) {
        stop(cli::format_error(
          c(
            "Impossibile to create a temporary directory for the model.",
            i = "Check your user permissions."
          )
        ))
      }
      
    } else {
      stop(cli::format_error(
        c(
          "Provided model is not a default model name.",
          i = "To see the available default models, run: template_model()"
        )
      ))
    }
    
    #----- copy files ----#
    for (f in c("model.cpp", "model.h", "templated_params")) {
      from_file <- file.path(models_dir, which, f)
      to_file <- file.path(new_dir, f)
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c(
            "Failed to setup internal files.",
            i = "Check your user permissions."
          )
        ))
      }
    }
  } else {
    stop(cli::format_error(
      c("Provided model is not a default model name.",
        i = "To see the available default models, run: template_model()")
    ))
  }
}

#' @title Household Model Class
#' @description
#' Compile household model, manage your data and parameters from a single class.
#' @export HouseholdModelClass
#' @docType class
HouseholdModelClass <- R6::R6Class(
  classname = "HouseholdModelClass",
  cloneable = FALSE,
  private = list(
    setup = NULL,
    param.frame = NULL,
    augm.frame = NULL,
    traces.frame = NULL,
    obj = NULL,
    environment_path = NULL,
    
    #' @description
    #' Associate data and definitions to the model.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param verbose (`logical(1)`)\cr
    #' If `TRUE` prints automatic build information.
    build = function(data = NULL, verbose = FALSE) {
      # change scientific notation
      store_scipen <- options()$scipen
      options(scipen = 500)
      
      #----- 1) SETUP -----#
      char_setup <- lapply(private$setup[!is.na(private$setup)], as.character)
      private$obj$import_setup(char_setup)
      #----- 2) PARAMETERS -----#
      char_params <- data.frame(lapply(private$param.frame, as.character))
      private$obj$import_params(char_params)
      #----- 3) DATA -----#
      if (is.null(data)) {
        private$obj$read_data("<NA>")
      } else if (is.character(data)) {
        private$obj$read_data(data)
      } else if (is.data.frame(data)) {
        char_data <- data.frame(lapply(data, as.character))
        private$obj$import_data(char_data)
      } else {
        stop(cli::format_error(
          c(
            "Wrong type for data", 
            i = "It should be data.frame or path/to/file."
          )
        ))
      }
      #----- 4) AUGMENTED VARIABLES -----#
      char_augm <- data.frame(lapply(private$augm.frame, as.character))
      private$obj$import_augm(char_augm)
      
      # restore scientific notation option
      options(scipen = store_scipen)
    }
  ),
  
  ############################################################################
  ############################################################################
  ############################################################################
  ############################################################################
  ############################################################################
  
  public = list(
    #' @description
    #' Creates a new instance of class [HouseholdModelClass][housepi::HouseholdModelClass].
    #' @param model (`character(1)`)\cr
    #' Template name or path to the model folder.
    #' @param quiet (`logical(1)`)\cr
    #' If `TRUE` hide compilation commands.
    #' @param force (`logical(1)`)\cr
    #' If `TRUE` compile even if a cached model already exists.
    #' @param debug (`logical(1)`)\cr
    #' If `TRUE` use debug compilation flags.
    initialize = function(model = NA,
                          quiet = TRUE,
                          force = FALSE,
                          debug = FALSE) {
      if (is.na(model)) {
        # print available templates
        self$template()
        
      } else {
        # system path to housepi
        housepi_dir <- system.file(package = "housepi")
        
        if (dir.exists(model)) {
          #----- custom user model -----#
          source_files <- file.path(model, c("model.cpp", "model.h"))
          # check that header and implementation files exist
          for (file_name in source_files) {
            if (!file.exists(file_name)) {
              stop(cli::format_error(c(
                paste0("A file '", file_name, "' was not found."),
                i = paste0("Check your model folder.")
              )))
            }
          }
          
        } else if (model %in% list.files(file.path(housepi_dir, "models"))) {
          #----- template model -----#
          source_files <- file.path(housepi_dir,
                                    "models",
                                    model,
                                    c("model.cpp", "model.h", "templated_params"))
          
        } else {
          #----- error -----#
          stop(cli::format_error(paste0(
            "'", model, "' is not an folder nor a template name."
          )))
        }
        
        # temporary folder
        env_name <- paste0("housepi", gsub("_", "", model))
        # FIXME :
        # env_name <- sprintf("housepi%x", as.integer(file.info(cpp_file)$mtime))
        # temporary session directory where to store model environment
        exec_dir <- file.path(tempdir(), env_name)
        
        # ----- Create model package ----- #
        
        # Returns the path/to/DLL file supposedly generated in this session
        executable <- file.path(exec_dir,
                                "src",
                                paste0("housepienvironment", .Platform$dynlib.ext))
        
        # option to force the removal of the package directory
        if ((file.exists(executable) && (force == TRUE))) {
          unlink(exec_dir, recursive = TRUE)
        }
        
        # deal with existing version
        if (file.exists(executable)) {
          # use cached version # TODO: to be tuned with phgbuild compile_dll
          cli::cli_alert_info(
            "Using cached {.pkg {pkgload::pkg_name(exec_dir)}}"
          )
        } else {
          # copy the package to the new directory
          dir.create(exec_dir,
                     recursive = TRUE,
                     showWarnings = FALSE)
          if (!dir.exists(exec_dir)) {
            stop(cli::format_error(
              c(
                "Impossibile to create a temporary directory for the model.",
                i = "Check your user permissions."
              )
            ))
          }
          # copy single package files
          for (f in c("NAMESPACE", "DESCRIPTION")) {
            from_file <- file.path(housepi_dir, "template", f)
            to_file <- file.path(exec_dir, f)
            if (!file.copy(from_file, to_file, overwrite = TRUE)) {
              stop(cli::format_error(
                c(
                  "Failed to setup internal files.",
                  i = "Check your user permissions."
                )
              ))
            }
          }
          # copy folders from template
          for (fldr in c("src", "R")) {
            from_path <- file.path(housepi_dir, "template", fldr)
            if (!file.copy(from_path,
                           exec_dir,
                           recursive = TRUE,
                           overwrite = TRUE)) {
              stop(cli::format_error(
                c(
                  "Failed to setup internal files.",
                  i = "Check your user permissions."
                )
              ))
            }
          }
          # copy source files in the new folder
          for (from_file in source_files) {
            to_file <- file.path(exec_dir, "src", basename(from_file))
            if (!file.copy(from_file, to_file, overwrite = TRUE)) {
              stop(cli::format_error(
                c(
                  "Failed to setup implementation files.", 
                  i = "Check your user permissions."
                )
              ))
            }
          }
          
          # ----- compile package ----- #
          # TODO : how to enhance pre-compilation?
          pkgbuild::compile_dll(exec_dir, quiet = quiet, debug = debug)
        }
        
        # load shared library
        dyn.load(executable)
        # expose internally the module
        m <- Rcpp::Module("model_housepi_rcpp", PACKAGE = "housepienvironment")
        
        #----- assign fields -----#
        # assign the path to model environment source
        private$environment_path <- exec_dir
        # assign the model object
        private$obj <- m$ModelHousepi$new()
        # assign the empty setup list
        self$setup.clear()
        # assign the empty parameter data.frame
        self$param.clear()
        # assign the empty augmented data.frame
        self$augm.clear()
      }
      
      invisible(self)
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #' @description
    #' Run Bayesian MCMC inference and store internally the chain traces.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param progbar (`logical(1)`)\cr
    #' If `FALSE` suppress the inference progress bar.
    #' @param ... further arguments passed to \code{setup.list}.
    infer = function(data = NULL, progbar = TRUE, ...) {
      # assign setup
      do.call(self$setup.list, list(...))
      # internal build function
      private$build(data = data)
      # infer method from cpp model object
      private$obj$infer(progbar)
      # load traces from last inference
      df <- private$obj$export_traces()
      # extract setup information for coda
      load_setup <- private$obj$export_setup()
      # Check if multiple chains
      if ("CHAIN__" %in% colnames(df)) {
        df_list <- lapply(split(df, df[["CHAIN__"]]), function(x) {
          coda::mcmc(
            data = x[, !names(x) %in% "CHAIN__"],
            start = load_setup[["burn_in"]] + load_setup[["thinning"]],
            end = load_setup[["num_iter"]],
            thin = load_setup[["thinning"]]
          )
        })
        private$traces.frame <- coda::as.mcmc.list(df_list)
      } else {
        private$traces.frame <- coda::mcmc(
          data = df,
          start = load_setup[["burn_in"]] + load_setup[["thinning"]],
          end = load_setup[["num_iter"]],
          thin = load_setup[["thinning"]]
        )
      }
      invisible(self)
    },
    
    #' @description
    #' Export the traces of last inference call as a coda object.
    #' @return
    #' (`coda::mcmc` | `coda::mcmc.list`)
    traces.export = function() {
      return(private$traces.frame)
    },

    #' @description
    #' Export the correlations between traces values of last inference call.
    #' @return
    #' A (`matrix` | 3d `array`)
    traces.crosscorr = function() {
      return(coda::crosscorr(private$traces.frame))
    },

    #' @description
    #' Export the traces of last inference call as a coda object.
    #' @return
    #' (`data.frame`)
    traces.summary = function() {
      res <- summary(private$traces.frame)
      imask <- private$param.frame["infer"] == TRUE
      inames <- private$param.frame[imask, "name"]
      output <- data.frame(
        Mean = res[["statistics"]][inames, "Mean"],
        Std.Dev = res[["statistics"]][inames, "SD"],
        Q2.5 = res[["quantiles"]][inames, "2.5%"],
        Q97.5 = res[["quantiles"]][inames, "97.5%"]
      )
      return(output)
    },
    
    #' @description
    #' Calls `coda` traceplot function to plot traces of a single chain.
    #' @param chain (`integer(1)`)
    #' Option to specifiy which chain to plot if multiple chains.
    #' @param col (`character(1)`)
    #' Name of the color or RGB hex value.
    #' @param ...
    #' Further arguments pass to [coda][coda::traceplot].
    traces.traceplot = function(chain = NA,
                                col = "dodgerblue",
                                ...) {
      if (is.null(private$traces.frame)) {
        stop("No traces found.")
      } else if (coda::is.mcmc.list(private$traces.frame)) {
        if (is.na(chain)) {
          coda::traceplot(private$traces.frame, col = col, ...)
        } else {
          coda::traceplot(private$traces.frame[chain], col = col, ...)
        }
      } else {
        coda::traceplot(private$traces.frame, col = col, ...)
      }
      
      invisible(self)
    },
    
    #' @description
    #' Run the simulation.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param progbar (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @param ... further arguments passed to \code{setup.list}.
    #' Named build options.
    simulate = function(data = NULL, progbar = TRUE, ...) {
      # assign setup
      do.call(self$setup.list, list(...))
      # internal build function
      private$build(data = data)
      # simulate method from internal Rcpp object
      private$obj$simulate(progbar)
      invisible(self)
    },
    
    #' @description
    #' Export the loaded database at current state.
    #' If the case, it also return the associated augmented variables.
    #' @return (`data.frame`)\cr
    #' A data.frame with individual data.
    data.export = function() {
      return(private$obj$export_data())
    },
    
    #' @description
    #' Print a summary of the loaded database.
    data.summary = function() {
      return(private$obj$show_summary())
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #' @description
    #' Use loo package to compute the leave-one-out cross-validation (PSIS-LOO CV).
    #' @param ... arguments passed to [loo][loo::loo].
    #' For arguments and further details type ?loo::loo.
    loo = function(...) {
      return(loo::loo(private$obj$export_loglike_matrix(), ...))
    },
    
    #' @description
    #' Use loo package to compute the widely applicable information criterion (WAIC).
    #' @param ... arguments passed to [loo][loo::waic].
    #' #' For arguments and further details type ?loo::waic.
    waic = function(...) {
      return(loo::waic(private$obj$export_loglike_matrix(), ...))
    },
    
    #' @description
    #' Estimate the model evidence (or marginal likelihood) with the method of power posteriors.
    #' @param data (`character(1)` | `data.frame`)\cr
    #' Path to a csv data file or a data.frame.
    #' @param num_rungs (`integer`)\cr
    #' The number of rungs in the thermodynamic ladder.
    #' @param adaptive (`logical(1)`)\cr
    #' If `TRUE` uses adaptive scheme for thermodynamic ladder. Else use power fractions.
    #' @param progbar (`logical(1)`)\cr
    #' If `FALSE` suppress the simulation progress bar.
    #' @return (`double(1)`)\cr
    #' The value of the model evidence.
    #' @param ... further arguments passed to  \code{setup.list}.
    evidence = function(data = NULL,
                        num_rungs = 40,
                        adaptive = TRUE,
                        progbar = TRUE,
                        ...) {
      # assign setup
      do.call(self$setup.list, list(...))
      # internal build function
      private$build(data = data)
      # evaluate marginal from the internal Rcpp object
      return(private$obj$evaluate_marginal(num_rungs, adaptive, progbar))
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #' @description
    #' Clear all current setup definitions.
    setup.clear = function() {
      private$setup <- list()
      invisible(self)
    },
    
    #' @description
    #' Setup the list of configurations for the different analysis.
    #' @param seed (`integer(1)`)\cr
    #' Seed of the random number generator.
    #' @param num_iter (`integer(1)`)\cr
    #' Number of iterations.
    #' @param burn_in (`integer(1)`)\cr
    #' Number of iterations to discard.
    #' @param thinning (`integer(1)`)\cr
    #' Name of the parameter.
    #' @param num_chains (`integer(1)`)\cr
    #' Name of the parameter.
    #' @param num_adapt (`integer(1)`)\cr
    #' Name of the parameter.
    #' @param num_threads (`integer(1)`)\cr
    #' Name of the parameter.
    #' @param rnd_init_params (`logical(1)`)\cr
    #' Name of the parameter.
    #' @param inv_temper (`double(1)`)\cr
    #' Name of the parameter.
    setup.list = function(seed = NA,
                          num_iter = 1000,
                          burn_in = 0,
                          thinning = 1,
                          num_chains = 1,
                          num_adapt = 0,
                          num_threads = 1,
                          rnd_init_params = TRUE,
                          inv_temper = 1.0) {
      private$setup[["seed"]] <- seed
      private$setup[["num_iter"]] <- num_iter
      private$setup[["burn_in"]] <- burn_in
      private$setup[["thinning"]] <- thinning
      private$setup[["num_chains"]] <- num_chains
      private$setup[["num_adapt"]] <- num_adapt
      private$setup[["num_threads"]] <- num_threads
      private$setup[["rnd_init_params"]] <- rnd_init_params
      private$setup[["inv_temper"]] <- inv_temper
      invisible(self)
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #' @description
    #' Export the parameter definitions at current state.
    #' @return (`data.frame`)\cr
    #' A data.frame with formatted parameter definitions.
    param.export = function() {
      return(private$param.frame)
    },
    
    #' @description
    #' Clear all current parameter definitions.
    param.clear = function() {
      private$param.frame <- data.frame(
        name = character(),
        value = numeric(),
        infer = logical(),
        prior = character(),
        args = character(),
        argsval = character(),
        rule = character(),
        rule_args = character(),
        rule_argsval = character(),
        min = numeric(),
        max = numeric(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },
    
    #' @description
    #' Append a list of parameter definitions.
    #' @param ... (named `list()`s)\cr
    #' Parameter names with respective list of definitions.
    param.append = function(...) {
      # convert definitions to list
      definitions <- list(...)
      
      if (length(definitions) == 0) {
        # reinitialize according to templated parameters
        self$param.template()
        cli::cli_alert_warning(cli::format_warning(
          c(
            "No definitions provided.", 
            i = "Falling back to template definitions."
          )
        ))
        
      } else {
        # user definitions
        for (p in names(definitions)) {
          this <- definitions[[p]]
          self$param.add(
            name = p,
            x = get_val_(this, "x"),
            infer = get_val_(this, "infer"),
            prior = get_val_(this, "prior"),
            args = get_val_(this, "args", default = list()),
            rule = get_val_(this, "rule"),
            rule_args = get_val_(this, "rule_args", default = list()),
            range = get_val_(this, "range", default = list())
          )
        }
        
        # check if definitions match the model
        self$param.check()
      }
      invisible(self)
    },
    
    #' @description
    #' Overwrite current parameter definitions with template definitions.
    param.template = function() {
      # retrieve the template parameters file
      templ_params_file <- file.path(private$environment_path, "src", "templated_params")
      if (!file.exists((templ_params_file))) {
        stop(cli::format_error("Impossible to automatically initialize parameters. No template parameters were found."))
      }
      # remove previous definitions
      self$param.clear()
      # append definitions
      do.call(self$param.append, dget(templ_params_file))
      invisible(self)
    },
    
    #' @description
    #' Check if parameter definitions are missing and pre-load them.
    param.check = function() {
      # retrieve the model script
      model_path <- file.path(private$environment_path, "src", "model.cpp")
      
      # read the script and extract parameters
      script <- readChar(model_path, file.info(model_path)$size)
      
      # from : get_param("string") or get_param('string')
      regular_expr <- 'get_param\\(\\s*["\'](.*?)["\']\\s*\\)'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- unique(unlist(matches[[1]][, 2]))
      
      # param from : get_param_pointer("string") or get_param_pointer('string')
      regular_expr <- 'get_param_pointer\\(\\s*["\'](.*?)["\']\\s*\\)'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- c(model_params_set, unique(unlist(matches[[1]][, 2])))
      
      # param from : assign_rate_per_contact(, , "string")
      regular_expr <- '\\.assign_rate_per_contact\\s*\\(.*?,.*?,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      aux <- unlist(lapply(matches, function(x) x[, 2]))
      model_params_set <- c(model_params_set, unique(aux[!is.na(aux)]))
      
      # param from : assign_rate_per_group(, "string")
      regular_expr <- '\\.assign_group_rate\\s*\\(.*?,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      aux <- unlist(lapply(matches, function(x) x[, 2]))
      model_params_set <- c(model_params_set, unique(aux[!is.na(aux)]))
      
      # re-order alphabetically
      model_params_set <- model_params_set[order(model_params_set)]
      
      # add missing parameters with standard options
      # FIXME : they should be added as templates maybe
      miss_params <- c()
      for (name in model_params_set) {
        if (!(name %in% private$param.frame[["name"]])) {
          miss_params <- c(miss_params, name)
          self$param.add(
            name = name,
            x = 1,
            infer = FALSE,
            prior = "unif",
            args = list(min = -500, max = 500),
            rule = "norm",
            rule_args = list(sd = 1)
          )
        }
      }
      if (length(miss_params) > 0) {
        oss <- paste0("Added missing parameter(s): '", miss_params[[1]], "'")
        for (name in miss_params[-1]) {
          oss <- paste0(oss, ", '", name, "'")
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }

      # point out parameters that are not used
      unused <- c()
      for (name in private$param.frame[["name"]]) {
        if (!(name %in% model_params_set)) {
          unused <- c(unused, name)
        }
      }
      if (length(unused) > 0) {
        oss <- paste0("Parameter(s) not used in model: '", unused[[1]], "'")
        for (name in unused[-1]) {
          oss <- paste0(oss, ", '", name, "'")
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }
      invisible(self)
    },
    
    #' @description
    #' Add a new parameter definition.
    #' @param name (`character(1)`)\cr
    #' Name of the parameter.
    #' @param x (`double(1)`)\cr
    #' Initial value of the parameter.
    #' @param infer (`logical(1)`)\cr
    #' If `FALSE` the parameter is treated as a constant.
    #' @param prior (`character(1)`)\cr
    #' Name of the prior distribution.
    #' @param args (named `list()` of `double(1)`)\cr
    #' Arguments for the prior distribution.
    #' @param rule (`character(1)`)\cr
    #' Name of the inference rule.
    #' @param rule_args (named `list()` of `double(1)`)\cr
    #' Arguments for the inference rule. They can be changed automatically during iterations with adaptation.
    #' @param range (named `list()` of `double(1)`)\cr
    #' The prior truncation 'min', 'max' or no bound for the parameter value.
    param.add = function(name = NA,
                         x = NA,
                         infer = NA,
                         prior = NA,
                         args = list(),
                         rule = NA,
                         rule_args = list(),
                         range = list()) {
      if (is.na(name)) {
        stop(cli::format_error(c("Missing parameter name.")))
      }
      
      if (name %in% private$param.frame[["name"]]) {
        cli::cli_alert_warning(cli::format_warning(c(
          paste0("Overwriting existing '", name, "' definition.")
        )))
        private$param.frame <- private$param.frame[private$param.frame$name != name, ]
      }
      
      #----- extrapolate infer value -----#
      if (is.na(infer)) {
        infer <- !is.na(rule)
        if (is.na(rule) && (!is.na(prior))) {
          cli::cli_alert_warning(cli::format_warning(
            c(paste0(
              "'", name, "' interpreted as if infer = FALSE."
            ), i = "Provide one between rule/infer to silence this warning.")
          ))
        }
      }
      
      if (infer == TRUE) {
        #----- prior -----#
        if (is.na(prior)) {
          stop(cli::format_error(c(
            paste0("Provide a prior name for '", name, "'.")
          )))
        }
        
        if (length(args) == 0) {
          args_ <- NA
          argsval_ <- NA
        } else {
          args_ <- str_join_(names(args))
          argsval_ <- str_join_(unname(unlist(args)))
        }
        
        #----- rule -----#
        if (is.na(rule)) {
          stop(cli::format_error(c(
            paste0("Provide a rule name for '", name, "'.")
          )))
        } else if (length(rule_args) == 0) {
          rule_args_ <- NA
          rule_argsval_ <- NA
        } else {
          rule_args_ <- str_join_(names(rule_args))
          rule_argsval_ <- str_join_(unname(unlist(rule_args)))
        }
        
        #----- range -----#
        min_ <- get_val_(names(range), "min")
        max_ <- get_val_(names(range), "max")
        
        #----- private$param.frame : append -----#
        private$param.frame <- rbind(
          private$param.frame,
          data.frame(
            name = name,
            value = x,
            infer = TRUE,
            prior = prior,
            args = args_,
            argsval = argsval_,
            rule = rule,
            rule_args = rule_args_,
            rule_argsval = rule_argsval_,
            min = min_,
            max = max_,
            stringsAsFactors = FALSE
          )
        )
        
      } else {
        if (is.na(x)) {
          stop(cli::format_error(
            paste0("Missing value for constant parameter '", name, "'.")
          ))
        }
        private$param.frame <- rbind(
          private$param.frame,
          data.frame(
            name = name,
            value = x,
            infer = FALSE,
            prior = NA,
            args = NA,
            argsval = NA,
            rule = NA,
            rule_args = NA,
            rule_argsval = NA,
            min = NA,
            max = NA,
            stringsAsFactors = FALSE
          )
        )
      }
      invisible(self)
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #' @description
    #' Export the augmented variable definitions at current state.
    #' @return (`data.frame`)\cr
    #' A data.frame with formatted augmented variable definitions.
    augm.export = function() {
      return(augm.frame)
    },
    
    #' @description
    #' Remove all current augmented variable definitions.
    augm.clear = function() {
      private$augm.frame <- data.frame(
        name = character(),
        status = character(),
        value = character(),
        valuetype = character(),
        model = character(),
        args = character(),
        argsval = character(),
        argstype = character(),
        rule = character(),
        rule_args = character(),
        rule_argsval = character(),
        rule_argstype = character(),
        rule_from = character(),
        rule_fromtype = character(),
        min = character(),
        mintype = character(),
        max = character(),
        maxtype = character(),
        init_rand = logical(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },
    
    #' @description
    #' Append a data frame of augmented variable definitions.
    #' @param ... (named `list()`s)\cr
    #' Augmented variable names with respective list of definitions.
    augm.append = function(...) {
      # convert definitions to list
      definitions <- list(...)
      # add each entry of the list
      for (a in names(definitions)) {
        this <- definitions[[a]]
        self$augm.add(
          name = a,
          status = get_val_(this, "status"),
          x = get_val_(this, "x"),
          model = get_val_(this, "model"),
          args = get_val_(this, "args", default = list()),
          rule = get_val_(this, "rule"),
          rule_args = get_val_(this, "rule_args", default = list()),
          rule_from = get_val_(this, "rule_from"),
          range = get_val_(this, "range", default = list()),
          init_rand = get_val_(this, "init_rand", default = FALSE)
        )
      }
      invisible(self)
    },
    
    #' @description Add a new augmented variable definition.
    #' @param name (`character(1)` )\cr
    #' Name of the augmented variable.
    #' @param x (`double(1)` | `character(1)`)\cr
    #' Scalar value or parameter name (followed by `@P`) or data name (followed by `@D`).
    #' @param model (`character(1)` )\cr
    #' Name of the model distribution.
    #' @param args - List of arguments for the prior distribution.
    #' @param rule (`character(1)` )\cr
    #' Name of the inference rule.
    #' @param rule_args (named `list()` of `double(1)` | `character(1)`)\cr
    #' Arguments for the inference rule.
    #' @param rule_from (`double(1)` | `character(1)`)\cr
    #' Advanced options for non random walk moves.
    #' @param range - List with 'min' and/or 'max' values (or none).
    #' @param init_rand (`character(1)` )\cr
    #' If `TRUE` draw initial vale from model distribution.
    #' @param status (`integer(1)` )\cr
    #' The associated status of the definition. Default NA is for all infectious.
    augm.add = function(name = NA,
                        x = NA,
                        model = NA,
                        args = list(),
                        rule = NA,
                        rule_args = list(),
                        rule_from = NA,
                        range = list(),
                        init_rand = FALSE,
                        status = NA) {
      #----- name -----#
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))
      } else {
        if (name %in% private$augm.frame[["name"]]) {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Existing '", name, "' definition has been overwritten.")
          )))
          private$augm.frame <- private$augm.frame[private$augm.frame$name != name, ]
        }
      }
      
      #----- value/init_rand -----#
      if (init_rand) {
        if (is.na(model)) {
          stop(cli::format_error("Unknown model for random initialization."))
        } else {
          value_ <- NA
          valuetype_ <- NA
        }
      } else {
        if (is.na(x)) {
          stop(cli::format_error(
            c(
              "Incomplete definition: missing value x.",
              i = "Alternatively provide a model and set init_rand = TRUE"
            )
          ))
        } else {
          value_ <- str_split_(x)[["values"]]
          valuetype_ <- str_split_(x)[["types"]]
        }
      }
      
      #----- method -----#
      args_ <- NA
      argsval_ <- NA
      argstype_ <- NA
      # assign method only if arguments are provided
      if (!is.na(model))  {
        if (length(args) == 0) {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Missing args for augmented variable '", name, "'.")
          )))
        } else {
          args_ <- str_join_(names(args))
          tmp <- str_split_(args)
          argsval_ <- str_join_(tmp[["values"]])
          argstype_ <- str_join_(tmp[["types"]])
        }
      }
      
      #----- rule -----#
      rule_args_ <- NA
      rule_argsval_ <- NA
      rule_argstype_ <- NA
      rule_from_ <- NA
      rule_fromtype_ <- NA
      if (!is.na(rule)) {
        if ((rule != "gibbs") && (length(rule_args) == 0)) {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Missing rule_args for augmented variable '", name, "'.")
          )))
        } else {
          rule_args_ <- str_join_(names(rule_args))
          tmp <- str_split_(rule_args)
          rule_argsval_ <- str_join_(tmp[["values"]])
          rule_argstype_ <- str_join_(tmp[["types"]])
          rule_from_ <- str_split_(rule_from)[["values"]]
          rule_fromtype_ <- str_split_(rule_from)[["types"]]
        }
      }
      
      #----- range -----#
      if (length(range) == 0) {
        # empty range
        min_ <- NA
        mintype_ <- NA
        max_ <- NA
        maxtype_ <- NA
      } else {
        # minimum value
        if ("min" %in% names(range)) {
          tmp <- str_split_(range[names(range) == "min"])
          min_ <- tmp[["values"]]
          mintype_ <- tmp[["types"]]
        } else {
          min_ <- NA
          mintype_ <- NA
        }
        # maximum value
        if ("max" %in% names(range)) {
          tmp <- str_split_(range[names(range) == "max"])
          max_ <- tmp[["values"]]
          maxtype_ <- tmp[["types"]]
        } else {
          max_ <- NA
          maxtype_ <- NA
        }
      }

      #----- append the new definitions -----#
      private$augm.frame <- rbind(
        private$augm.frame,
        data.frame(
          name = name,
          status = status,
          value = value_,
          valuetype = valuetype_,
          model = model,
          args = args_,
          argsval = argsval_,
          argstype = argstype_,
          rule = rule,
          rule_args = rule_args_,
          rule_argsval = rule_argsval_,
          rule_argstype = rule_argstype_,
          rule_from = rule_from_,
          rule_fromtype = rule_fromtype_,
          min = min_,
          mintype = mintype_,
          max = max_,
          maxtype = maxtype_,
          init_rand = init_rand,
          stringsAsFactors = FALSE
        )
      )
      invisible(self)
    }
  )
)