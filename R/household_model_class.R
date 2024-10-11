#
# model_class.R
#
# Copyright (C) 2024 Francesco Camaglia, Institut Pasteur, Paris (France).
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
  vals <- sapply(strsplit(as.character(lst), "\\@"), function(x)
    x[1])
  types <- sapply(strsplit(as.character(lst), "\\@"), function(x)
    x[2])
  
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

#' @title Household Model Class
#' @description 
#' Manage your compiled household model, your data and definitions from a single class.
#' @export HouseholdModelClass
#' @docType class
HouseholdModelClass <- R6::R6Class(
  "HouseholdModelClass",
  private = list(
    setup = NULL,
    param.frame = NULL,
    augm.frame = NULL,
    obj = NULL
  ),

  public = list(

    #' @description
    #' Creates a new instance of class [ModelClass][housepi::HouseholdModelClass].
    #' @param model (`character(1)`)\cr 
    #' Template name or path to model folder.
    #' @param quiet (`logical(1)`)\cr 
    #' Hide compilation commands.
    #' @param force (`logical(1)`)\cr
    #' Compile even if a cached model already exists.
    #' @param debug (`logical(1)`)\cr
    #' Use debug compilation flags.
    initialize = function(model = NA,
                          quiet = TRUE,
                          force = FALSE,
                          debug = FALSE) {

      if (is.na(model) ) {
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
          # use cached version
          cli::cli_alert_info("Using cached {.pkg {pkgload::pkg_name(exec_dir)}}")
          
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
                c("Failed to setup internal files.", i = "Check your user permissions.")
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
                c("Failed to setup internal files.", i = "Check your user permissions.")
              ))
            }
          }
          
          # copy source files in the new folder
          for (from_file in source_files) {
            to_file <- file.path(exec_dir, "src", basename(from_file))
            if (!file.copy(from_file, to_file, overwrite = TRUE)) {
              stop(cli::format_error(
                c("Failed to setup implementation files.", i = "Check your user permissions.")
              ))
            }
          }
          
          # ----- compile package ----- #
          # FIXME : how to enhance pre-compilation?
          pkgbuild::compile_dll(exec_dir, quiet = quiet, debug = debug)
        }
        
        # load shared library
        dyn.load(executable)
        
        # expose internally the module
        m <- Rcpp::Module("model_housepi_rcpp", PACKAGE = "housepienvironment")
        
        #----- assign fields -----#
        
        # assign the model object
        private$obj <- m$ModelHousepi$new()
        
        # assign the empty parameter data.frame
        self$param.clear()
        
        # assign the empty augmented data.frame
        self$augm.clear()
      }

      invisible(self)
    },
  
    #' @description 
    #' Copy a default model to a given directory.  
    #' @param which - String with default model name.
    #' @param path - String with path to where to copy template scripts.
    template = function(which = NULL, path = NULL) {
  
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
            c("Provided model is not a default model name.", i = "To see the available default models, run: template_model()")
          ))
        }
        
        # copy files
        for (f in c("model.cpp", "model.h")) {
          from_file <- file.path(models_dir, which, f)
          to_file <- file.path(new_dir, f)
          if (!file.copy(from_file, to_file, overwrite = TRUE)) {
            stop(cli::format_error(
              c("Failed to setup internal files.", i = "Check your user permissions.")
            ))
          }
        }
      } else {
        stop(cli::format_error(
          c("Provided model is not a default model name.", i = "To see the available default models, run: template_model()")
        ))
      }
      invisible(self)
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
  
    #' @description
    #' Associate data and definitions to the model.
    #' @param data - Dataframe or String with path to data file.
    #' @param verbose - Boolean flag for automatic prints.
    build = function(data = NULL, verbose = FALSE) {
      
      # change scientific notation
      store_scipen <- options()$scipen
      options(scipen = 500)
      
      #----- 1) SETUP -----#
      private$obj$import_setup(lapply(private$setup, as.character))
      
      #----- 2) PARAMETERS -----#
      private$obj$import_params(data.frame(lapply(private$param.frame, as.character)))
      
      #----- 3) DATA -----#
      if (is.null(data)) {
        private$obj$read_data("<NA>")
        
      } else if (is.character(data)) {
        private$obj$read_data(data)
        
      } else if (is.data.frame(data)) {
        private$obj$import_data(data.frame(lapply(data, as.character)))
        
      } else {
        stop(cli::format_error(
          c("Wrong type for data", i = "It should be data.frame or path/to/file.")
        ))
      }
      
      #----- 4) AUGMENTED VARIABLES -----#
      private$obj$import_augm(data.frame(lapply(private$augm.frame, as.character)))
      
      options(scipen = store_scipen)
      invisible(self)
    },
    
    #' @description
    #' Run the inference.
    #' @param progbar - Boolean flag to suppress progress bar.
    #' @return A coda mcmc object.
    infer = function(progbar = TRUE) {
      # call infer method from cpp model object
      private$obj$infer(progbar)      
      # load traces from last inference
      df <- private$obj$traces()
      
      # extract setup information for coda
      # FIXME : this should be the internal one
      # FIXME : to be safe, build should be run automatically before infer?
      S <- private$obj$setup()
      
      # Check if multiple chains
      if ("CHAIN__" %in% colnames(df)) {
        df_list <- lapply(split(df, df[["CHAIN__"]]), function(x) {
          coda::mcmc(
            data = x[, !names(x) %in% "CHAIN__"],
            start = S[["burn_in"]] + S[["thinning"]],
            end = S[["num_iter"]],
            thin = S[["thinning"]]
          )
        })
        output <- coda::as.mcmc.list(df_list)
        
      } else {
        output <- coda::mcmc(
          data = df,
          start = S[["burn_in"]] + S[["thinning"]],
          end = S[["num_iter"]],
          thin = S[["thinning"]]
        )
      }
      
      return(output)
      invisible(self)
    },
  
    #' @description
    #' Run the simulation.  
    #' @param progbar - Boolean flag to suppress progress bar.
    simulate = function(progbar = TRUE) {
  
      private$obj$simulate(progbar)
      invisible(self)
    },
    
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    
    #' @description
    #' Remove all current definitions.
    param.clear = function() {
      private$param.frame <- data.frame(
        name = character(),
        value = numeric(),
        infer = logical(),
        prior = character(),
        args = character(),
        argsval = character(),
        rule = character(),
        rargs = character(),
        rargsval = character(),
        min = numeric(),
        max = numeric(),
        stringsAsFactors = FALSE
      )
      invisible(self)
    },
  
    #' @description
    #' Append a data frame of parameter definitions.
    #' @param definitions - List of parameter names with respective list of definitions.
    param.append = function(definitions = list()) {
      
      if (length(definitions) == 0) {
        # reinitialize according to templated parameters
        self$param.template(private$obj)
        
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
        self$param.check(private$obj)
      }
      invisible(self)
    },
    
    #' @description
    #' Overwrite current parmater definitions with template definitions.
    param.template = function() {
      # retrieve the template parameters file
      hemi_dir <- system.file(package = attr(class(private$obj), "package"))
      templ_params_file <- file.path(hemi_dir, "src", "templated_params")
      # remove previous definitions
      self$param.clear()
      # append definitions
      self$param.append(dget(templ_params_file))
      invisible(self)
    },
    
    #' @description
    #' Check if parameter definitions are missing and pre-load them.
    param.check = function() {
      
      hemi_dir <- system.file(package = attr(class(private$obj), "package"))
      model_path <- file.path(hemi_dir, "src", "model.cpp")
      
      # read the script and extract parameters
      script <- readChar(model_path, file.info(model_path)$size)
      
      # from : get_param("string") or get_param('string')
      regular_expr <- 'get_param\\(\\s*["\'](.*?)["\']\\s*\\)'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- unique(unlist(matches[[1]][, 2]))
      
      # from : get_param_pointer("string") or get_param_pointer('string')
      regular_expr <- 'get_param_pointer\\(\\s*["\'](.*?)["\']\\s*\\)'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- c(model_params_set, unique(unlist(matches[[1]][, 2])))
      
      # from : assign_rate_per_contact(, , "string")
      regular_expr <- '\\.assign_rate_per_contact\\s*\\(.*?,.*?,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      aux <- unlist(lapply(matches, function(x)
        x[, 2]))
      model_params_set <- c(model_params_set, unique(aux[!is.na(aux)]))
      
      # from : assign_rate_per_group(, "string")
      regular_expr <- '\\.assign_group_rate\\s*\\(.*?,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      aux <- unlist(lapply(matches, function(x)
        x[, 2]))
      model_params_set <- c(model_params_set, unique(aux[!is.na(aux)]))
      
      # re-order alphabetically
      model_params_set <- model_params_set[order(model_params_set)]
      
      # add missing parameters with standard options
      for (name in model_params_set) {
        if (!(name %in% private$param.frame[["name"]])) {
          self$param.add(
            name = name,
            x = 1,
            infer = FALSE,
            prior = "unif",
            args = list(min = 0, max = 10),
            rule = "norm",
            rule_args = list(sd = 1)
          )
        }
      }
      # point out parameters that are not used
      unused <- c()
      for (name in private$param.frame[["name"]]) {
        if (!(name %in% model_params_set)) {
          unused <- c(unused, name)
        }
      }
      if (length(unused) > 0) {
        oss <- paste0("Parameter(s) not declared within model: '", unused[[1]], "'")
        for (name in unused[-1]) {
          oss <- paste0(oss, ", '", name, "'")
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }
      invisible(self)
    },
    
    #' @description
    #' Add a new parameter definition.
    #' @param name - String with name of the parameter.
    #' @param x - Scalar with value of the parameter.
    #' @param infer - Logical value indicating whether to infer.
    #' @param prior - String with the prior distribution name.
    #' @param args - List of arguments for the prior distribution.
    #' @param rule - String with the inference rule.
    #' @param rule_args - List of arguments for the inference rule.
    #' @param range - List with 'min' and/or 'max' values (or none).
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
          rargs_ <- NA
          rargsval_ <- NA
        } else {
          rargs_ <- str_join_(names(rule_args))
          rargsval_ <- str_join_(unname(unlist(rule_args)))
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
            rargs = rargs_,
            rargsval = rargsval_,
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
            rargs = NA,
            rargsval = NA,
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
    #' Remove all current definitions.
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
        rargs = character(),
        rargsval = character(),
        rargstype = character(),
        rfrom = character(),
        rfromtype = character(),
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
    #' @param definitions - List of augmented variable names with respective list of definitions.
    augm.append = function(definitions = list()) { 
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
  
    #' @description Add a new augmented variable to definitions.
    #' @param name - String with name of the parameter.
    #' @param x - Scalar/string with value of the parameter.
    #' @param model - String with the model distribution name.
    #' @param args - List of arguments for the prior distribution.
    #' @param rule - String with the inference rule.
    #' @param rule_args - List of arguments for the inference rule.
    #' @param rule_from - Scalar/string with addtional rule options.
    #' @param range - List with 'min' and/or 'max' values (or none).
    #' @param init_rand - Logical to randomly initialize with model.
    #' @param status - Integer of the associated status.
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
        if (name %in% private$augm.frame[["name"]])
        {
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
            c("Incomplete definition: missing value x.", i = "Alternatively provide a model and set init_rand = TRUE")
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
        if (length(args) == 0)
        {
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
      rargs_ <- NA
      rargsval_ <- NA
      rargstype_ <- NA
      rfrom_ <- NA
      rfromtype_ <- NA
      if (!is.na(rule)) {
        if ((rule != "gibbs") && (length(rule_args) == 0)) {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Missing rule_args for augmented variable '", name, "'.")
          )))
        } else {
          rargs_ <- str_join_(names(rule_args))
          tmp <- str_split_(rule_args)
          rargsval_ <- str_join_(tmp[["values"]])
          rargstype_ <- str_join_(tmp[["types"]])
          
          rfrom_ <- str_split_(rule_from)[["values"]]
          rfromtype_ <- str_split_(rule_from)[["types"]]
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
          rargs = rargs_,
          rargsval = rargsval_,
          rargstype = rargstype_,
          rfrom = rfrom_,
          rfromtype = rfromtype_,
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