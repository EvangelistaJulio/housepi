
# Define a custom get function (like in Python dictionaries)
get_val_ <- function(lst, key, default = NA) {
  return(if (key %in% names(lst))
    lst[[key]]
    else
      default)
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

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' @title Parameter definitions
#' @description A class for managing a data frame of parameter definitions.
#' @field df A data frame that stores parameter definitions.
#' @examples
#' # example for the parameters of a gamma distribution
#' params <- ParameterFrame(list(
#'   # define the `shape` parameter
#'   shape = list(
#'     x = 1,                          # initial value
#'     prior = "uniform",              # choose a uniform prior
#'     args = list(min = 0,            # minimum of the uniform prior
#'                 max = 5),           # maximum of the uniform prior
#'     rule = "normal",                # Metropolis normal move
#'     rule_args = list(sd = 0.5)      # initial sd of the move (later adapted)
#'   ),
#'   # define the `scale` parameter
#'   scale = list(
#'     x = 0.5,                        # initial value
#'     infer = FALSE,                  # do to infer scale
#'     prior = "lognorm",              # choose a lognormal prior
#'     args = list(meanlog = 0,        # meanlog of lognormal prior
#'                 sdlog = 1),         # sdlog of lognormal prior
#'     rule = "prior"                  # sample from the prior
#'  )
#' ))
#'
#' # to add or overwrite a definition, e. g. a model constant
#' params$add(name = "rule_sd", x = 0.125, infer = FALSE)
#'
#' # to visualize the definitions data frame
#' params$df
#' @export ParameterFrame
#' @exportClass ParameterFrame
#' @docType class
ParameterFrame <- setRefClass(
  #===== NAME =====#
  "ParameterFrame",
  
  #===== ATTRIBUTES =====#
  fields = list(df = "data.frame"),
  
  #===== METHODS =====#
  methods = list(
    #===== Method : initialize/new =====#
    initialize = function(params = list()) {
      df <<- data.frame(
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
      
      for (p in names(params)) {
        this <- params[[p]]
        add(
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
      
    },
    
    #===== from_model =====#
    from_model = function(housepi_model) {
      "Pre-load missing parameter definitions needed by the model.
      \\itemize{
      \\item Parameter: housepi_model - housepi model object.
      }"
      
      hemi_dir <- system.file(package = attr(class(housepi_model), "package"))
      model_path <- file.path(hemi_dir, "src", "model.cpp")
      
      # read the script
      script <- readChar(model_path, file.info(model_path)$size)
      
      # extract from get_param("string") or get_param('string')
      regular_expr <- 'get_param\\(\\s*["\'](.*?)["\']\\s*\\)'
      matches <- stringr::str_match_all(script, regular_expr)
      model_params_set <- unique(unlist(matches[[1]][, 2]))
      
      # extract from : assign_rate_per_contact(, , "string")
      regular_expr <- '\\.assign_rate_per_contact\\s*\\(.*?,.*?,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      aux <- unlist(lapply(matches, function(x)
        x[, 2]))
      model_params_set <- c(model_params_set, unique(aux[!is.na(aux)]))
      
      # extract from : assign_rate_per_group(, "string")
      regular_expr <- '\\.assign_rate_per_group\\s*\\(.*?,\\s*["\'](.*?)["\']'
      matches <- stringr::str_match_all(script, regular_expr)
      aux <- unlist(lapply(matches, function(x)
        x[, 2]))
      model_params_set <- c(model_params_set, unique(aux[!is.na(aux)]))
      
      # re-order alphabetically
      model_params_set <- model_params_set[order(model_params_set)]
      
      # add missing parameters with standard options
      for (name in model_params_set) {
        if (!(name %in% df[["name"]])) {
          add(
            name = name,
            x = 1,
            infer = FALSE,
            prior = "uniform",
            args = list(min = 0, max = 10),
            rule = "normal",
            rule_args = list(sd = 1)
          )
        }
      }
      
      # point out parameters that are not used
      unused <- c()
      for (name in df[["name"]]) {
        if (!(name %in% model_params_set)) {
          unused <- c(unused, name)
        }
      }
      if (length(unused) > 0) {
        oss <- paste0("The following parameter(s) are not present in the model: `",
                      unused[[1]],
                      "`")
        for (name in unused[-1]) {
          oss <- paste0(oss, ", `", name, "`")
        }
        cli::cli_alert_warning(cli::format_warning(c(paste0(oss, "."))))
      }
    },
    
    #===== Method : add =====#
    add = function(name = NA,
                   x = NA,
                   infer = NA,
                   prior = NA,
                   args = list(),
                   rule = NA,
                   rule_args = list(),
                   range = list()) {
      "Add a new parameter to data frame
      \\itemize{
      \\item Parameter: name - String with name of the parameter.
      \\item Parameter: x - Scalar with value of the parameter.
      \\item Parameter: infer - Logical value indicating whether to infer.
      \\item Parameter: prior - String with the prior distribution name.
      \\item Parameter: args - List of arguments for the prior distribution.
      \\item Parameter: rule - String with the inference rule.
      \\item Parameter: rule_args - List of arguments for the inference rule.
      \\item Parameter: range - List with `min` and/or `max` values (or none).
      }"
      
      if (is.na(name) || is.na(x)) {
        stop(cli::format_error(c(
          "Incomplete definition: provide both name and value x."
        )))
      }
      
      if (name %in% df[["name"]]) {
        cli::cli_alert_warning(cli::format_warning(c(
          paste0("Overwriting existing `", name, "` definition.")
        )))
        df <<- df[df$name != name, ]
      }
      
      if (is.na(infer)) {
        infer <- !is.na(rule)
        if (is.na(rule) && (!is.na(prior))) {
          cli::cli_alert_warning(cli::format_warning(
            c(
              paste0("`", name, "` choice interpreted as if infer = FALSE."),
              i = "Provide one between rule and infer to silence this warning."
            )
          ))
        }
      }
      
      if (infer == TRUE) {
        #----- prior -----#
        if (is.na(prior)) {
          stop(cli::format_error(c(
            paste0("Provide a prior name for `", name, "`.")
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
            paste0("Provide a rule name for `", name, "`.")
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
        
        df <<- rbind(
          df,
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
        df <<- rbind(
          df,
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
    }
  )
)

#' @title Augmented variable definitions
#' @description A class for managing a data frame of augmented variable definitions.
#' @field df A data frame that stores augmented variable definitions.
#' @examples
#' # associate to each individual an exposure and an onset time
#' # the definitions relies on quantities that will be imported from data
#' augm <- AugmentedFrame(list(
#'   # define augmented variable 'exposure' for exposure time
#'   exposure = list(
#'     x = "exposure_min@D", # start from exposure_min imported from data
#'     model = "uniform",    # associate a uniform model
#'     args = list(min = "exposure_min@D",  # 'exposure_min' imported from data
#'                 max = "exposure_max@D"), # 'exposure_max' imported from data
#'     rule = "normal",      # rule for Metropolis move is normal
#'     rule_args = list(sd = 0.125)         # std. dev. of normal rule
#'   ),
#'   # define augmented variable 'onset' for onset time
#'   onset = list(
#'     x = "onset_max@D",    # start from maximum # ...
#'     model = "uniform",
#'     args = list(min = "onset_min@D",
#'                 max = "onset_max@D"),
#'     rule = "normal",
#'     rule_args = list(sd = "rule_sd@P")   # 'rule_sd' imported from parameters
#'   )
#' ))
#'
#' # to visualize the definitions data frame
#' augm$df
#' @export AugmentedFrame
#' @exportClass AugmentedFrame
#' @docType class
AugmentedFrame <- setRefClass(
  #===== NAME =====#
  "AugmentedFrame",
  
  #===== ATTRIBUTES =====#
  fields = list(df = "data.frame"),
  
  #===== METHODS =====#
  methods = list(
    #===== Method : initialize/new =====#
    initialize = function(augm = list()) {
      df <<- data.frame(
        name = character(),
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
      
      for (a in names(augm)) {
        this <- augm[[a]]
        add(
          name = a,
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
      
    },
    
    #===== Method : add =====#
    add = function(name = NA,
                   x = NA,
                   model = NA,
                   args = list(),
                   rule = NA,
                   rule_args = list(),
                   rule_from = NA,
                   range = list(),
                   init_rand = FALSE) {
      "Add a new parameter to data frame
      \\itemize{
      \\item Parameter: name - String with name of the parameter.
      \\item Parameter: x - Scalar/string with value of the parameter.
      \\item Parameter: prior - String with the prior distribution name.
      \\item Parameter: args - List of arguments for the prior distribution.
      \\item Parameter: rule - String with the inference rule.
      \\item Parameter: rule_args - List of arguments for the inference rule.
      \\item Parameter: rule_from - Scalar/string with addtional rule options.
      \\item Parameter: range - List with `min` and/or `max` values (or none).
      \\item Parameter: init_rand - Logical to randomly initialize with model.
      }"
      
      if (is.na(name) || is.na(x)) {
        stop(cli::format_error(c(
          "Incomplete definition: provide both name and value x."
        )))
      }
      
      #----- value -----#
      value_ <- str_split_(x)[["values"]]
      valuetype_ <- str_split_(x)[["types"]]
      
      if (name %in% df[["name"]])
      {
        cli::cli_alert_warning(cli::format_warning(c(
          paste0("Existing `", name, "` definition has been overwritten.")
        )))
        df <<- df[df$name != name, ]
      }
      
      # assign method only if arguments are provided
      args_ <- NA
      argsval_ <- NA
      argstype_ <- NA
      
      if (!is.na(model))  {
        if (length(args) == 0)
        {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Missing args for augmented variable `", name, "`.")
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
        if (length(rule_args) == 0) {
          cli::cli_alert_warning(cli::format_warning(c(
            paste0("Missing rule_args for augmented variable `", name, "`.")
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
        min_ <- NA
        mintype_ <- NA
        max_ <- NA
        maxtype_ <- NA
        
      } else {
        
        if ("min" %in% names(range)) {
          tmp <- str_split_(range[names(range) == "min"])
          min_ <- tmp[["values"]]
          mintype_ <- tmp[["types"]]
        } else {
          min_ <- NA
          mintype_ <- NA
        }
        
        if ("max" %in% names(range)) {
          aux <- str_split_(range[names(range) == "max"])
          max_ <- tmp[["values"]]
          maxtype_ <- tmp[["types"]]
        } else {
          max_ <- NA
          maxtype_ <- NA
        }
      }
      
      df <<- rbind(
        df,
        data.frame(
          name = name,
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
    }
  )
)
