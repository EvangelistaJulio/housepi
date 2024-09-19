
# Define a custom get function (like in Python dictionaries)
get_val_ <- function(lst, key, default = NA) {
  return(if (key %in% names(lst)) lst[[key]] else default)
}

# Format NA types with defaults according to values
default_types_ <- function(typs, vals) {
  for (i in which(is.na(typs))) {
    y <- vals[[i]]
    
    if (is.numeric(y)) {
      typs[i] <- "S"
      
    } else if (is.character(y)) {
      if(is.na(suppressWarnings(as.numeric(y)))) {
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
  vals <- sapply(strsplit(as.character(lst), "\\@"), function(x) x[1])
  types <- sapply(strsplit(as.character(lst), "\\@"), function(x) x[2])

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
        add(name = p, x = get_val_(this, "x"),
            infer = get_val_(this, "infer"),
            prior = get_val_(this, "prior"),
            args = get_val_(this, "args", default = list()),
            rule = get_val_(this, "rule"),
            rule_args = get_val_(this, "rule_args", default = list()),
            range = get_val_(this, "range", default = list()))
      }
      
    },
    
    #===== Method : add =====#
    add = function(name = NA, x = NA, infer = NA,
                   prior = NA, args = list(),
                   rule = NA, rule_args = list(),
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
        cli::cli_alert_warning(cli::format_warning(
          c(paste0("Existing `", name, "` definition has been overwritten."))))
        df <<- df[df$name != name, ]
      }
      
      if (is.na(infer)) {
        infer <- !is.na(rule)
        if (is.na(rule) && (!is.na(prior))) {
          cli::cli_alert_warning(cli::format_warning(
            c(paste0("`", name, "` choice interpreted as if infer = FALSE."),
            i = "Provide one between rule and infer to silence this warning.")))
        }
      }
      
      if (infer == TRUE) {
        #----- prior -----#
        if (is.na(prior)){
          stop(cli::format_error(
            c(paste0("Provide a prior name for `", name, "`."))))
        }
        
        if (length(args) == 0) {
          args_ <- NA
          argsval_ <- NA
        } else {
          args_ <- str_join_(names(args))
          argsval_ <- str_join_(unname(unlist(args)))
        }
        
        #----- rule -----#
        if (is.na(rule)){
          stop(cli::format_error(
            c(paste0("Provide a rule name for `", name, "`."))))
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
        
        df <<- rbind(df, data.frame(
          name = name, value = x, infer = TRUE,
          prior = prior, args = args_, argsval = argsval_,
          rule = rule, rargs = rargs_, rargsval = rargsval_,
          min = min_, max = max_,
          stringsAsFactors = FALSE)
        )
        
      } else {
        df <<- rbind(df, data.frame(
          name = name, value = x, infer = FALSE,
          prior = NA, args = NA, argsval = NA,
          rule = NA, rargs = NA, rargsval = NA,
          min = NA, max = NA,
          stringsAsFactors = FALSE)
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
        infer = logical(),
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
        add(name = a, x = get_val_(this, "x"), 
            infer = get_val_(this, "infer"), 
            model = get_val_(this, "model"), 
            args = get_val_(this, "args", default = list()),
            rule = get_val_(this, "rule"), 
            rule_args = get_val_(this, "rule_args", default = list()),
            range = get_val_(this, "range", default = list()))
      }
      
    },
    
    #===== Method : add =====#
    add = function(name = NA, x = NA, infer = NA, 
                   model = NA, args = list(), 
                   rule = NA, rule_args = list(), 
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
      
      #----- value -----#
      value_ <- str_split_(x)[["values"]]
      valuetype_ <- str_split_(x)[["types"]]
      
      if (name %in% df[["name"]])
      {
        cli::cli_alert_warning(cli::format_warning(
          c(paste0("Existing `", name, "` definition has been overwritten."))))
        df <<- df[df$name != name, ]
      }
      
      if (is.na(infer)) {
        infer <- !is.na(rule)
        if (is.na(rule)) {
          cli::cli_alert_warning(cli::format_warning(
            c(paste0("`", name, "` choice interpreted as if infer = FALSE."),
              i = "Provide a rule to silence this warning.")))
        }
      }
      
      if(infer == TRUE) {

        if(length(args) == 0)
        {
          args_ <- NA
          argsval_ <- NA
          argstype_ <- NA
        } else {
          args_ <- str_join_(names(args))
          tmp <- str_split_(args)
          argsval_ <- str_join_(tmp[["values"]])
          argstype_ <- str_join_(tmp[["types"]])
        }
        
        #----- rule -----#
        if(is.na(rule)){
          stop(cli::format_error(
            c(paste0("Provide a rule name for `", name, "`."))))
          
        } else if(length(rule_args) == 0) {
          rargs_ <- NA
          rargsval_ <- NA
          rargstype_ <- NA
          
        } else {
          rargs_ <- str_join_(names(rule_args))
          tmp <- str_split_(rule_args)
          rargsval_ <- str_join_(tmp[["values"]])
          rargstype_ <- str_join_(tmp[["types"]])
        }
        
        #----- range -----#
        if (length(range) == 0) {
          min_ <- NA
          mintype_ <- NA
          max_ <- NA
          maxtype_ <- NA
          
        } else {
          
          tmp <- str_split_(args)
          rargsval_ <- str_join_(tmp[["values"]])
          rargstype_ <- str_join_(tmp[["types"]])
          
          if("min" %in% names(range)) {
            tmp <- str_split_(range[names(range) == "min"])
            min_ <- tmp[["values"]]
            mintype_ <- tmp[["types"]]
          } else {
            min_ <- NA
            mintype_ <- NA 
          }
          
          if("max" %in% names(range)) {
            aux <- str_split_(range[names(range) == "max"])
            max_ <- tmp[["values"]]
            maxtype_ <- tmp[["types"]]
          } else {
            max_ <- NA
            maxtype_ <- NA 
          }
        }
        
        df <<- rbind(df, data.frame(
          name = name, 
          value = value_, valuetype = valuetype_,
          infer = TRUE,
          model = model, 
          args = args_, argsval = argsval_, argstype = argstype_,
          rule = rule, 
          rargs = rargs_, rargsval = rargsval_, rargstype = rargstype_,
          min = min_, mintype = mintype_,
          max = max_, maxtype = maxtype_,
          stringsAsFactors = FALSE)
        )
        
      } else {
        df <<- rbind(df, data.frame(
          name = name, value = x, infer = FALSE,
          model = NA, args = NA, argsval = NA,
          rule = NA, rargs = NA, rargsval = NA,
          min = NA, max = NA,
          stringsAsFactors = FALSE)
        )
      }
    }
  )
)
