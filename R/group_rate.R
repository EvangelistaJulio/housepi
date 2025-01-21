#
# group_rate.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

# formatting options
TAB <- "\t"
TAB2 <- strrep(TAB, 2)
REFERENCE <- "REFERENCE"

#' @title
#' Group rate logic table
#' @description
#' Return logic table of group rates
#' @param dataset
#' the dataset from which the etract group rates
#' @param by
#' Features as list of dataset column names
#' @return
#' data.frame with logic table
#' @noRd
group_rate_logic_table <- function(dataset, by) {
  # convert features
  features <- as.character(by)

  # assert that features are present
  missing <- c()
  for (name in features) {
    if (!(name %in% names(dataset))) {
      missing <- c(missing, name)
    }
  }
  if (length(missing) > 0) {
    oss <- sprintf("Unknown features(s): \"%s\"", missing[[1]])
    for (name in missing[-1]) {
      oss <- sprintf("%s, \"%s\"", oss, name)
    }
    stop(cli::format_error(c(
      paste0(oss, "."),
      i = "Check your features and dataset headers."
    )))
  }

  # extract unique features and sort by list length
  uni_features <- lapply(features, function(f) sort(unique(dataset[[f]])))
  sorted_indices <- order(sapply(uni_features, length), decreasing = TRUE)
  features <- features[sorted_indices]
  uni_features <- uni_features[sorted_indices]
  names(uni_features) <- features

  # create a dataframe of all possible multi-features
  logic_table <- expand.grid(uni_features)

  # check that the multi-feature exists
  multifeat_count <- apply(logic_table, 1, function(row) {
    sum(apply(dataset[features], 1, function(d_row) all(d_row == row)))
  })

  # create the final truth table
  logic_table <- logic_table[multifeat_count > 0, , drop = FALSE]

  # reverse table to ease transcript
  logic_table <- rev(logic_table)

  # add counts
  logic_table$..counts <- multifeat_count[multifeat_count > 0]

  return(logic_table)
}

#' @title 
#' Associate parameter names
#' @description
#' Associate parameter name to group rate logic table
#' @param logic_table
#' Data.frame containing conditions for the group rate.
#' @param prefix
#' Prefix to parameter names associated to the group rate.
#' @param reference
#' Choose reference providing corresponding features.
#' Default is the group with the largest sample size.
#' @return
#' input data.frame with logic table with parameter names.
#' @noRd
associate_param_names <- function(logic_table,
                                  prefix = NA,
                                  reference = NULL) {
  # extract reference index
  if (is.null(reference)) {
    idx_refer <- which.max(logic_table$..counts)
    # drop counts
    logic_table$..counts <- NULL

  } else {
    # drop counts
    logic_table$..counts <- NULL

    reference <- as.list(reference)[names(logic_table)]
    idx_refer <- unname(which(
      apply(logic_table, 1, function(row) all(row == unlist(reference)))
    ))
    if (length(idx_refer) != 1) {
      stop("Provided reference group is not valid.")
    }
  }

  # format prefix and associate paramter names
  if (is.na(prefix) || !is.character(prefix)) {
    prefix <- ""
  } else {
    prefix <- paste0(prefix, "_")
  }

  # associate parameter names
  logic_table$.param_name <- apply(
    logic_table, 1, function(row) {
      paste0(prefix, paste(paste0(names(logic_table), row), collapse = "_"))
    }
  )
  logic_table$.param_name[idx_refer] <- REFERENCE

  return(logic_table)
}

#' @description
#' construct the C++ if clauses
#' @param r
#' row of the group rate
#' @param is_ref
#' if `TRUE` the if clause is "REFERENCE"
#' @noRd
cpp_ifclause_line <- function(r,
                              is_ref = FALSE) {
  mask <- setdiff(names(r), ".param_name")
  # if condition
  condition <- paste(
    paste0(
      sprintf("house.data(i, \"%s\") == ", names(r[mask])),
      r[mask]), collapse = ") && ("
  )
  # wrap all together
  c(sprintf("%sif ((%s))", ifelse(is_ref, "", "else "), condition),
    sprintf("%semplace_group_rate(r, %s);",
            TAB,
            ifelse(is_ref, REFERENCE, sprintf("\"%s\"", r[[".param_name"]])))
  )
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#' @title Housepi Group Rate R6 Class
#' @description
#' Help with group rate definitions of model Rcpp object.
#' @export HousepiGroupRateClass
#' @noRd
HousepiGroupRateClass <- R6::R6Class(
  classname = "HousepiGroupRateClass",
  cloneable = FALSE,

  ###################
  ##### PUBLIC ######
  ###################

  public = list(
    defs = list(),
    param_defs = list(),
    logic_table = list(),

    #' @description
    #' Assign the empty parameter data.frame
    #' @noRd
    initialize = function() {
      self$drop_all()
      invisible(self)
    },

    #' @description
    #' Drop all current group rate definitions.
    #' @noRd
    drop_all = function() {
      self$defs <- list()
      self$param_defs = list()
      self$logic_table <- list()
      invisible(self)
    },

    #' @description
    #' Remove a given group_rate definition.
    #' @noRd
    remove = function(name,
                      show_warnings = TRUE) {

      if (any(names(self$defs) == name)) {
        # drop previous definition
        self$defs <- self$defs[names(self$defs) != name, ]

        # store corresponding logic talbe
        this_table <- self$logic_table[[name]]
        # drop previous logic table
        self$logic_table <- self$logic_table[names(self$logic_table) != name, ]

        # store list of necesary parameter definitions
        reserve_param <- c()
        for (rate_id in names(self$logic_table)) {
          reserve_param <- c(reserve_param,
                             unlist(self$logic_table[[rate_id]]$.param_name))
        }
        # drop corresponding parameters not used for other group rates
        for (par_id in this_table$.param_name) {
          if (!(par_id %in% reserve_param)) {
            definitions[[par_id]] <- self$param_defs[[rate_id]]
          }
        }

        return(TRUE)

      } else {
        if (isTRUE(show_warnings)) {
          cli::cli_alert_warning(cli::format_warning(c(
            sprintf("Unknown group_rate \"%s\" to remove.", name
          ))))
        }
      }

      return(FALSE)
    },

    #' @description
    #' Add a new group rate definition.
    #' @param name (`character(1)`)\cr
    #' Name of the group rate.
    #' @param by (`list()`)\cr
    #' Feature(s) associated to group rate as dataset headers.
    #' Corresponding dataset column have to be categorical.
    #' @param prefix (`character(1)`)\cr
    #' Prefix to parameter names associated to the group rate.
    #' @param reference (named `list()`)\cr
    #' Choose reference group providing corresponding features.
    #' Default is the group with the largest sample size.
    #' @param independent
    #' If `TRUE` features are considered separately.
    #' @param prior (`character(1)`)\cr
    #' Name of the prior distribution.
    #' @param args (named `list()`)\cr
    #' Arguments for the prior distribution.
    #' @param rule (`character(1)`)\cr
    #' Name of the inference proposal rule.
    #' @param rule_scale (`double(1)`)\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat`>0.
    #' @param min (`double(1)`)\cr
    #' Value of the minimum.
    #' @param max (`double(1)`)\cr
    #' Value of the maximum.
    #' @reference
    #' Cauchemez, et al. (2009)
    #' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
    #' New England Journal of Medicine. 361: 2619–2627. doi:10.1056/NEJMoa0905498
    #' @noRd
    add = function(name,
                   by,
                   prefix = NA,
                   reference = NULL,
                   independent = TRUE,
                   prior = "lnorm",
                   args = list(meanlog = 0, sdlog = 1),
                   rule = "lnorm",
                   rule_scale = 1,
                   min = NA,
                   max = 10) {

      # UPDATE: controls are missing
      if (is.na(name)) {
        stop(cli::format_error("Incomplete definition: missing name."))
      }

      # remove previous definition if it exists
      if (self$remove(name, FALSE)) {
        cli::cli_alert_warning(cli::format_warning(c(
          sprintf(
            "Existing \"%s\" definition has been overwritten.",
            name
        ))))
      } 

      # add new definition
      self$defs[[name]] <- list(
        by = by,
        prefix = prefix,
        reference = reference,
        independent = isTRUE(independent))

      # add parameter definitions
      self$param_defs[[name]] <- list(
        prior = prior,
        args = args,
        rule = rule,
        rule_scale = rule_scale,
        min = min,
        max = max
      )
    },

    #' @description
    #' Assert group rate definitions from the implementation script.
    #' @param refer_file (`character(1)`)\cr
    #' path/to/.refer containing the model source.
    #' @return (`logical(1)`)\cr
    #' @noRd
    assert = function(refer_file) {
      # logical state for complete definitions
      is_def_complete <- TRUE

      # vector of regular expression to use for assertion
      regular_expr <- c(
        # group rate from : .group_rate(..., "string"
        '\\.group_rate\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # retrieve required group rates from regular expressions
      require_set <- template_find_requirement(refer_file, regular_expr)

      if (is.null(require_set)) {
        # if requirement set is NULL it means nothing can be said
        is_def_complete <- NA

      } else {
        # find missing parameters
        missing_set <- template_find_missing(require_set,
                                        names(self$defs))
        if (length(missing_set) > 0) {
          oss <- paste0("Missing group rate(s): ",
                        character_to_string_list(missing_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))

          is_def_complete <- FALSE
        }
          
        # point out parameters that are not used
        unused_set <- template_find_unused(require_set, 
                                      names(self$defs))
        # drop unused parameters
        if (length(unused_set) > 0) {
          oss <- paste0("Dropping not used group rate(s): ",
                        character_to_string_list(unused_set), ".")
          cli::cli_alert_warning(cli::format_warning(oss))

          for (name in unused_set) {
            self$remove(name)
          }
        }
      }
      
      return(is_def_complete)
    },

    #' @description
    #' Generate a script to be printed within source code model.cpp
    #' @param dataset
    #' Dataset to be considered.
    #' @noRd
    cpp_implement = function(dataset) {

      # initialize lines of the script
      rate_script <- c("rateID r;")
      self$logic_table <- list()

      for (i in seq_along(self$defs)){
        rate_id <- names(self$defs[i])
        this_def <- self$defs[[i]]

        # default prefix
        prefix <- as.character(this_def[["prefix"]])
        if (is.na(prefix)) {
          prefix <- format_string(tolower(rate_id))
          prefix <- substring(prefix, 1, 3)
        }

        # features
        features <- as.list(this_def[["by"]])
        if (isFALSE(this_def[["independent"]])) {
          features <- list(features)
        }

        # initialize script
        script <- c()

        for (this_feat in features) {
          # create logic table with group rates
          this_table <- group_rate_logic_table(dataset, this_feat)
          # associate names for parameters
          this_table <- associate_param_names(this_table,
                                              prefix,
                                              this_def[["reference"]])
          # mask reference case
          mask_refer <- this_table$.param_name == REFERENCE
          script <- c(
            script,
            # start with reference case
            cpp_ifclause_line(this_table[mask_refer, ], TRUE),
            # other non-reference cases
            apply(this_table[!mask_refer, ], 1, cpp_ifclause_line),
            # close with missing case
            "else",
            sprintf("%semplace_group_rate(r, MISSING_RATE);", TAB)
          )
        }

        # append the new rate chunk of script
        rate_script <- c(rate_script, 
                         "",
                         sprintf("//----- \"%s\" -----//", rate_id),
                         sprintf("r = \"%s\";", rate_id),
                         script)
        # store the new rate logic table
        self$logic_table[[rate_id]] <- this_table
      }

      # return the script collapsed
      return(paste(c(
          "/* Generated automatically -> do not edit by hand */",
          "void Model::assign_group_rate_to_individual(Household &house, indiv i)",
          "{",
          sapply(rate_script, function(x) { paste0(TAB, x) }),
          "}"
        ),
        collapse = "\n"))
    },

    #' @description
    #' Assert if group rates are implementated.
    #' @param dataset
    #' Dataset to be considered.
    #' @param source_file
    #' Path/to/model.cpp source code.
    #' @return (`logical(1)`)\cr
    #' `TRUE` if definitions exist and are up to date. 
    #' @noRd
    cpp_assert = function(dataset,
                          source_file) {
      # check if there are user defined group rates
      if (length(self$defs) == 0) {
        # default is TRUE (i.e. no definitions)
        return(TRUE)

      } else {
        # assert if group rates are complete
        if (isFALSE(self$assert(source_file))) {
          cli::cli_alert_warning(cli::format_warning(c(
            "Ignoring incomplete group rate definitions.",
            i = "To implement changes in group rates, fill missing definitions."
          )))
          return(TRUE)

        } else {
          # produce associated script
          rate_script <- self$cpp_implement(dataset)
          # read source file
          lines <- readLines(source_file)
          # check if the definitions are already implemented
          return(stringr::str_detect(paste(lines, collapse = "\n"),
                                      stringr::fixed(rate_script)))
        }
      }
    },

    #' @description
    #' Overwrite script to source code model.cpp
    #' @param dataset
    #' Dataset to be considered.
    #' @return (`logical(1)`)\cr
    #' Wheter code requires recompilation or not.
    #' @noRd
    cpp_write = function(dataset,
                         source_file) {

      # check that source file exists
      if (!file.exists(source_file)) {
        stop(cli::format_error("Impossible to find source file."))
      }

      # produce associated script
      rate_script <- self$cpp_implement(dataset)
      # read source file
      lines <- readLines(source_file)

      # find the start of the method definition
      start_line <- grep("^void\\s+Model::assign_group_rate_to_individual\\s*\\(", lines)

      if (length(start_line) == 0) {
        cli::cli_alert_warning("Model does not accommodate for group rates.")

      } else {

        # locate method definition opening brace '{'
        opening_brace_line <- NA
        for (i in start_line:length(lines)) {
          if (grepl("\\{", lines[i])) {
            opening_brace_line <- i
            break
          }
        }
        if (is.na(opening_brace_line)) {
          stop("Unmatched curly bracket(s) found in source code.")
        }

        # find method definition closing brace '}'
        brace_level <- 0
        end_line <- NA
        for (i in opening_brace_line:length(lines)) {
          # update brace level based on opening and closing braces
          brace_level <- brace_level +
            stringr::str_count(lines[i], "\\{") -
            stringr::str_count(lines[i], "\\}")
          # closing brace corresponds to level 0
          if (brace_level == 0) {
            end_line <- i
            break
          }
        }
        if (is.na(end_line)) {
          stop("Unmatched curly bracket(s) found in source code.")
        }

        # define the new implementation
        lines <- c(
          lines[1:(start_line - 1)],
          rate_script,
          lines[(end_line + 1):length(lines)]
        )

        # replace the source code
        writeLines(lines, source_file)
      }

      invisible(self)
    }
  )
)
