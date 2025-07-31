#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# group_rates.R
#
#' @title
#' Assert dataset features
#' @description
#' Assert that features are present in dataset,
#'
#' @param features
#' Features as list of dataset column names
#' @param dataset_headers
#' Headers of the dataset.
#'
#' @keywords internal
assert_dataset_features <- function(features, dataset_headers) {
  # assert that features are present
  missing <- c()
  for (this in features) {
    if (!(this %in% dataset_headers)) missing <- c(missing, this)
  }
  if (length(missing) > 0) stop_h(
    sprintf(
      "Unknown feature%s %s",
      ifelse(length(missing) > 1, "s:", ""),
      character_to_string_list(missing)
    ),
    "Check your features and dataset headers."
  )
}

#' @title
#' Group rate logic table
#' @description
#' Return logic table of group rates
#'
#' @param dataset
#' the dataset from which the etract group rates
#' @param by
#' Features as list of dataset column names
#'
#' @return
#' data.frame with logic table
#'
#' @keywords internal
group_rate_logic_table <- function(dataset, by) {
  # sort alphabetically the features
  features <- sort(as.character(by))

  # convert features
  assert_dataset_features(features, colnames(dataset))

  # extract unique features and sort by list length
  uni_features <- lapply(features, function(f) sort(unique(dataset[[f]])))
  sorted_indices <- order(sapply(uni_features, length), decreasing = TRUE)
  features <- features[sorted_indices]
  uni_features <- uni_features[sorted_indices]
  names(uni_features) <- features

  # create a dataframe of all possible multi-features
  y <- expand.grid(uni_features)

  # check that the multi-feature exists
  multifeat_count <- apply(y, 1, function(row) {
    sum(apply(dataset[features], 1, function(d_row) all(d_row == row)))
  })

  # create the final truth table and reverse it
  y <- rev(y[multifeat_count > 0, , drop = FALSE])

  # add counts
  y$counts <- multifeat_count[multifeat_count > 0]

  return(y)
}

#' @title
#' Associate parameter names
#' @description
#' Associate parameter name to group rate logic table.
#'
#' @param this_table
#' Data.frame containing conditions for the group rate.
#' @param prefix
#' Prefix to parameter names associated to the group rate.
#' @param reference
#' Choose reference providing corresponding features.
#' Default is the group with the largest sample size.
#'
#' @return
#' input data.frame with logic table with parameter names.
#'
#' @keywords internal
associate_param_names <- function(this_table,
                                  prefix = NA,
                                  reference = NULL) {
  # column names to consider
  colmask <- setdiff(names(this_table),
                     c("counts", "param_name"))

  # extract reference index
  if (length(reference) == 0) {
    idx_refer <- which.max(this_table$counts)

  } else {
    # find row
    idx_refer <- unname(which(
      apply(this_table[colmask], 1,
            function(row) all(row == unlist(reference[colmask])))
    ))

    if (length(idx_refer) != 1) stop_h(
      "Provided reference group is not valid."
    )
  }

  # format prefix and associate parameter names
  prefix <- ifelse((!is.na(prefix)) && is.character(prefix),
                   paste0(prefix, "_"), "")

  # associate parameter names
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      paste0(prefix, paste(paste(colmask, row, sep = "_"), collapse = "_"))
    }
  )
  this_table$param_name[idx_refer] <- "reference_"

  return(this_table)
}

#' @title
#' Tidy-up Group Rate Tables
#' @description
#' Combine list with group rate definition into a single dataframe
#'
#' @param all_tables `list`\cr
#' List containing definitions for each group rate.
#'
#' @return
#' `data.frame`
#'
#' @keywords internal
tidyup_group_rate_tables <- function(all_tables) {
  # tidy version of group rate tables
  drop_counts <- function(x) x[, colnames(x) != "counts"]
  cbind_rate <- function(x, rate) cbind(group_rate = rate, x)
  df_list <- c()
  for (rate in names(all_tables)) {
    this_table <- all_tables[[rate]]
    if (length(this_table) == 0) {
      df <- list(data.frame(group_rate = rate, param_name = "reference_"))
    } else {
      df <- lapply(this_table, function(x) cbind_rate(drop_counts(x), rate))
    }
    df_list <- c(df_list, df)
  }
  
  # all unique column names reordered
  all_cols <- unique(c("group_rate", "param_name",
                       sort(unlist(lapply(df_list, names)))))
  
  # Add missing columns with NA
  combined_df <- do.call(rbind, lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    df[missing_cols] <- NA
    return(df[all_cols])
  }))
  rownames(combined_df) <- NULL
  
  return(combined_df)
}

############################################################################
############################################################################
############################# Housepi Group Rate ###########################
############################################################################
############################################################################

#' @title
#' Housepi Group Rate R6 Class
#' @description
#' Help with group rate definitions of model Rcpp object.
#'
#' @keywords internal
HousepiGroupRate <- R6::R6Class(
  classname = "HousepiGroupRate",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE ##################################
  ############################################################################

  private = list(
    # @field param_defs `list`\cr
    # List containing parameters definitions.
    .param_defs = list()
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field param_defs `list`\cr
    #' List containing parameter definitions associated to group rates.
    param_defs = function(value) {
      if (missing(value)) {
        new_defs <- list()
        for (rate_id in names(self$logic_tables)) {
          for (this_table in self$logic_tables[[rate_id]]) {
            for (par_id in setdiff(this_table$param_name, c("reference_"))) {
              new_defs[[par_id]] <- private$.param_defs[[rate_id]]
            }
          }
        }
        return(new_defs)
      }
      stop_h("'param_defs' is read only.")
    },

    #' @field param_required `character()`\cr
    #' List of parameter required by group rates.
    param_required = function(value) {
      if (missing(value)) {
        required <- c()
        for (rate_id in names(self$logic_tables)) {
          for (this_table in self$logic_tables[[rate_id]]) {
            required <- c(required, this_table$param_name)
          }
        }
        output <- setdiff(unique(required), c("reference_"))
        return(set_source_name_(output, "group rate"))
      }
      stop_h("'param_required' is read only.")
    },

    #' @field tidy_table `data.frame()`\cr
    #' Data.Frame containg the combined list with group rate and parameter definitions.
    tidy_table = function(value) {
      if (missing(value)) return(tidyup_group_rate_tables(self$logic_tables))
      stop_h("'tidy_table' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `list`\cr
    #' List containing current group rate definitions.
    defs = list(),

    #' @field logic_tables `list`\cr
    #' List containing table of parameters associated to each group rate.
    logic_tables = list(),

    #' @description
    #' New instance of class
    #' [HousepiGroupRate][housepi::HousepiGroupRate].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Find group rates required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regular_expr <- c(
        # group rate from : .group_rate("string"
        '\\.group_rate\\s*\\(\\s*["\'](.*?)["\']'
      )
      # find required group rates from regular expressions
      output <- template_required_find(refer_file, regular_expr)
      return(set_source_name_(output, "model"))
    },

    #' @description
    #' Remove a given group_rate definition.
    #'
    #' @param rate_id `character(1)`\cr
    #' The group rate definition to remove.
    #' @param show_warnings `logical(1)`\cr
    #' Wheter to show warnings.
    #' 
    #' @return
    #' `logical(1)`
    remove = function(rate_id, show_warnings = TRUE) {

      if (any(names(self$defs) == rate_id)) {
        # drop previous definition
        self$defs[rate_id] <- NULL
        private$.param_defs[rate_id] <- NULL
        self$logic_tables[rate_id] <- NULL
        return(TRUE)

      } else if (isTRUE(show_warnings)) {
        warning_h(
          sprintf("Attempt to remove unknown group_rate \"%s\".", this)
        )
      }

      return(FALSE)
    },

    #' @description
    #' Add a new group rate definition.
    #'
    #' @param name `character(1)`\cr
    #' Name of the group rate.
    #' @param by `list`\cr
    #' Feature(s) associated to group rate as dataset headers.
    #' Corresponding dataset column have to be categorical.
    #' @param prefix `character(1)`\cr
    #' Prefix to parameter names associated to the group rate.
    #' @param reference named `list`\cr
    #' Choose reference group providing corresponding features.
    #' Default is the group with the largest sample size.
    #' @param independent
    #' If `TRUE` features are considered separately.
    #' @param prior `character(1)`\cr
    #' Name of the prior distribution.
    #' @param args named `list`\cr
    #' Arguments for the prior distribution.
    #' @param rule `character(1)`\cr
    #' Name of the inference proposal rule.
    #' @param scale `double(1)`\cr
    #' Scale parameter for the proposal rule.
    #' Changed automatically during inference if `n_adpat`>0.
    #' @param min `double(1)`\cr
    #' Value of the minimum.
    #' @param max `double(1)`\cr
    #' Value of the maximum.
    #' @param show_warnings `logical(1)`\cr
    #' If `TRUE` (default) show warnings about overwriting definitions.
    #'
    #' @references
    #' Cauchemez, et al. (2009)
    #' Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
    #' New England Journal of Medicine. 361: 2619–2627. doi:10.1056/NEJMoa0905498
    add = function(name,
                   by,
                   prefix = NA,
                   reference = NULL,
                   independent = TRUE,
                   prior = "lnorm",
                   args = list(meanlog = 0, sdlog = 3),
                   min = NA,
                   max = 10,
                   rule = "lnorm",
                   scale = 1) {

      # UPDATE: controls are missing
      if (is.na(name)) stop_h(
        "Incomplete definition: missing name."
      )

      # remove previous definition if it exists
      if (self$remove(name, FALSE)) warning_h(
        sprintf("Overwriting existing \"%s\" definition.", name)
      )

      forbidden <- c("counts", "param_name", "group_rate")
      if (!all(is.na(by))) {
        if (isFALSE(is.character(by))) stop_h(
          "Argument 'by' must be of type character."
        ) 
        else if (any(by %in% forbidden)) stop_h(
          sprintf("Feature name \"%s\" is a reserved word.", by),
          "Please, rename feature and corresponding data column."
        )
      }

      # load reference
      if (!is.null(reference)) {
        if (length(by) != length(reference)) stop_h(
          "Arguments 'by' and 'reference' should have the same dimension."
        )
        else if (length(setdiff(names(reference), by)) != 0) stop_h(
          "Unrecognized names in reference list."
        )
        else reference <- setNames(as.list(reference), by)
      }

      # add new definition
      self$defs[[name]] <- list(
        by = by, prefix = prefix,
        reference = reference, independent = isTRUE(independent)
      )

      # add parameter definitions
      private$.param_defs[[name]] <- list(
        value = 1, prior = prior, args = args,
        min = min, max = max, rule = rule, scale = scale
      )

      # reorder
      y <- sort(names(self$defs))
      self$defs <- self$defs[y]
      private$.param_defs <- private$.param_defs[y]

      invisible(self)
    },

    #' @description
    #' Load group rate definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from named `list`\cr
    #' Named list of parameter definitions.
    #' @param append `logical(1)`\cr
    #' If `TRUE` appends to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_list = function(...,
                         from = NULL,
                         append = FALSE) {
      definitions <- load_definition_list(..., from = from)

      # remove previos definitions if not append
      if (isFALSE(append)) self$defs <- NULL

      for (i in seq_along(definitions)){
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])
        if (length(names(this_def)) == 0) {
          if (length(this_def) == 1) names(this_def) <- "by"
          else stop_h(
            sprintf("Wrong definition for parameter \"%s\"", def_name),
            "Provide a definition list or a scalar (if constant)."
          )
        }
        this_def$name <- def_name
        do.call(self$add, this_def)
      }

      invisible(self)
    },

    #' @description
    #' Implement group rate definition according to the dataset.
    #'
    #' @param dataset `data.frame`\cr
    #' Dataset to be considered.
    #' @param factor_info `list`\cr
    #' Levels for each factorized feature.
    implement = function(dataset, factor_info) {
      self$logic_tables <- list()

      for (rate_id in names(self$defs)) {
        this_def <- self$defs[[rate_id]] # definition of the rate
        features <- as.list(this_def[["by"]]) # data features

        new_logic_table <- list()
        if (!all(is.na(features))) {
          features <- features[!is.na(features)]

          if (isFALSE(this_def[["independent"]])) features <- list(features)

          # default prefix
          prefix <- as.character(this_def[["prefix"]])
          if (is.na(prefix)) {
            prefix <- substring(format_string(tolower(rate_id)), 1, 3)
          }

          for (this_feat in features) {
            # 1. create logic table with group rates
            this_table <- group_rate_logic_table(dataset, this_feat)
            # associate names for parameters
            this_table <- associate_param_names(
              this_table, prefix, this_def$reference
            )
            new_logic_table[[length(new_logic_table) + 1]] <- this_table

            # 2. create relative cpp script chunk
            # factorize according to external factor information
            for (col in names(factor_info)) {
              lv <- factor_info[[col]]$levels
              if (is.null(lv)) stop_h(
                "Wrong group rate definition: \"", col, "\" has no levels."
              )

              if (col %in% colnames(this_table)) {
                this_table[[col]] <- factor(this_table[[col]], 
                                            levels = lv, 
                                            labels = seq_along(lv))
              }
            }
            # concatenate new lines to the script
            # mask reference case
            mask <- this_table$param_name == "reference_"
          }
        }

        # store the new rate logic table
        self$logic_tables[[rate_id]] <- new_logic_table
      }

      invisible(self)
    },

    #' @description
    #' Assert that group rate definitions match the requirements.
    #'
    #' @param required `character()`\cr
    #' Vector containing required group rates.
    assert = function(required) {
      # load definitions
      defined <- names(self$defs)

      # UPDATE
      # point out group rates that are not used
      notused <- template_not_used_find(required, defined)
      # drop unused group rates
      if (length(notused) > 0) {
        for (this in notused) self$remove(this)
        stop_h(
          sprintf("Dropped unused group rate%s %s",
                  ifelse(length(notused) > 1, "s:", ""),
                  character_to_string_list(notused)),
          find_closest_match(notused, required, min_dist = 3)
        )
      }

      # find missing group rates and set as default
      for (g in template_missing_find(required, defined)) {
        self$add(g, by = NA)
      }

      invisible(self)
    }
  )
)
