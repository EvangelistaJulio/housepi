#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
# See the COPYRIGHT file for details.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the COPYING file for details.
#
# group_rates.R
#

#' @title
#' Assert data features
#' @description
#' Assert that features are present in data.
#'
#' @param features
#' Features as list of data column names.
#' @param headers
#' All column names of the data.
#'
#' @keywords internal
#' @noRd
assert_data_features <- function(features, headers) {
  # assert that features are present
  missing <- c()
  for (this in features) {
    if (!(this %in% headers)) missing <- c(missing, this)
  }
  if (length(missing) > 0) stop_h(
    sprintf(
      "Unknown feature%s %s",
      ifelse(length(missing) > 1, "s:", ""),
      character_to_string_list(missing)
    ),
    "Please check your data headers."
  )
}

#' @title
#' Group rate logic table
#' @description
#' Return logic table of group rates
#'
#' @param data_wide
#' Data from which to extract group rates
#' @param by
#' Features as list of data column names.
#'
#' @return
#' data.frame with logic table
#'
#' @keywords internal
#' @noRd
group_rate_logic_table <- function(data_wide, by) {
  # sort alphabetically the features
  features <- sort(as.character(by))

  # convert features
  assert_data_features(features, colnames(data_wide))

  # extract unique features and sort by list length
  uni_features <- lapply(features, function(f) sort(unique(data_wide[[f]])))
  sorted_indices <- order(sapply(uni_features, length), decreasing = TRUE)
  features <- features[sorted_indices]
  uni_features <- uni_features[sorted_indices]
  names(uni_features) <- features

  # create a dataframe of all possible multi-features
  y <- expand.grid(uni_features)

  # check that the multi-feature exists
  multifeat_count <- apply(y, 1, function(row) {
    sum(apply(data_wide[features], 1, function(d_row) all(d_row == row)))
  })

  # create the final truth table and reverse it
  y <- y[multifeat_count > 0, , drop = FALSE]

  # add counts
  y$n <- multifeat_count[multifeat_count > 0]

  y
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
#' @noRd
associate_param_names <- function(this_table, prefix, reference) {
  # column names to consider
  colmask <- setdiff(names(this_table), c("n", "param_name"))

  # extract reference index
  if (length(reference) == 0) {
    idx_refer <- which.max(this_table$n)

  } else {
    # find row
    idx_refer <- unname(which(
      apply(this_table[colmask], 1,
            function(row) all(row == unlist(reference[colmask])))
    ))

    if (length(idx_refer) != 1) stop_h("Provided reference group is not valid.")
  }

  # format prefix and associate parameter names
  prefix <- ifelse(
    (prefix == "") || (!is.character(prefix)),
    "", paste0(prefix, "_")
  )

  # associate parameter names
  this_table$param_name <- apply(
    this_table[colmask], 1, function(row) {
      paste0(prefix, paste(paste(colmask, row, sep = "_"), collapse = "_"))
    }
  )
  this_table$param_name[idx_refer] <- ".reference"

  this_table
}

#' @title
#' Tidy-up Group Rate Tables
#' @description
#' Combine list with group rate definition into a single dataframe
#'
#' @param all_tables
#' A `list` containing definitions for each group rate.
#'
#' @return
#' `data.frame`
#'
#' @keywords internal
#' @noRd
tidyup_group_rate_tables <- function(all_tables) {
  # tidy version of group rate tables
  cbind_rate <- function(x, rate) cbind(group_rate = rate, x)
  df_list <- c()
  for (rate in names(all_tables)) {
    this_table <- all_tables[[rate]]
    if (length(this_table) == 0) {
      df <- list(data.frame(group_rate = rate, param_name = ".reference"))
    } else {
      df <- lapply(this_table, function(x) cbind_rate(x, rate))
    }
    df_list <- c(df_list, df)
  }

  # all unique column names reordered
  all_cols <- unique(c("group_rate", "param_name", "n",
                       sort(unlist(lapply(df_list, names)))))

  # Add missing columns with NA
  combined_df <- do.call(rbind, lapply(df_list, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    df[missing_cols] <- NA_character_
    df[all_cols]
  }))
  rownames(combined_df) <- NULL

  combined_df
}

############################################################################
############################################################################
############################# Housepi Group Rate ###########################
############################################################################
############################################################################

#' @title
#' Housepi Group Rate
#' @description
#' Help with group rate definitions of model Rcpp object.
#'
#' @keywords internal
#' @noRd
HousepiGroupRate <- R6::R6Class(
  classname = "HousepiGroupRate",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE ##################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # @field parameters `list`\cr
    # List containing parameters definitions.
    .parameters = list(),

    # @title
    # Default MCMC Setup
    # @description
    # Return list with default MCMC setup if user choice is missing.
    #
    # @return `list`
    #
    # @references
    # Cauchemez, et al. (2009)
    # Household Transmission of 2009 Pandemic Influenza A (H1N1) Virus in the United States.
    # New England Journal of Medicine. 361: 2619-2627. doi:10.1056/NEJMoa0905498
    #
    # @keywords internal
    # @noRd
    .get_default_mcmc_setup = function(prior, ...,
                                       min, max,
                                       proposal, proposal_scale) {
      args <- unlist(list(...))
      # check wheter default configuration
      is_default <- sapply(c(prior, args, min, max, proposal, proposal_scale),
                           function(x) is.na(x) || is.null(x))
      # if all defaults, provide default choice
      if (all(is_default)) return(
        list(value = 1, prior = "lognormal",
             meanlog = 0, sdlog = 3,
             min = NA_real_, max = 10,
             proposal = "lognormal", proposal_scale = 1)
      )
      # else go with user choice
      list(value = 1, prior = prior, ...,
           min = min, max = max,
           proposal = proposal, proposal_scale = proposal_scale)
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field parameters `list`\cr
    #' List containing parameter definitions associated to group rates.
    parameters = function(value) {
      if (!missing(value)) stop_h("'parameters' is read only.")
      new_params <- list()
      for (rate_id in names(self$logic_tables)) {
        for (this_table in self$logic_tables[[rate_id]]) {
          for (par_id in setdiff(this_table$param_name, c(".reference"))) {
            new_params[[par_id]] <- private$.parameters[[rate_id]]
          }
        }
      }
      new_params
    },

    #' @field param_required `character()`\cr
    #' List of parameter required by group rates.
    param_required = function(value) {
      if (!missing(value)) stop_h("'param_required' is read only.")
      required <- c()
      for (rate_id in names(self$logic_tables)) {
        for (this_table in self$logic_tables[[rate_id]]) {
          required <- c(required, this_table$param_name)
        }
      }
      output <- setdiff(unique(required), c(".reference"))
      set_source_name_(output, "group rate")
    },

    #' @field tidy_table `data.frame()`\cr
    #' Data.Frame containg the combined list
    #' with group rate and parameter definitions.
    tidy_table = function(value) {
      if (!missing(value)) stop_h("'tidy_table' is read only.")
      tidyup_group_rate_tables(self$logic_tables)
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
    #' New object of class
    #' [HousepiGroupRate][housepi::HousepiGroupRate].
    initialize = function() {
      self$defs <- list()
      self$logic_tables <- list()
      .parameters <- list()
      invisible(self)
    },

    #' @description
    #' Find group rates required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regexp <- c(
        # group rate from : .group_rate("string"
        '\\.group_rate\\s*\\(\\s*["\'](.*?)["\']'
      )
      # find required group rates from regular expressions
      output <- template_required_find(refer_file, regexp)
      set_source_name_(output, "model")
    },

    #' @description
    #' Remove a given group_rate definition.
    #'
    #' @param rate_id
    #' A `character(1)` with the group rate definition to remove.
    #' @param show_warnings
    #' `logical(1)` to show warnings.
    #' Default is `TRUE`.
    #'
    #' @return
    #' `logical(1)`
    remove = function(rate_id, show_warnings = FALSE) {
      if (any(names(self$defs) == rate_id)) {
        # check if old definition was just default
        is_old_default <- all(is.na(self$defs[rate_id]$by))

        # drop previous definition
        self$defs[rate_id] <- NULL
        self$logic_tables[rate_id] <- NULL
        private$.parameters[rate_id] <- NULL
        return(!is_old_default)

      } else if (isTRUE(show_warnings)) {
        warning_h(sprintf("Attempt to remove unknown group_rate '%s'.", this))
      }

      FALSE
    },

    #' @inherit add_group_rate
    add = function(name, by, prefix = NA_character_,
                   reference = NULL, independent = TRUE,
                   prior = NA_character_, ...,
                   min = NA_real_, max = NA_real_,
                   proposal = NA_character_,
                   proposal_scale = NA_real_) {
      # control valid name
      assert_valid_definition_name(name, "group rate")

      # remove previous definition if it exists
      if (self$remove(name)) warning_h(
        sprintf("Overwriting group rate definition for '%s'.", name)
      )

      forbidden <- c("n", "param_name", "group_rate")
      if (!all(is.na(by))) {
        if (isFALSE(is.character(by))) stop_h(
          "Argument 'by' must be of type character."
        )
        else if (any(by %in% forbidden)) stop_h(
          sprintf("Feature name '%s' is a reserved word.", by),
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
        reference <- as.list(reference)
        reference <- stats::setNames(reference, by)
      }

      # add new definition
      self$defs[[name]] <- list(
        by = by, prefix = prefix,
        reference = reference, independent = isTRUE(independent)
      )

      # add parameter definitions
      private$.parameters[[name]] <- private$.get_default_mcmc_setup(
        prior = prior, ..., min = min, max = max,
        proposal = proposal, proposal_scale = proposal_scale
      )

      # reorder
      y <- sort(names(self$defs))
      self$defs <- self$defs[y]
      private$.parameters <- private$.parameters[y]

      invisible(self)
    },

    #' @description
    #' Load group rate definitions either as direct arguments
    #' or from predefined list.
    #'
    #' @param ...
    #' Named arguments with respective list of definition details.
    #' @param from
    #' A named `list` with parameter definitions.
    #' @param append
    #' `logical(1)` to append to current definitions.
    #' If `FALSE` drop all pre-existing definitions (default).
    load_list = function(..., from = NULL, append = FALSE) {
      definitions <- helper_load_definition_list(..., from = from)

      # remove previous definitions if not append
      if (isFALSE(append)) self$initialize()
      all_removed <- c()

      for (i in seq_along(definitions)){
        this_def <- as.list(definitions[[i]])
        def_name <- names(definitions[i])
        if (length(names(this_def)) == 0) {
          if (length(this_def) == 1) names(this_def) <- "by"
          else stop_h(
            sprintf("Wrong definition for parameter '%s'", def_name),
            "Provide a definition list or a scalar (if constant)."
          )
        }
        this_def$name <- def_name
        if (self$remove(def_name)) all_removed <- c(all_removed, def_name)
        do.call(self$add, this_def)
      }

      if (length(all_removed) > 0) warning_h(sprintf(
        "Overwriting group rate definition%s%s",
        ifelse(length(all_removed) > 1, "s for:\n", " for "),
        character_to_string_list(all_removed, NA)
      ))

      invisible(self)
    },

    #' @description
    #' Implement group rate definition according to data.
    #'
    #' @param data_wide
    #' A `data.frame` with the dataset to be considered.
    #' @param factor_info
    #' A `list` with the levels for each factorized feature.
    implement = function(data_wide, factor_info) {
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
            prefix <- substring(tolower(rate_id), 1, 3)
          }

          for (this_feat in features) {
            # 1. create logic table with group rates
            this_table <- group_rate_logic_table(data_wide, this_feat)
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
                "Wrong group rate definition: '", col, "' has no levels."
              )

              if (col %in% colnames(this_table)) {
                this_table[[col]] <- factor(this_table[[col]],
                                            levels = lv,
                                            labels = seq_along(lv))
              }
            }
            # concatenate new lines to the script
            # mask reference case
            mask <- this_table$param_name == ".reference"
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
    #' @param required `character()` containing the required group rates.
    #' @param capital `logical(1)` to raise error.
    assert = function(required, capital) {
      # load definitions
      defined <- names(self$defs)

      # UPDATE
      # point out group rates that are not used
      notused <- template_not_used_find(required, defined)
      # drop unused group rates
      if (length(notused) > 0) {
        sms <- c(
          sprintf("Dropped unused group rate%s %s",
                  ifelse(length(notused) > 1, "s:", ""),
                  character_to_string_list(notused)),
          find_closest_match(notused, required, min_dist = 3)
        )
        if (isTRUE(capital)) {
          for (which in notused) self$remove(which)
          stop_h(sms)
        }
        warning_h(sms)
      }

      # find missing group rates and set as default
      if (isTRUE(capital)) for (g in template_missing_find(required, defined)) {
        self$add(g, by = NA_character_)
      }
      template_missing_stop(required, defined, "group rates")
    },

    #' @description
    #' Tidy table of group rates definition to be passed to Cpp.
    #'
    #' @param info
    #' A `list` of `list`s with info about class and levels of data columns.
    #'
    expose_to_cpp = function(info) {
      defs <- self$tidy_table
      defs$n <- NULL
      for (k in intersect(colnames(defs), names(info))) {
        defs[[k]] <- as.integer(factor(defs[[k]], levels = info[[k]]$levels))
      }
      defs
    }
  )
)
