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
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# augmented.R
#

# Parse augmented implicit function.
helper_augmented_function <- function(opts, func, func_args) {
  if (is.numeric(showcaser_extract_numeric(func))) {
    # implicit numeric definition
    opts$value <- func
    if (is.na(opts$const)) {
      opts$const <- all(is.na(c(opts$model, opts$proposal)))
    }
  } else {
    # distinguish purely implicit function to plain name
    opts$func <- ifelse(
      is.null(func_args), func,
      housepi_combine_to_call(func, func_args)
    )
    # note: a priori a function is not constant
    if (is.na(opts$const)) opts$const <- FALSE
  }

  opts
}

# process model arguments
helper_augmented_model <- function(opts, args) {
  if (is.na(opts$model)) stop_h(
    sprintf("Spare arguments for augmented '%s'.%s", opts$name,
            ifelse(is.na(opts$func), "", "Did you mean 'func_args'?"))
  )
  if (!is.null(args)) {
    opts$model <- housepi_combine_to_call(opts$model, args)
  }
  opts
}

# Interpret user choice for Monte Carlo proposal.
augmented_helper_proposal <- function(opts, proposal_scale) {
  # sampling proposal but no constant parameter provided
  if (is.na(opts$const) && !(is.na(opts$proposal) && is.na(opts$model))) {
    opts$const <- FALSE
  }

  if (isFALSE(opts$const)) {
    if (is.na(opts$proposal)) {
      if (is.na(opts$model)) opts$const <- TRUE
      else opts$proposal <- "sample"
    } else if (opts$proposal == "sample") {
      if (is.na(opts$model)) stop_h(
        sprintf("No model to sample from for augmented '%s'.", opts$name)
      )
    } else {
      if (is.na(proposal_scale)) opts$proposal_scale <- 1
      else if (is.numeric(proposal_scale) &&
                 proposal_scale > 0) opts$proposal_scale <- proposal_scale
      else stop_h("Proposal scale must be a numeric greater than 0.")
    }
  } else {
    opts$proposal <- NA_character_
  }
  opts
}

# Process constant or random.
augmented_helper_const_vs_rand <- function(opts, value) {
  if (!is.null(value)) opts$value <- value

  if (isFALSE(opts$const) && all(is.na(c(opts$value, opts$model)))) stop_h(
    sprintf("Impossible to randomly initialize undefined augmented '%s'.",
            opts$name),
    "Provide one between value and model."
  )
  else if (is.na(opts$const)) opts$const <- TRUE

  opts
}

#' @title Housepi Augmented
#' @description
#' Help with augmented variable definitions of model Rcpp object.
#'
#' @keywords internal
#' @noRd
HousepiAugmented <- R6::R6Class(
  classname = "HousepiAugmented",
  cloneable = TRUE,

  ############################################################################
  ################################# PRIVATE ##################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value)
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field augm_required `character()`\cr
    #' Find augmented required by augmented variable definitions.
    augm_required = function(value) {
      if (!missing(value)) stop_h("'augm_required' is read only.")
      augmented_find_required_def(self$defs, "A")
    },

    #' @field data_required `character()`\cr
    #' Find data required by augmented variable definitions.
    data_required = function(value) {
      if (!missing(value)) stop_h("'data_required' is read only.")
      augmented_find_required_def(self$defs, "D")
    },

    #' @field group_rates_required `character()`\cr
    #' Find group rates required by augmented variable definitions.
    group_rates_required = function(value) {
      if (!missing(value)) stop_h("'group_rates_required' is read only.")
      augmented_find_required_def(self$defs, "G")
    },

    #' @field param_required `character()`\cr
    #' Find parameters required by augmented variable definitions.
    param_required = function(value) {
      if (!missing(value)) stop_h("'param_required' is read only.")
      augmented_find_required_def(self$defs, "P")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field defs `data.frame`\cr
    #' Augmented variable formatted definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiAugmented][housepi::HousepiAugmented].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Find augmented variables required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regexp <- c(
        # from : .aug("string"
        '.aug\\s*\\(\\s*["\'](.*?)["\']',
        # from : get_augmented("string"
        '.get_augm_variable\\s*\\(\\s*["\'](.*?)["\']',
        # from : .set_augm_value("string"
        '.set_augm_value\\s*\\(\\s*["\'](.*?)["\']'
      )
      # find required augmented variables from regular expressions
      output <- template_required_find(refer_file, regexp)

      c(set_source_name_(output, "model"), self$augm_required)
    },

    #' @description
    #' Remove a given augmented variable.
    #'
    #' @param name_lst
    #' A `character()` with the name of augmented definitions to remove.
    #' @param status
    #' An `integer()` with the associated status of the definition.
    #' @param warn_removal
    #' `logical(1)`to show warnings for failed removal.
    #' @param warn_overwrite
    #' `logical(1)`to show warnings in case of overwriting.
    remove = function(name_lst, status_lst,
                      warn_removal = FALSE,
                      warn_overwrite = FALSE) {
      stopifnot(length(name_lst) == length(status_lst))
      stopifnot(is.list(name_lst))
      stopifnot(is.list(status_lst))

      all_rmv <- c()
      # loop over names
      for (i in seq_along(name_lst)) {
        nm <- name_lst[[1]]
        st <- status_lst[[1]]
        mask <- self$defs[["name"]] == nm
        if (any(mask)) {
          # check if same name and status
          mask <- mask & sapply(self$defs[["status"]], function(x) {
            identical(x, as.integer(st))
          })

          if (any(mask)) {
            self$defs <- self$defs[!mask, ]
            all_rmv <- c(all_rmv, nm)
          } else if (isTRUE(warn_removal)) {
            warning_h(sprintf(
              "Attempt to remove unknown augmented '%s' for status %d.", nm, st
            ))
          }
        } else if (isTRUE(warn_removal)) {
          warning_h(sprintf("Attempt to remove unknown augmented '%s'.", nm))
        }
      }

      # warning about overwriting.
      if (warn_overwrite && (length(all_rmv) > 0)) warning_h(sprintf(
        "Overwriting augmented definition%s%s",
        ifelse(length(all_rmv) > 1, "s for:\n", " for "),
        character_to_string_list(all_rmv, NA)
      ))

      invisible(self)
    },

    #' @inherit add_augmented
    add_quo = function(name = NA_character_, value = NULL,
                       status = NA_integer_, const = NA, discr = FALSE,
                       func = NULL, func_args = NULL,
                       model = NA_character_, ...,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       proposal_scale = NA_real_) {

      # load non-quosures
      name <- housepi_smart_unwrap(name, "name")
      status <- housepi_smart_unwrap(status, "status")
      const <- housepi_smart_unwrap(const, "const")
      discr <- housepi_smart_unwrap(discr, "discr")
      model <- housepi_smart_unwrap(model, "model")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      proposal_scale <- housepi_smart_unwrap(proposal_scale, "proposal_scale")
      # load quosures
      value <- housepi_extractor(!!value)
      func <- housepi_extractor(!!func)
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)
      # load dots # UPDATE: not super elegant...
      dots <- list(...)
      # flatten dots if single quosures list (from enquos)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      args <- housepi_extractor(!!!dots)

      # flatten function arguments if single quosures list (from enquos)
      # NOTE: lazy fix of rlang warning...
      if (rlang::is_quosures(func_args)) {
        func_args <- housepi_extractor(!!!func_args)
      } else {
        func_args <- housepi_extractor(!!func_args)
      }

      # control valid name
      assert_valid_definition_name(name, "augmented", !is.na(func))
      # control valid status
      assert_valid_status(status)
      # control meaningful definition
      if (all(is.null(c(value, func)), is.na(model))) stop_h(
        sprintf("Incomplete definition for augmented '%s'.", name),
        "Provide at least one between 'value', 'model' and 'func'."
      )

      # auxiliary typed options for Rcpp build method
      opts <- list(name = as.character(name),
                   value = NA_real_,
                   status = status,
                   const = as.logical(const),
                   discr = isTRUE(discr),
                   func = NA_character_,
                   fcumul = FALSE, # UPDATE: can I make this an option?
                   model = as.character(model),
                   min = NA_real_,
                   max = NA_real_,
                   proposal = as.character(proposal),
                   proposal_scale = NA_real_)

      # drop previous definition eventually
      self$remove(list(opts$name), list(opts$status),
                  warn_removal = FALSE, warn_overwrite = TRUE)

      # model arguments
      if (!is.null(args)) {
        opts <- helper_augmented_model(opts, args)
      }
      # augmented function
      if (!is.null(func)) {
        opts <- helper_augmented_function(opts, func, func_args)
      }
      # augmented variable
      if (is.na(opts$func)) {
        opts <- augmented_helper_proposal(opts, proposal_scale)
        opts <- augmented_helper_const_vs_rand(opts, value)
      }

      # domain
      if (!is.null(min)) opts$min <- min
      if (!is.null(max)) opts$max <- max

      # append definition(s)
      self$defs <- rbind(self$defs, data.frame(opts))

      invisible(self)
    },

    #' @inherit load_augmented_list
    load_quosures = function(xquos, append = FALSE) {
      # UPDATE: append requires checking existing definitions.
      # I can store old ones, compare them and then drop at the end
      # remove previous definitions if append is false
      if (isFALSE(append)) self$initialize()
      self <- housepi_load_defin_quos(self, "func", xquos)

      invisible(self)
    },

    #' @description
    #' Assert augmented definitions from the model implementation script.
    #'
    #' @param required `character()` containing required augmented
    #' @param status_set `integer()` with all known statuses in data.
    #' @param capital `logical(1)` to raise error.
    assert = function(required, status_set, capital) {
      status_set <- sort(status_set)
      defined <- self$defs[["name"]]
      # point out augmented variables that are not used
      template_not_used_warn(required, defined, "augmented")

      # stop if augmented variables are missing
      sms <- template_missing_stop(required, defined, "augmented")
      if (isTRUE(capital) && length(sms) > 0) stop_h(sms)
      else sms
    }
  )
)
