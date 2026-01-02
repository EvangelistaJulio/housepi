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
# Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# parameter.R
#

# Helper to associate prior definition.
parameter_helper_prior <- function(opts, prior, args) {
  # prior to tell if constant
  if (is.na(prior)) stop_h(
    sprintf("No prior provided for parameter '%s'.", opts$name)
  )
  opts$const <- FALSE
  # associate prior arguments, if any
  if (is.null(args)) opts$prior <- prior
  else opts$prior <- housepi_combine_to_call(prior, args)
  opts
}

# Helper with Monte Carlo proposal.
parameter_helper_proposal <- function(opts, discr, proposal, proposal_scale) {
  opts$discr <- isTRUE(discr)

  # default proposal
  opts$proposal <- "sample"
  if (!is.na(proposal)) {
    opts$proposal <- proposal
    opts$proposal_scale <- ifelse(is.na(proposal_scale), 1, proposal_scale)
  }
  opts
}

############################################################################
############################################################################
############################# Housepi Parameter ############################
############################################################################
############################################################################

#' @title Housepi Parameter
#' @description
#' Help with parameter definitions of model Rcpp object.
#'
#' @keywords internal
#' @noRd
HousepiParameter <- R6::R6Class(
  classname = "HousepiParameter",
  cloneable = TRUE,

  private = list(
    deep_clone = function(name, value) helper_clone(value)
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)` with the module family.
    family = "parameter",
    #' @field store_warn_over `character()` with cached overwritten names.
    store_warn_over = character(),
    #' @field defs `data.frame` with formatted parameter definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiParameter][housepi::HousepiParameter].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description
    #' Find parameters required by the model.
    #'
    #' @param refer_file
    #' Path/to/.refer source file.
    require = function(refer_file) {
      # vector of regular expression to use for assertion
      regexp <- c(
        # from : par("string" or get_param('string'
        'par\\s*\\(\\s*["\'](.*?)["\']',
        # from : param_pointer("string" or param_pointer('string'
        'param_pointer\\s*\\(\\s*["\'](.*?)["\']',
        # from : emplace_contact_rate("string"
        'emplace_contact_rate\\s*\\(\\s*["\'](.*?)["\']',
        # from : emplace_group_rate(..., "string"
        'emplace_group_rate\\s*\\(\\s*[^,]*\\s*,\\s*["\'](.*?)["\']'
      )
      # retrieve required parameters from regular expressions
      set_source_name_(template_required_find(refer_file, regexp),
                       "model")
    },

    #' @description
    #' Remove a given parameter.
    #' @param name
    #' A `character(1)` with the name of parameter definitions to remove.
    #' @param warn_removal
    #' `logical(1)`to show warnings for failed removal.
    remove = function(name, warn_removal = FALSE) {
      mask <- self$defs[["name"]] == name
      if (any(mask)) {
        self$defs <- self$defs[!mask, ]
        return(name)
      } else if (isTRUE(warn_removal)) {
        warning_h(sprintf("Attempt to remove unknown parameter '%s'.", name))
      }
      c()
    },

    #' @inherit add_parameter
    add_quo = function(name, value = NULL,
                       const = NA, discr = FALSE,
                       prior = NA_character_, ...,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       proposal_scale = NA_real_,
                       .warn_over = TRUE) {
      # load non-quosures
      name <- housepi_smart_unwrap(name, "name")
      const <- housepi_smart_unwrap(const, "const")
      discr <- housepi_smart_unwrap(discr, "discr")
      prior <- housepi_smart_unwrap(prior, "prior")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      proposal_scale <- housepi_smart_unwrap(proposal_scale, "proposal_scale")
      # load quosures
      value <- housepi_extractor(!!value)
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)
      # load dots # UPDATE: not super elegant...
      dots <- list(...)
      # flatten if single quosures list (from enquos)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      args <- housepi_extractor(!!!dots)

      # control valid name
      assert_valid_definition_name(name, self$family, !is.na(value))
      if (all(is.null(value), is.na(prior))) stop_h(
        sprintf("Incomplete definition for parameter '%s'.", name),
        "Provide at least one between 'value' and 'prior'."
      )

      opts <- list(name = as.character(name), value = NA_real_,
                   const = TRUE, discr = FALSE,
                   prior = NA_character_,
                   min = NA_real_, max = NA_real_,
                   proposal = NA_character_,
                   proposal_scale = NA_real_)

      # drop previous definition eventually
      rmv_name <- self$remove(opts$name)
      if (.warn_over) template_overwrite_warn(rmv_name, self$family)
      else self$store_warn_over <- c(self$store_warn_over, rmv_name)

      # assume it's a constant if a prior is not provided
      if (is.na(const)) const <- is.na(prior)

      # parse initial value
      if (is.null(value)) opts$value <- NA_real_
      else opts$value <- value

      if (isTRUE(const)) {
        # check that a value is provided
        if (is.na(opts$value)) stop_h(
          sprintf("Missing value for constant parameter '%s'.", opts$name)
        )
      } else {
        opts <- parameter_helper_prior(opts, prior, args)
        if (!is.null(min)) opts$min <- min
        if (!is.null(max)) opts$max <- max
        opts <- parameter_helper_proposal(opts, discr, proposal, proposal_scale)
      }

      # append as a row
      self$defs <- rbind(self$defs, data.frame(opts))

      invisible(self)
    },

    #' @inherit load_parameter_list
    load_quosures = function(xquos, append = FALSE) {
      if (isFALSE(append)) self$initialize()
      housepi_load_defin_quos(self, "value", xquos)
      invisible(self)
    },

    #' @description
    #' Assert parameter definitions match the requirement list.
    #'
    #' @param required `character()` containing required parameters.
    #' @param capital `logical(1)` to raise error.
    assert = function(required, capital) {
      defined <- self$defs[["name"]]
      # point out parameters that are not used
      template_not_used_warn(required, defined, self$family)
      # stop if parameters are missing
      sms <- template_missing_stop(required, defined, self$family)
      if (isTRUE(capital) && length(sms) > 0) stop_h(sms)
      else sms
    }
  )
)
