## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @title Internal R6 holder for housepi augmented definitions
#' @description
#' Helper class managing augmented definitions of the model Rcpp object.
#' @keywords internal
#' @noRd
HousepiAugmented <- R6::R6Class(
  classname = "HousepiAugmented",
  cloneable = TRUE,

  ############################################################################
  ################################# PRIVATE ##################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Definitions required by augmented definitions
    # according to the "type" of the definition.
    find_required_def_ = function(type) {
      find_required_definitions(
        self$defs, type, c("value", "func", "model", "min", "max"), self$family
      )
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field augm_required `character()`.
    #' Find augmented required by augmented variable definitions.
    augm_required = function(value) {
      if (!missing(value)) stop_read_only("augm_required")
      private$find_required_def_("A")
    },

    #' @field data_required `character()`.
    #' Find data required by augmented variable definitions.
    data_required = function(value) {
      if (!missing(value)) stop_read_only("data_required")
      private$find_required_def_("D")
    },

    #' @field group_rates_required `character()`.
    #' Find group rates required by augmented variable definitions.
    group_rates_required = function(value) {
      if (!missing(value)) stop_read_only("group_rates_required")
      private$find_required_def_("G")
    },

    #' @field param_required `character()`.
    #' Find parameters required by augmented variable definitions.
    param_required = function(value) {
      if (!missing(value)) stop_read_only("param_required")
      private$find_required_def_("P")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Module family.
    family = "augmented",
    #' @field store_warn_over `character()`. Cached overwritten names.
    store_warn_over = character(),
    #' @field defs `data.frame`. Formatted augmented definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiAugmented][housepi::HousepiAugmented].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description Find augmented variables required by model.
    #' @param refer_file Path to reference file.
    require = function(refer_file) {
      output <- template_extract_keys(refer_file, "AugmKey")
      c(set_source_name_(output, "model"), self$augm_required)
    },

    #' @description Remove a given augmented definition.
    #' @param name `character(1)`. Augmented definitions to remove.
    #' @param status `integer(1)`. Associated status of definition.
    #' @param warn_removal `logical(1)`. Show warnings for failed removal.
    remove = function(name, status, warn_removal = FALSE) {
      mask <- self$defs[["name"]] == name
      if (any(mask)) {
        # check if same name and status
        mask <- mask & sapply(self$defs[["status"]], function(x) {
          identical(x, as.integer(status))
        })
        if (any(mask)) {
          self$defs <- self$defs[!mask, ]
          return(name)
        } else if (isTRUE(warn_removal)) {
          warn_h(sprintf(
            "Attempt to remove unknown augmented '%s' for status %d.",
            name, status
          ))
        }
      } else if (isTRUE(warn_removal)) {
        warn_h(sprintf("Attempt to remove unknown augmented '%s'.", name))
      }

      c()
    },

    #' @inherit add_augmented
    add_quo = function(name = NA_character_, value = NULL,
                       status = NA_integer_,
                       const = NA, discr = FALSE,
                       model = NA_character_, ...,
                       func = NULL, func_args = NULL,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       proposal_scale = NA_real_,
                       .warn_over = TRUE,
                       .assert = FALSE) {
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

      # flatten function arguments if quosures
      if (rlang::is_quosure(func_args)) {
        func_args <- housepi_extractor(!!func_args)
      } else {
        func_args <- housepi_extractor(!!!func_args)
      }

      # control valid name
      assert_valid_definition_name(name, self$family, !is.na(value))
      # control valid status
      assert_valid_status(status)
      # control meaningful definition
      if (all(is.null(c(value, func)), is.na(model))) stop_h(
        sprintf("Incomplete definition for augmented '%s'.", name),
        "Provide at least one between 'value', 'model' and 'func'."
      )
      if (!is.null(value) && !is.null(func)) stop_h(
        sprintf("Wrong definition for augmented '%s'.", name),
        "Provide only one between 'value' and 'func'."
      )

      # auxiliary typed options for Rcpp build method
      opts <- list(name = as.character(name),
                   value = NA_real_,
                   status = as.integer(status),
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
      rmv_name <- self$remove(opts$name, opts$status)
      if (.warn_over) template_overwrite_warn(rmv_name, self$family)
      else self$store_warn_over <- c(self$store_warn_over, rmv_name)

      opts <- helper_variable_core(opts, self$family, value, func, func_args)
      opts <- helper_variable_model(opts, self$family, args)
      opts <- helper_variable_move(opts, self$family, proposal_scale)

      # domain
      if (!is.null(min)) opts$min <- min
      if (!is.null(max)) opts$max <- max

      # append definition(s)
      self$defs <- rbind(self$defs, data.frame(opts))
      # assert last definition is correct
      if (isTRUE(.assert)) housepi_assert_args(self, last = TRUE)

      invisible(self)
    },

    #' @inherit add_augmented_list
    add_quo_list = function(xquos, append = FALSE) {
      if (isFALSE(append)) self$initialize()
      housepi_load_defin_quos(self, "value", xquos)
      housepi_assert_args(self)
      invisible(self)
    },

    #' @description
    #' Assert augmented definitions from the model implementation script.
    #' @param required `character()`. Required definitions.
    #' @param status_set `integer()`. All known statuses in data.
    #' @param capital `logical(1)`. Whether to raise error.
    #' @param warn_unused `logical(1)`. Warn if unused are found.
    assert = function(required, status_set, capital, warn_unused) {
      status_set <- sort(status_set)
      defined <- self$defs[["name"]]
      # point out augmented variables that are not used
      if (isTRUE(warn_unused)) {
        template_not_used_warn(required, defined, self$family)
      }
      # stop if augmented variables are missing
      sms <- template_missing_stop(required, defined, self$family)
      if (isTRUE(capital) && length(sms) > 0) stop_h(sms)
      else sms
    }
  )
)
