## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @title Internal R6 holder for housepi parameter definitions
#' @description
#' Helper class managing parameter definitions of the model object.
#' @keywords internal
#' @noRd
HousepiParameter <- R6::R6Class(
  classname = "HousepiParameter",
  cloneable = TRUE,

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    # Definitions required by parameter definitions
    # according to the "type" of the definition.
    find_required_def_ = function(type) {
      find_required_definitions(
        self$defs, type, c("value", "func", "min", "max"), self$family
      )
    }
  ),

  active = list(
    #' @field param_required `character()`.
    #' Find parameters required by parameter function definitions.
    param_required = function(value) {
      if (!missing(value)) stop_read_only("param_required")
      private$find_required_def_("P")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Holder family.
    family = "parameter",

    #' @field store_warn_over `character()`. Cached overwritten names.
    store_warn_over = character(),

    #' @field defs `data.frame`. Formatted parameter definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiParameter][housepi::HousepiParameter].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description Find parameters required by model.
    #' @param refer_file Path to reference file.
    require = function(refer_file) {
      set_source_name_(template_extract_keys(refer_file, "ParamKey"),
                       "model")
    },

    #' @description Remove a given parameter.
    #' @param name `character(1)`. Parameter definition to remove.
    #' @param warn_removal `logical(1)`. Show warnings for failed removal.
    remove = function(name, warn_removal = FALSE) {
      mask <- self$defs[["name"]] == name
      if (any(mask)) {
        self$defs <- self$defs[!mask, ]
        return(name)
      } else if (isTRUE(warn_removal)) {
        warn_h(sprintf("Attempt to remove unknown parameter '%s'.", name))
      }
      c()
    },

    #' @inherit add_parameter
    add_quo = function(name, value = NULL,
                       const = NA, discr = FALSE,
                       prior = NA_character_, ...,
                       func = NULL, func_args = NULL,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       proposal_scale = NA_real_,
                       .warn_over = TRUE,
                       .assert = FALSE) {
      # load non-quosures
      name <- housepi_smart_unwrap(name, "name")
      const <- housepi_smart_unwrap(const, "const")
      discr <- housepi_smart_unwrap(discr, "discr")
      prior <- housepi_smart_unwrap(prior, "prior")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      proposal_scale <- housepi_smart_unwrap(proposal_scale, "proposal_scale")

      # load quosures
      value <- housepi_extractor(!!value)
      func <- housepi_extractor(!!func)
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)

      # load dots # UPDATE: not super elegant...
      dots <- list(...)
      # flatten if single quosures list (from enquos)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      args <- housepi_extractor(!!!dots)

      # flatten function arguments if quosures
      if (rlang::is_quosure(func_args)) {
        func_args <- housepi_extractor(!!func_args)
      } else {
        func_args <- housepi_extractor(!!!func_args)
      }

      # control valid name
      assert_valid_definition_name(name, self$family)
      if (all(is.null(c(value, func)), is.na(prior))) stop_h(
        sprintf("Incomplete definition for parameter '%s'.", name),
        "Provide at least one between 'value', 'prior', and 'func'."
      )
      if (!is.null(func) && !is.na(prior)) stop_h(
        sprintf("Wrong definition for parameter '%s'.", name),
        "Provide only one between 'func' and 'prior'."
      )
      if (!is.null(value) && !is.null(func)) stop_h(
        sprintf("Wrong definition for parameter '%s'.", name),
        "Provide only one between 'value' and 'func'."
      )

      opts <- list(name = as.character(name),
                   value = NA_real_,
                   const = as.logical(const),
                   discr = isTRUE(discr),
                   func = NA_character_,
                   prior = as.character(prior),
                   min = NA_real_,
                   max = NA_real_,
                   proposal = as.character(proposal),
                   proposal_scale = NA_real_)

      # drop previous definition eventually
      rmv_name <- self$remove(opts$name)
      if (.warn_over) template_overwrite_warn(rmv_name, self$family)
      else self$store_warn_over <- c(self$store_warn_over, rmv_name)

      opts <- helper_variable_core(opts, self$family, value, func, func_args)
      opts <- helper_variable_model(opts, self$family, args)
      opts <- helper_variable_move(opts, self$family, proposal_scale)

      if (!isTRUE(const)) {
        if (!is.null(min)) opts$min <- min
        if (!is.null(max)) opts$max <- max
      }

      # append definition(s)
      self$defs <- rbind(self$defs, data.frame(opts))
      # assert last definition is correct
      if (isTRUE(.assert)) housepi_assert_args(self, last = TRUE)

      invisible(self)
    },

    #' @inherit add_parameter
    add_quo_list = function(xquos, append = FALSE) {
      housepi_load_defin_quos(self, "value", xquos, append)
      housepi_assert_args(self)
      invisible(self)
    },

    #' @description
    #' Assert parameter definitions match the requirement list.
    #' @param required `character()`. Required definitions.
    #' @param capital `logical(1)`. Whether to raise error.
    #' @param warn_unused `logical(1)`. Warn if unused are found.
    assert = function(required, capital, warn_unused = TRUE) {
      defined <- self$defs[["name"]]
      # point out parameters that are not used
      if (isTRUE(warn_unused)) {
        template_not_used_warn(required, defined, self$family)
      }
      # stop if parameters are missing
      sms <- template_missing_stop(required, defined, self$family)
      if (isTRUE(capital) && length(sms) > 0) stop_h(sms)
      else sms
    }
  )
)
