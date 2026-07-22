## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

############################################################################
############################################################################
########################### Housepi Infec Profile ##########################
############################################################################
############################################################################

#' @title Internal R6 holder for housepi infectivity-profile definitions
#' @description
#' Helper class managing infectivity-profile definitions of the model.
#' @keywords internal
#' @noRd
HousepiInfecProfile <- R6::R6Class(
  classname = "HousepiInfecProfile",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Definitions required by infectivity profile definitions
    # according to the "type" of the definition
    find_required_def_ = function(type) {
      find_required_definitions(
        self$defs, type, c("func", "max"), self$family
      )
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field augm_required `character()`.
    #' Augmented variables required by infectivity profile definition.
    augm_required = function(value) {
      private$find_required_def_("A")
    },

    #' @field data_required `character()`.
    #' Dataset entries required by infectivity profile definition.
    data_required = function(value) {
      private$find_required_def_("D")
    },

    #' @field group_risks_required `character()`.
    #' Find group risks required by infectivity profile definition.
    group_risks_required = function(value) {
      private$find_required_def_("G")
    },

    #' @field param_required `character()`.
    #' Parameters required by infectivity profile definition.
    param_required = function(value) {
      private$find_required_def_("P")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field family `character(1)`. Holder family.
    family = "infectivity profile",
    #' @field store_warn_over `character()`. Cached overwritten names.
    store_warn_over = character(),
    #' @field defs `data.frame`. Formatted infectivity profile definitions.
    defs = NULL,

    #' @description
    #' New object of class
    #' [HousepiInfecProfile][housepi::HousepiInfecProfile].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description Remove a given infectivity profile definition.
    #' @param status `integer(1)`. Associated status of the definition.
    remove = function(status) {
      mask <- sapply(self$defs[["status"]], function(x) {
        identical(x, as.integer(status))
      })
      if (any(mask)) {
        self$defs <- self$defs[!mask, ]
        return(status)
      }
      c()
    },

    #' @inherit add_infectivity_profile
    add_quo = function(func, ...,
                       status = NA_integer_, discr = FALSE,
                       min = NULL, max = NULL,
                       .warn_over = TRUE, .assert = FALSE) {
      # load non-quosures
      func <- housepi_smart_unwrap(func, "func")
      status <- housepi_smart_unwrap(status, "status")
      discr <- housepi_smart_unwrap(discr, "discr")
      # load quosures
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)
      # load dots # UPDATE: not super elegant...
      dots <- list(...)
      # flatten if single quosures list (from enquos)
      if (length(dots) == 1 && rlang::is_quosures(dots[[1]])) dots <- dots[[1]]
      func_args <- housepi_extractor(!!!dots)
      # ensure the argument name doesn't mistakenly be interpreted as argument
      func_args <- func_args[names(func_args) != "name"]

      # check infectivity profile function
      if (is.na(func)) stop_h(
        "Incomplete definition: missing function name {.arg func}."
      )
      assert_valid_status(status)

      # drop previous definition eventually
      rmv_name <- self$remove(status)
      if (.warn_over) template_overwrite_warn(rmv_name, self$family)
      else self$store_warn_over <- c(self$store_warn_over, rmv_name)

      # append last definition
      self$defs <- rbind(self$defs, data.frame(list(
        func = housepi_combine_to_call(func, func_args),
        status = as.integer(status),
        min = ifelse(is.null(min), ifelse(isTRUE(discr), "1|N", "0|N"), min),
        max = ifelse(is.null(max), NA_character_, max),
        discr = isTRUE(discr)
      )))
      # assert last definition is correct
      if (isTRUE(.assert)) housepi_assert_args(self, last = TRUE)

      invisible(self)
    },

    #' @inherit add_infectivity_profile
    add_quo_list = function(xquos, append = FALSE) {
      housepi_load_defin_quos(self, "func", xquos, append)
      housepi_assert_args(self)
      invisible(self)
    },

    #' @description
    #' Returns a readable definition object.
    showcase = function() {
      if (length(self$defs) > 0) {
        output <- housepi_defin_showcaser(self$defs)

        # append func_args if not present
        i <- which(names(output) == "func")
        if (names(output)[i + 1] != "func_args") {
          output <- append(output, list(func_args = NA), after = i)
        }
        return(output)
      }
      NULL
    }
  )
)
