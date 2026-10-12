## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

########  HousepiInfecProfile  ########

#' Hold the infectivity-profile definitions of a model.
#' @noRd
HousepiInfecProfile <- R6::R6Class(
  classname = "HousepiInfecProfile",
  cloneable = TRUE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value)
  ),

  ####  Public  ####

  public = list(
    # `character(1)` Holder family.
    family = "infectivity profile",
    # `logical(1)` Whether the model is discrete-time.
    discr_model = FALSE,
    # `character()` Cached overwritten names.
    store_warn_over = character(),
    # `data.frame` Formatted infectivity profile definitions.
    defs = NULL,

    # Create an empty holder.
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    # Remove the definition specific to `status`.
    # @return `integer(1)` | `NULL` Status removed, `NULL` when none matched.
    remove = function(status) {
      mask <- vapply(self$defs[["status"]], function(x) {
        identical(x, as.integer(status))
      }, logical(1))
      if (!any(mask)) return(NULL)
      self$defs <- self$defs[!mask, ]
      return(status)
    },

    # Add one definition, arguments as in add_infectivity_profile().
    add_quo = function(func, ...,
                       status = NA_integer_, discr = NA,
                       min = NULL, max = NULL,
                       .warn_over = TRUE, .assert = FALSE) {
      # load non-quosures
      func <- housepi_smart_unwrap(func, "func")
      status <- housepi_smart_unwrap(status, "status")
      discr <- housepi_smart_unwrap(discr, "discr")
      assert_valid_law_name(func, "function")
      # load quosures
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)
      # load dots: a stray 'name' of a listed definition is no argument
      dots <- flatten_definition_dots(list(...))
      func_args <- housepi_extractor(!!!dots[names(dots) != "name"])

      # check infectivity profile function
      if (is.na(func)) stop_h(
        "Incomplete definition: missing function name {.arg func}."
      )
      assert_valid_status(status)

      # drop any previous definition
      rmv_name <- self$remove(status)
      warn_or_store_overwrite_(self, rmv_name, .warn_over)

      # a discrete-time model requires a discrete profile
      if (isFALSE(discr) && self$discr_model) stop_h(
        "Infectivity profile of a discrete-time model must be discrete.",
        "Drop {.code discr = FALSE} or compile the model in continuous time."
      )
      if (is.na(discr)) discr <- self$discr_model

      # `min` not provided:
      # * discrete time, no `loc` -> 1;
      # * otherwise -> 0, `loc` already placing the origin.
      if (is.null(min))
        min <- if (discr && !"loc" %in% names(func_args)) "1|N" else "0|N"
      if (is.null(max)) max <- NA_character_

      # append new definition
      new_defs <- data.frame(list(
        func = housepi_combine_to_call(func, func_args),
        status = as.integer(status), discr = discr, min = min, max = max
      ))
      self$defs <- rbind(self$defs, new_defs)
      # assert appended definitions are correct
      if (isTRUE(.assert)) housepi_assert_args(self, last = nrow(new_defs))

      invisible(self)
    },

    # Add definitions from a list of quosures.
    # * not `append` -> previous definitions dropped.
    add_quo_list = function(xquos, append = FALSE) {
      housepi_load_defin_quos(self, "func", xquos, append)
      housepi_assert_args(self)
      invisible(self)
    },

    # Read the definitions in a readable form.
    # @param defs `data.frame()` Definitions to read,
    # the model passing the table with its sources resolved.
    showcase = function(defs = self$defs) {
      if (!length(defs)) return(NULL)
      output <- housepi_defin_showcaser(defs)

      # append func_args if not present
      i <- which(names(output) == "func")
      if (names(output)[i + 1] != "func_args")
        output <- append(output, list(func_args = NA), after = i)
      output
    }
  )
)
