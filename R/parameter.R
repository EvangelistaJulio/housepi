## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Hold the parameter definitions of a model.
#' @noRd
HousepiParameter <- R6::R6Class(
  classname = "HousepiParameter",
  cloneable = TRUE,

  ####  Public  ####

  public = list(
    # `character(1)` Holder family.
    family = "parameter",

    # `character()` Cached overwritten names.
    store_warn_over = character(),

    # `data.frame` Formatted parameter definitions.
    defs = NULL,

    # Create an empty holder.
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    # Remove the parameter `name`.
    # @return `character(1)` | `NULL` Name removed, `NULL` when none matched.
    remove = function(name) {
      mask <- self$defs[["name"]] == name
      if (any(mask)) {
        self$defs <- self$defs[!mask, ]
        return(name)
      }
      c()
    },

    # Add one definition, arguments as in add_parameter().
    add_quo = function(name, ...,
                       const = NULL, discr = NA,
                       prior = NA_character_,
                       func = NULL,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       .warn_over = TRUE,
                       .assert = FALSE) {
      # load non-quosures
      name <- housepi_smart_unwrap(name, "name")
      discr <- housepi_smart_unwrap(discr, "discr")
      prior <- housepi_smart_unwrap(prior, "prior")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      assert_valid_law_name(prior, "prior")
      assert_valid_law_name(proposal, "proposal")

      # load quosures
      assert_const_value_(const, self$family, name)
      const <- housepi_extractor(!!const)
      func <- helper_load_func_(func)
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)

      # extract args and proposal from the dots
      def_args <- helper_definition_args(proposal, list(...), self$family, name)
      args <- def_args$args
      proposal <- def_args$proposal

      # check the name
      assert_valid_definition_name(name, self$family)
      assert_no_sibling_dist_arg(args, self$family, name)
      # a numeric func is a constant
      if (is.null(const) && is.na(prior) &&
            is.numeric(showcaser_extract_numeric(func))) {
        const <- func
        func <- NULL
      }
      assert_definition_core_(self$family, name, const, func, prior, proposal,
                              bounded = !all(is.null(c(min, max))))

      opts <- list(name = as.character(name),
                   const = if (is.null(const)) NA_real_ else const,
                   discr = as.logical(discr),
                   func = if (is.null(func)) NA_character_ else func,
                   prior = as.character(prior),
                   min = if (is.null(min)) NA_real_ else min,
                   max = if (is.null(max)) NA_real_ else max,
                   proposal = proposal)

      # drop any previous definition
      rmv_name <- self$remove(opts$name)
      warn_or_store_overwrite_(self, rmv_name, .warn_over)

      opts <- helper_variable_args(opts, self$family, args)
      opts <- helper_variable_move(opts, self$family)

      # append definition(s)
      new_defs <- data.frame(opts)
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

    # Assert the parameters the model requires are defined.
    # * `wall` -> error on missing ones, else their messages returned;
    # * `warn_unused` -> warning on parameters the model never reads.
    assert = function(required, wall, warn_unused = TRUE) {
      defined <- self$defs[["name"]]
      # unused parameters -> warning
      if (isTRUE(warn_unused))
        definition_not_used_warn(required, defined, self$family)

      # missing parameters -> error, or messages returned
      sms <- definition_missing_stop(required, defined, self$family)
      wall_or_return_(sms, wall)
    }
  )
)
