## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Select augmented definitions of one kind, column `kind` dropped.
#' @param kind `character(1)` One of "latent", "member", "house".
#' @return `data.frame()` Definitions of that kind, `NULL` when none.
#' @noRd
augm_defs_of_kind <- function(defs, kind = definition_kinds_) {
  kind <- match.arg(kind)
  if (is.null(defs)) return(defs)
  out <- defs[defs[["kind"]] == kind, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out[["kind"]] <- NULL
  out
}

#' Hold the latent and observation definitions of a model.
#' @noRd
HousepiAugmented <- R6::R6Class(
  classname = "HousepiAugmented",
  cloneable = TRUE,

  ####  Private  ####

  private = list(
    deep_clone = function(name, value) helper_clone(value),

    # Kind of definitions managed, for display:
    # * "latent": latent variables;
    # * "member": member-level observations;
    # * "house": household-level observations.
    kind_ = "latent"
  ),

  ####  Active  ####

  active = list(
    # `character(1)` Family of the managed kind.
    family = function(value) kind2family(private$kind_),

    # `character()` Names of member-level observations.
    observation_names = function(value) {
      if (is.null(self$defs)) return(character())
      self$defs[["name"]][self$defs[["kind"]] == "member"]
    },

    # `character()` Observations producing a data column:
    # * household level -> always;
    # * member level -> when carrying no model.
    surrogate_names = function(value) {
      if (is.null(self$defs)) return(character())
      mask <- self$defs[["kind"]] == "house" |
        (self$defs[["kind"]] == "member" & is.na(self$defs[["model"]]))
      self$defs[["name"]][mask]
    },

    # `character()` Character `status =` labels across latent definitions,
    # in order of first appearance.
    status_labels = function(value) {
      if (is.null(self$defs)) return(character())
      labels <- unlist(self$defs[["status"]], use.names = FALSE)
      labels <- labels[!is.na(labels)]
      # literal character statuses only (integers are already codes)
      labels <- labels[is.na(suppressWarnings(as.integer(labels)))]
      unique(as.character(labels))
    }
  ),

  ####  Public  ####

  public = list(
    # `character()` Cached overwritten names.
    store_warn_over = character(),
    # `data.frame` Formatted latent definitions.
    defs = NULL,

    # Create an empty holder.
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    # Remove the definition matching `name` and `status`.
    # @return `character(1)` | `NULL` Name removed, `NULL` when none matched.
    remove = function(name, status) {
      mask <- self$defs[["name"]] == name
      if (any(mask)) {
        # same name and status, raw status an integer code or a label
        mask <- mask & vapply(self$defs[["status"]], function(x) {
          identical(as.character(x), as.character(status))
        }, logical(1))
        if (any(mask)) {
          self$defs <- self$defs[!mask, ]
          return(name)
        }
      }

      c()
    },

    # Whether each infected status in `required` has a definition block.
    # * generic definition, without status -> covers them all;
    # * household observation -> needs no block.
    # @param defs `data.frame()` Definitions with status codes resolved.
    # @param required `integer()` Infected statuses the model can carry:
    # every dictionary entry with a status process, else the observed ones.
    covers_statuses = function(defs, required) {
      if (is.null(defs)) return(TRUE)

      # household observations need no status block
      block <- defs[defs[["kind"]] != "house", , drop = FALSE]
      if (!nrow(block)) return(TRUE)
      declared <- unlist(block[["status"]], use.names = FALSE)
      if (anyNA(declared)) return(TRUE)

      required <- required[!is.na(required)]
      length(required) && all(required %in% as.integer(declared))
    },

    # Add one definition, arguments as in add_latent() and add_observation().
    add_quo = function(name, ...,
                       status = NA_integer_,
                       const = NULL, discr = NA,
                       model = NA_character_,
                       func = NULL,
                       min = NULL, max = NULL,
                       proposal = NA_character_,
                       level = NULL,
                       .warn_over = TRUE,
                       .assert = FALSE,
                       .kind = private$kind_) {
      # kind defaults to managed (see private$kind_, public$add_quo_list)
      .kind <- match.arg(.kind, definition_kinds_)
      # load non-quosures
      name <- housepi_smart_unwrap(name, "name")
      level <- housepi_smart_unwrap(level, "level", null_ok = TRUE)
      status <- housepi_smart_unwrap(status, "status")
      discr <- housepi_smart_unwrap(discr, "discr")
      model <- housepi_smart_unwrap(model, "model")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      assert_valid_law_name(model, "model")
      assert_valid_law_name(proposal, "proposal")

      # load quosures
      # an observation ignores const, see below
      if (identical(.kind, "latent")) assert_const_value_(const, "latent", name)
      const <- housepi_extractor(!!const)
      func <- helper_load_func_(func)
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)

      # extract args and proposal from the dots
      def_args <- helper_definition_args(
        proposal, list(...), kind2family(.kind), name
      )
      args <- def_args$args
      proposal <- def_args$proposal

      # a per-definition 'level' overrides the managed kind for this entry only
      if (!is.null(level) && !all(is.na(level))) {
        if (!.kind %in% c("member", "house")) stop_h(
          "{.arg level} is reserved to observations.",
          "Declare {.val {name}} with {.fn add_observation}."
        )
        .kind <- match.arg(level, c("member", "house"))
      }

      # status on a household observation -> error, one value per household
      if (identical(.kind, "house") && !all(is.na(status))) stop_h(
        "Household observation {.val {name}} cannot be status specific.",
        "Drop {.arg status} or declare {.val {name}} at member level."
      )

      if (!identical(.kind, "latent")) {
        # a household observation cannot be sampled
        if (identical(.kind, "house") && !is.na(model)) stop_h(
          "Household observation {.val {name}} cannot be sampled.",
          "Provide {.arg func} instead of {.arg model}."
        )
        # observation, never sampled -> sampling options dropped, dots included
        in_dots <- intersect(c("const", "proposal"), names(args))
        if (length(in_dots)) args <- args[!names(args) %in% in_dots]
        ignored <- c("const", "proposal")[c(!is.null(const), !is.na(proposal))]
        ignored <- union(ignored, in_dots)
        if (length(ignored)) warn_h(
          "Observation {.val {name}} cannot be sampled.",
          `!` = "Ignored {.arg {ignored}} definition."
        )
        proposal <- NA_character_
        # observation pinned to its datum, the like-named column
        const <- if (!is.na(model)) paste0(name, "|D")
      }

      # managed kind read from private$kind_, restored on exit
      old_kind <- private$kind_
      on.exit(private$kind_ <- old_kind, add = TRUE)
      private$kind_ <- .kind

      assert_valid_definition_name(name, self$family)
      # check status
      assert_valid_status(status)
      # check the definition is meaningful
      assert_no_sibling_dist_arg(args, self$family, name)
      # a numeric func is a constant
      if (is.null(const) && is.na(model) &&
            is.numeric(showcaser_extract_numeric(func))) {
        const <- func
        func <- NULL
      }
      assert_definition_core_(self$family, name, const, func, model, proposal,
                              bounded = !all(is.null(c(min, max))))

      # auxiliary typed options for Rcpp build method
      opts <- list(name = as.character(name),
                   status = as.character(status),
                   const = if (is.null(const)) NA_real_ else const,
                   discr = as.logical(discr),
                   func = if (is.null(func)) NA_character_ else func,
                   fcumul = FALSE, # @todo expose as an option
                   model = as.character(model),
                   kind = as.character(.kind),
                   min = if (is.null(min)) NA_real_ else min,
                   max = if (is.null(max)) NA_real_ else max,
                   proposal = proposal)

      # a name stands for one quantity of the model, whatever kind declares it
      assert_no_kind_clash(self$defs, opts$name, .kind)

      # drop any previous definition
      rmv_name <- self$remove(opts$name, opts$status)
      warn_or_store_overwrite_(self, rmv_name, .warn_over)

      opts <- helper_variable_args(opts, self$family, args)
      opts <- helper_variable_move(opts, self$family)
      opts <- helper_variable_level(opts, self$defs)

      # append definition(s)
      new_defs <- data.frame(opts)
      self$defs <- rbind(self$defs, new_defs)
      # assert appended definitions are correct
      if (isTRUE(.assert)) housepi_assert_args(self, last = nrow(new_defs))

      invisible(self)
    },

    # Add definitions from a list of quosures.
    # * not `append` -> previous definitions of the managed kind dropped.
    add_quo_list = function(xquos, append = FALSE, status = NA_integer_,
                            .kind = definition_kinds_) {
      managed <- match.arg(.kind)

      # managed kind read from private$kind_, restored on exit
      old_kind <- private$kind_
      on.exit(private$kind_ <- old_kind, add = TRUE)
      private$kind_ <- managed

      if (isFALSE(append)) {
        # kinds share one store and reset independently:
        # * observation levels -> one family, reset together.
        drop <- if (managed %in% c("member", "house")) c("member", "house")
        else managed
        if (is.null(self$defs)) self$initialize()
        else self$defs <- self$defs[!self$defs[["kind"]] %in% drop, ]
      }

      # a bulk-level 'status =' applies to every entry not setting its own
      defaults <- if (all(is.na(status))) list() else list(status = status)
      housepi_load_defin_quos(self, "func", xquos, append = TRUE,
                              defaults = defaults)
      housepi_assert_args(self)
      invisible(self)
    },

    # Assert the latents the model requires are defined.
    # * `wall` -> error on missing ones, else their messages returned;
    # * `warn_unused` -> warning on latents the model never reads.
    assert = function(required, wall, warn_unused) {
      defined <- self$defs[["name"]]

      # unused latents -> warning, observations excluded as never required
      if (isTRUE(warn_unused)) {
        non_obs <- defined[self$defs[["kind"]] == "latent"]
        definition_not_used_warn(required, non_obs, self$family)
      }
      # missing latents -> error, or messages returned
      sms <- definition_missing_stop(required, defined, self$family)
      wall_or_return_(sms, wall)
    }
  )
)
