## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @description Augmented definitions of a single kind.
#' Column 'kind' is dropped.
#' @param defs `data.frame()` Definitions to filter.
#' @param kind `character(1)` One of "latent", "member", "house", "status".
#' @details Used to read a resolved copy of the declarations.
#' @return `data.frame()` or `NULL` when no definition of that kind exists.
#' @keywords internal
#' @noRd
augm_defs_of_kind <- function(defs, kind = definition_kinds_) {
  kind <- match.arg(kind)
  if (is.null(defs)) return(defs)
  out <- defs[defs[["kind"]] == kind, , drop = FALSE]
  if (!nrow(out)) return(NULL)
  out[["kind"]] <- NULL
  out
}

#' @title Internal R6 holder for housepi latent definitions
#' @description
#' Helper class managing latent definitions of the model object.
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

    # Kind of definitions currently being managed. Used for display.
    # One between:
    # "latent": latent variables;
    # "member": member-level observation;
    # "house": household-level observations;
    # "status": status process (resolves status in simulation).
    kind_ = "latent"
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field family `character(1)`. Holder family currently being managed.
    family = function(value) kind2family(private$kind_),

    #' @field observation_names `character()`.
    #' Names of member-level observations.
    observation_names = function(value) {
      if (is.null(self$defs)) return(character())
      self$defs[["name"]][self$defs[["kind"]] == "member"]
    },

    #' @field surrogate_names `character()`.
    #' Names of deterministic observations producing a data column:
    #' member- and household-level observations, carrying no model.
    surrogate_names = function(value) {
      if (is.null(self$defs)) return(character())
      mask <- self$defs[["kind"]] == "house" |
        (self$defs[["kind"]] == "member" & is.na(self$defs[["model"]]))
      self$defs[["name"]][mask]
    },

    #' @field status_names `character()`.
    #' Reserved node name of the status process, if defined.
    status_names = function(value) {
      if (is.null(self$defs)) return(character())
      self$defs[["name"]][self$defs[["kind"]] == "status"]
    },

    #' @field status_dictionary `character()`.
    #' Ordered infected-status labels (code `k` is the `k`-th entry). Taken from
    #' the status process when present, else factorized from literal character
    #' `status =` labels declared across the latent definitions.
    status_dictionary = function(value) {
      if (length(self$status_levels)) return(self$status_levels)
      if (is.null(self$defs)) return(character())
      labels <- unlist(self$defs[["status"]], use.names = FALSE)
      labels <- labels[!is.na(labels)]
      # literal character statuses only (integers are already codes)
      labels <- labels[is.na(suppressWarnings(as.integer(labels)))]
      unique(as.character(labels))
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @field store_warn_over `character()`. Cached overwritten names.
    store_warn_over = character(),
    #' @field defs `data.frame`. Formatted latent definitions.
    defs = NULL,
    #' @field status_levels `character()`. Ordered status-process labels.
    status_levels = character(),

    #' @description
    #' New object of class
    #' [HousepiAugmented][housepi::HousepiAugmented].
    initialize = function() {
      self$defs <- NULL
      invisible(self)
    },

    #' @description Record the ordered status-process labels (the dictionary).
    #' @param levels `character()`. Status labels in declaration order,
    #' one per status, empty when unlabelled.
    set_status_levels = function(levels) {
      # the statuses are coded 1...K, labelled or not
      # K is capped by the engine status storage
      if (length(levels) > max_status_level) stop_h(
        paste("Status process implies {length(levels)} statuses, exceeding",
              "the maximum of {max_status_level}.")
      )
      levels <- levels[!is.na(levels) & nzchar(levels)]
      self$status_levels <- as.character(levels)
      invisible(self)
    },

    #' @description Remove a given latent definition.
    #' @param name `character(1)`. Latent definitions to remove.
    #' @param status `integer(1)`. Associated status of definition.
    remove = function(name, status) {
      mask <- self$defs[["name"]] == name
      if (any(mask)) {
        # check if same name and status (raw: integer code or character label)
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

    #' @description Whether every infected status lands on a definition block.
    #' A generic (unstatused) definition covers them all, otherwise each status
    #' the model can carry needs one of its own: every dictionary entry when a
    #' status process draws them, the observed ones when data supplies them.
    #' @param defs `data.frame()`. Definitions with status codes resolved.
    #' @param status_set `integer()`. Infected statuses observed in data.
    #' @return `logical(1)`
    covers_statuses = function(defs, status_set) {
      if (is.null(defs)) return(TRUE)

      # household observations and the status process need no status block
      block <- defs[!defs[["kind"]] %in% c("house", "status"), , drop = FALSE]
      if (!nrow(block)) return(TRUE)
      declared <- unlist(block[["status"]], use.names = FALSE)
      if (anyNA(declared)) return(TRUE)

      # a drawn status is any dictionary code, whatever the data carries
      required <- if (length(self$status_names)) {
        seq_along(self$status_dictionary)
      } else {
        status_set[!is.na(status_set)]
      }
      length(required) && all(required %in% as.integer(declared))
    },

    #' @inherit add_latent
    add_quo = function(name, ...,
                       value = NULL, status = NA_integer_,
                       const = NA, discr = NA,
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
      const <- housepi_smart_unwrap(const, "const")
      discr <- housepi_smart_unwrap(discr, "discr")
      model <- housepi_smart_unwrap(model, "model")
      proposal <- housepi_smart_unwrap(proposal, "proposal")
      assert_valid_law_name(model, "model")
      assert_valid_law_name(proposal, "proposal")

      # load quosures
      value <- housepi_extractor(!!value)
      func <- housepi_extractor(!!func)
      min <- housepi_extractor(!!min)
      max <- housepi_extractor(!!max)

      # extract args and proposal from the dots
      def_args <- helper_definition_args(
        proposal, list(...), kind2family(.kind),
        if (identical(.kind, "status")) NULL else name
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

      # a household value is one per household: no status block selects it
      if (identical(.kind, "house") && !all(is.na(status))) stop_h(
        "Household observation {.val {name}} cannot be status specific.",
        "Drop {.arg status} or declare {.val {name}} at member level."
      )

      if (!identical(.kind, "latent")) {
        # a household observation cannot be sampled
        if (identical(.kind, "house") && !is.na(model)) stop_h(
          "Household observation {.val {name}} cannot be sampled.",
          "Provide {.arg value} instead of {.arg model}."
        )
        # an observation is evaluated at its datum, never sampled:
        # sampling options are dropped, those landing in the dots included
        in_dots <- intersect(c("const", "proposal"), names(args))
        if (length(in_dots)) args <- args[!names(args) %in% in_dots]
        ignored <- c("const", "proposal")[!is.na(c(const, proposal))]
        ignored <- union(ignored, in_dots)
        if (length(ignored)) warn_h(
          "Observation {.val {name}} cannot be sampled.",
          `!` = "Ignored {.arg {ignored}} definition."
        )
        const <- NA
        proposal <- NA_character_

        if (!is.na(model)) {
          const <- TRUE
          # observation reads its datum from a like-named column
          # status process requires no data column
          if (is.null(value)) {
            value <- if (identical(.kind, "status")) "0|N"
            else paste0(name, "|D")
          }
        }
      }

      # NOTE: single source of truth for currently managed kind
      old_kind <- private$kind_
      on.exit(private$kind_ <- old_kind, add = TRUE)
      private$kind_ <- .kind

      # add_status_process names its node, no user name may claim the prefix
      assert_valid_definition_name(
        name, self$family, allow_reserved = identical(.kind, "status")
      )
      # control valid status
      assert_valid_status(status)
      # control meaningful definition
      assert_no_sibling_dist_arg(args, self$family, name)
      # a latent without a model is sampled within its bounds
      is_bounded <- identical(.kind, "latent") && !all(is.null(c(min, max)))
      if (all(is.null(c(value, func)), is.na(model), !is_bounded)) stop_h(
        "Incomplete definition for latent {.val {name}}.",
        "Provide at least one between {.arg {c('value', 'model', 'func')}}."
      )
      if (is_bounded && !is.null(value) && is.na(model) &&
            is.na(proposal) && is.na(const)) stop_h(
        "Latent {.val {name}} has value and bounds, but no model or proposal.",
        "Add {.arg proposal} to sample it, or drop {.arg min} and {.arg max}."
      )
      if (!is.null(func) && !is.na(model)) stop_h(
        "Wrong definition for latent {.val {name}}.",
        "Provide only one between {.arg {c('func', 'model')}}."
      )
      if (!is.null(value) && !is.null(func)) stop_h(
        "Wrong definition for latent {.val {name}}.",
        "Provide only one between {.arg {c('value', 'func')}}."
      )

      # auxiliary typed options for Rcpp build method
      opts <- list(name = as.character(name),
                   value = NA_real_,
                   status = as.character(status),
                   const = as.logical(const),
                   discr = as.logical(discr),
                   func = NA_character_,
                   fcumul = FALSE, # UPDATE: can I make this an option?
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

      opts <- helper_variable_core(opts, self$family, value, func)
      opts <- helper_variable_args(opts, self$family, args, func)
      opts <- helper_variable_move(opts, self$family)
      opts <- helper_variable_level(opts, self$defs)

      # append definition(s)
      new_defs <- data.frame(opts)
      self$defs <- rbind(self$defs, new_defs)
      # assert appended definitions are correct
      if (isTRUE(.assert)) housepi_assert_args(self, last = nrow(new_defs))

      invisible(self)
    },

    #' @inherit add_latent
    add_quo_list = function(xquos, append = FALSE, status = NA_integer_,
                            .kind = definition_kinds_) {
      managed <- match.arg(.kind)

      # NOTE: single source of truth for currently managed kind
      old_kind <- private$kind_
      on.exit(private$kind_ <- old_kind, add = TRUE)
      private$kind_ <- managed

      if (isFALSE(append)) {
        # NOTE: kinds share one store, but are treated as independent;
        # observations are a single family, so both levels reset together
        drop <- if (managed %in% c("member", "house")) c("member", "house")
        else managed
        if (is.null(self$defs)) self$initialize()
        else self$defs <- self$defs[!self$defs[["kind"]] %in% drop, ]
      }

      # a bulk-level 'status =' applies to every entry not setting its own
      defaults <- if (all(is.na(status))) list() else list(status = status)
      housepi_load_defin_quos(self, "value", xquos, append = TRUE,
                              defaults = defaults)
      housepi_assert_args(self)
      invisible(self)
    },

    #' @description
    #' Assert latent definitions from the model implementation script.
    #' @param required `character()`. Required definitions.
    #' @param wall `logical(1)`. Whether to raise error.
    #' @param warn_unused `logical(1)`. Warn if unused are found.
    assert = function(required, wall, warn_unused) {
      defined <- self$defs[["name"]]

      # point out latent variables that are not used
      # NOTE: observations are never model-required, so excluded from check
      if (isTRUE(warn_unused)) {
        non_obs <- defined[self$defs[["kind"]] == "latent"]
        template_not_used_warn(required, non_obs, self$family)
      }
      # stop if latent variables are missing
      sms <- template_missing_stop(required, defined, self$family)
      wall_or_return_(sms, wall)
    }
  )
)
