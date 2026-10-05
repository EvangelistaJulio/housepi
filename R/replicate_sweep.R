## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

## Replicate sweep shared by the checks that simulate surveys and fit them
## back: one simulation call per household count, one fit per replicate,
## each fit reduced to one row per parameter.

## Roxygen text shared by the checks, read at the use site as `r <name>`.

rd_sweep_nrep_ <- "`integer(1)` Number of surveys simulated and fitted back
at each household count."
rd_sweep_sampler_ <- "`function` Sampler of the household structure the surveys
are generated on, called as `structure_sampler(data_wide, houses)`.
`NULL` reuses the attached structure verbatim and ignores `houses`."
rd_sweep_seed_ <- "Replicates are independent draws: call `set.seed()` before
the check to replay a whole sweep, and only at `threads = 1`."

#' Columns the surveys are generated on: the attached data without the
#' observations and latent variables, which simulation regenerates and a
#' definition reading one of them must not see observed.
#' Status is kept as simulation receives it: a transmission model draws it anew.
#' @noRd
sweep_base_structure_ <- function(hhmodel) {
  data <- hhmodel$data_wide
  if (!NROW(data)) stop_h(
    "No data is attached to the model.",
    paste("The check generates its surveys on the structure of the",
          "attached data: add one with {.fn add_data_wide}.")
  )
  written <- c(hhmodel$observations[["name"]], hhmodel$latent[["name"]])
  data[setdiff(names(data), written)]
}

#' Simulation value of every parameter, `NA` where it has none.
#' @noRd
sweep_pins_ <- function(hhmodel) {
  defs <- get_parameters(hhmodel)
  stats::setNames(defs[["value"]], defs[["name"]])
}

#' Simulation table: the parameters of `pins` holding a value, one row per
#' replicate. A parameter without a column is drawn from its prior.
#' @noRd
sweep_param_table_ <- function(pins, nrep) {
  pins <- pins[!is.na(pins)]
  table <- data.frame(row.names = seq_len(nrep))
  table[names(pins)] <- as.list(pins)
  table
}

#' The parameters a check is about: estimated, hence carrying a prior and
#' free.
#' @noRd
sweep_estimated_ <- function(hhmodel) {
  defs <- get_parameters(hhmodel)
  estimated <- defs[["name"]][!is.na(defs[["prior"]]) &
                                !defs[["const"]] %in% TRUE]
  if (!length(estimated)) stop_h(
    "The model holds no parameter to estimate.",
    "The check needs a parameter carrying a prior and free to move."
  )
  estimated
}

#' Diagnostics and verdict the fit reached, one row per parameter of
#' `variables`.
#' @noRd
sweep_diagnostics_ <- function(draws, variables) {
  diag <- draws$.__enclos_env__$private$diagnostics_
  diag <- diag[match(variables, diag[["variable"]]), ]
  rownames(diag) <- NULL
  diag[c("rhat", "ess_bulk", "flagged")]
}

#' Sweep cost per household count, averaged over the replicates.
#' @param replicates `data.frame` One row per replicate.
#' @noRd
sweep_overview_ <- function(replicates) {
  sweep_by_(replicates, "houses", function(d) {
    ok <- !is.na(d[["inference_secs"]]) & !is.na(d[["min_ess_bulk"]])
    data.frame(
      houses = d[["houses"]][[1L]], nrep = nrow(d),
      n_crashed = sum(!is.na(d[["error"]])),
      n_flagged = sum(d[["flagged"]] %in% TRUE),
      mean_inference_secs = mean(d[["inference_secs"]], na.rm = TRUE),
      secs_per_ess =
        sum(d[["inference_secs"]][ok]) / sum(d[["min_ess_bulk"]][ok])
    )
  })
}

#' One row of `fun` per group of `cols`, stacked back.
#' @noRd
sweep_by_ <- function(x, cols, fun) {
  sweep_stack_(lapply(split(x, row_keys_(x, cols)), fun), cols)
}

#' Bind rows back into one frame, ordered by the grouping columns.
#' @noRd
sweep_stack_ <- function(parts, cols) {
  out <- do.call(rbind, parts)
  out <- out[do.call(order, out[cols]), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Run `fn` on the model re-attached to `data`, with the given arguments.
#' @noRd
sweep_run_fn_ <- function(fn, hhmodel, data, args) {
  do.call(fn, c(list(add_data_wide(hhmodel, data)), args))
}

#' Household counts of the sweep, one cell each.
#' @param one `logical(1)` Demand a single count.
#' @noRd
sweep_houses_ <- function(structure, houses, verbatim, one = FALSE) {
  observed <- length(unique(structure[["hhid"]]))
  if (!verbatim) return(unique(sweep_house_counts_(houses, observed, one)))
  if (!all(is.na(houses))) warn_h(
    "Argument {.arg houses} is ignored.",
    "A {.code NULL} sampler reuses the attached structure verbatim."
  )
  observed
}

#' Household counts as whole positive numbers, `NA` standing for `observed`.
#' @param one `logical(1)` Demand a single count.
#' @noRd
sweep_house_counts_ <- function(houses, observed, one = FALSE) {
  n <- suppressWarnings(as.integer(houses))
  n[is.na(houses)] <- observed
  ok <- length(n) && !anyNA(n) && all(n >= 1L) && !(one && length(n) > 1L)
  if (!ok) stop_h("Number of households should be a positive integer.")
  n
}

#' Progress bar over the replicate grid, or nothing when it is switched off.
#' The bar of each fit is suppressed, so a sweep prints one bar.
#' @noRd
sweep_progress_bar_ <- function(progress, total, label) {
  if (!isTRUE(as.logical(progress[[1L]]))) return(NULL)
  cli::cli_progress_bar(label, total = total, .envir = parent.frame())
}

#' Distinct warning messages raised over a sweep, with their counts.
#' @param raised `list` Per-fit warnings, as a fit's `warnings` field holds.
#' @noRd
sweep_warning_counts_ <- function(raised) {
  held <- do.call(rbind, raised)
  if (!NROW(held)) {
    return(data.frame(message = character(), count = integer()))
  }
  counts <- tapply(held[["count"]], held[["message"]], sum)
  data.frame(message = names(counts), count = as.integer(counts),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' One warning per distinct message, with the count of replicates that
#' raised it. Foreign text carries no markup.
#' @param what `character(1)` Name of the check, heading the warning.
#' @noRd
sweep_report_warnings_ <- function(warnings, what) {
  if (!NROW(warnings)) return()
  # nolint next: object_usage_linter. read by the cli glue below
  n_warn <- nrow(warnings)
  lines <- mapply(\(n, msg) inline_cli_("{n} replicate{?s}: {msg}"),
                  warnings[["count"]], warnings[["message"]])
  header <- inline_cli_("{what} raised {n_warn} distinct warning{?s}.")
  bullets <- stats::setNames(lines, rep("!", length(lines)))
  emit_cli_condition(c(list(header), as.list(bullets)), "warn")
}

#' The sweep: one simulation call per household count, then one fit per
#' replicate, reduced by `reducer(draws, truth)` to one row per parameter
#' carrying `ess_bulk` and `flagged`. The fitted model is built from the given
#' model and never from the one that simulated, which holds the truth of the
#' last replicate.
#' @param sim_args,fit_args `list` Arguments of the simulation and of each fit.
#' @param label `character(1)` Label of the progress bar.
#' @noRd
sweep_run_ <- function(hhmodel, structure, sampler, houses, nrep, sim_args,
                       fit_args, reducer, progress, label) {
  estimates <- replicates <- raised <- list()
  bar <- sweep_progress_bar_(progress, length(houses) * nrep, label)
  for (n_houses in houses) {
    struct <- structure
    if (!is.null(sampler)) struct <- sampler(structure, n_houses)
    surveys <- sweep_run_fn_(simulate_with_parameters_, hhmodel, struct,
                             sim_args)
    cell <- lapply(seq_len(nrep), function(rep) {
      truth <- unlist(surveys$parameters[rep, , drop = FALSE])
      one <- sweep_one_replicate_(hhmodel, surveys$wide[[rep]], truth,
                                  n_houses, rep, fit_args, reducer)
      if (!is.null(bar)) cli::cli_progress_update(id = bar)
      one
    })
    rm(surveys)
    done <- Filter(Negate(is.null), lapply(cell, `[[`, "estimates"))
    # an engine message carries no markup: interpolate, never inline
    # nolint next: object_usage_linter. read by the cli glue below
    reason <- cell[[1L]]$replicate[["error"]]
    if (!length(done)) stop_h(
      "Every replicate crashed at {.val {n_houses}} households.",
      `!` = "{reason}"
    )
    estimates <- c(estimates, done)
    replicates <- c(replicates, lapply(cell, `[[`, "replicate"))
    raised <- c(raised, lapply(cell, `[[`, "warnings"))
  }
  list(estimates = sweep_stack_(estimates, c("houses", "rep")),
       replicates = sweep_stack_(replicates, c("houses", "rep")),
       warnings = sweep_warning_counts_(raised))
}

#' One replicate: the fit of the survey simulated at `truth`.
#' A replicate that crashes is recorded and the sweep carries on.
#' @noRd
sweep_one_replicate_ <- function(hhmodel, simulated, truth, houses, rep,
                                 fit_args, reducer) {
  crash <- NA_character_
  fit <- tryCatch(
    sweep_fit_one_(hhmodel, simulated, truth, houses, rep, fit_args, reducer),
    error = function(e) {
      crash <<- conditionMessage(e)
      NULL
    }
  )
  list(
    estimates = fit$estimates,
    warnings = fit$warnings,
    replicate = data.frame(
      houses = houses, rep = rep,
      inference_secs = if (is.null(fit)) NA_real_ else fit$secs,
      min_ess_bulk = if (is.null(fit)) NA_real_ else
        min(fit$estimates$ess_bulk),
      flagged = if (is.null(fit)) NA else any(fit$estimates$flagged),
      error = crash
    )
  )
}

#' Fit one survey back and reduce its draws.
#' Kept whole and free of sweep state, so it can be mapped.
#' @noRd
sweep_fit_one_ <- function(hhmodel, simulated, truth, houses, rep, fit_args,
                           reducer) {
  latent <- intersect(hhmodel$latent[["name"]], names(simulated))
  hhdata <- simulated[setdiff(names(simulated), latent)]
  draws <- sweep_run_fn_(run_inference, hhmodel, hhdata, fit_args)
  list(estimates = cbind(houses = houses, rep = rep, reducer(draws, truth)),
       warnings = draws$warnings,
       secs = draws$inference_secs)
}
