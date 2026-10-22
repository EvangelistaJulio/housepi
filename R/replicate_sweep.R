## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Replicate sweep of the checks that simulate surveys and fit them back:
# one simulation call per household count, one fit per replicate,
# each fit reduced to one row per parameter.

# Roxygen text shared by the checks, read at the use site as `r <name>`.

rd_sweep_nrep_ <- "`integer(1)`
Number of surveys simulated and fitted back at each household count."
rd_sweep_sampler_ <- "`function`
Sampler of the household structure the surveys are generated on,
called as `structure_sampler(data_wide, houses)`.
If `NULL`, the attached structure is reused verbatim and `houses` ignored."
rd_sweep_seed_ <- "Replicates are independent draws.
Call `set.seed()` before the check to replay a whole sweep,
at any `threads` up to `nrep`."
rd_sweep_reporting_ <- "Each survey goes through the `wide` reporting process
of the model before it is fitted back, see [add_reporting_process()].
Inference does not model it,
so a process losing information makes the check measure the misspecification."
rd_sweep_threads_ <- "`integer(1)` Number of threads.
Up to `threads` replicates are fitted at once.
Beyond `nrep`, the households of each fit share the threads left over."

#' Columns the surveys are generated on: the attached data without
#' the observations and latent variables, which simulation regenerates,
#' and which a definition reading one of them must not see observed.
#' Status is kept as simulation receives it, drawn anew by transmission.
#' @noRd
sweep_base_structure_ <- function(hh_model) {
  data <- get_data_wide(hh_model)
  if (!NROW(data)) stop_h(
    "No data is attached to the model.",
    paste("The check generates its surveys on the structure of the",
          "attached data: add one with {.fn add_data_wide}.")
  )
  written <- c(get_observations(hh_model)[["name"]],
               get_latent(hh_model)[["name"]])
  data[setdiff(names(data), written)]
}

#' Simulation table: the parameters of `pins` holding a value,
#' one row per replicate.
#' * parameter without a column -> drawn from its prior.
#' @noRd
sweep_param_table_ <- function(pins, nrep) {
  pins <- pins[!is.na(pins)]
  table <- data.frame(row.names = seq_len(nrep))
  table[names(pins)] <- as.list(pins)
  table
}

#' Parameters a check is about: the estimated ones, with a prior.
#' @noRd
sweep_estimated_ <- function(hh_model) {
  defs <- get_parameters(hh_model)
  estimated <- defs[["name"]][!is.na(defs[["prior"]])]
  if (!length(estimated)) stop_h(
    "The model holds no parameter to estimate.",
    "The check needs a parameter carrying a prior and free to move."
  )
  estimated
}

#' Diagnostics and verdict of the fit, one row per parameter of `variables`.
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
    ok <- !is.na(d[["min_ess_bulk"]])
    data.frame(
      houses = d[["houses"]][[1L]], nrep = nrow(d),
      n_flagged = sum(d[["flagged"]]),
      mean_inference_secs = mean(d[["inference_secs"]]),
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

#' Household counts of the sweep, one cell each.
#' @param one `logical(1)` Demand a single count.
#' @noRd
sweep_houses_ <- function(structure, houses, verbatim, one = FALSE) {
  observed <- length(unique(structure[["hh_id"]]))
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

#' Distinct warning messages raised over a sweep, with their counts.
#' @param raised `list` Per-fit warnings, as a fit's `warnings` field holds.
#' @noRd
sweep_warning_counts_ <- function(raised) {
  held <- do.call(rbind, raised)
  if (!NROW(held))
    return(data.frame(message = character(), count = integer()))

  counts <- tapply(held[["count"]], held[["message"]], sum)
  data.frame(message = names(counts), count = as.integer(counts),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' One warning per distinct message, with the count of replicates raising it.
#' Foreign text carries no markup.
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

#' Workers of a sweep: `outer` replicates fitted at once,
#' each fit's households sharing `inner` threads.
#' @noRd
sweep_split_ <- function(nrep, threads) {
  threads <- assert_count_(threads, "Number of threads")
  # @todo threads idle on the remainder of `threads %/% outer`,
  # and through a last wave holding fewer than `outer` replicates
  outer <- min(nrep, threads)
  list(outer = outer, inner = threads %/% outer)
}

#' Value of `expr`.
#' * `expr` fails -> error naming replicate `rep` at `houses` households,
#'   followed by the lines of the error raised.
#' @noRd
sweep_guard_ <- function(expr, rep, houses) {
  tryCatch(expr, error = function(e) {
    # rendered text: shown verbatim, never parsed as markup
    header <- paste(rlang::cnd_header(e), collapse = "\n")
    cause <- escape_cli_(c(`!` = header, rlang::cnd_body(e)))
    stop_h("Replicate {rep} at {.val {houses}} households failed.", cause)
  })
}

#' Run the sweep: one simulation call per household count,
#' then one fit per replicate, reduced by `reducer(draws, truth)`
#' to one row per parameter carrying `ess_bulk` and `flagged`.
#' The fitted model is built from the given model,
#' never from the simulating one, holding the truth of the last replicate.
#' Fits run in waves of replicates sampled at once,
#' each wave reduced before the next one runs.
#' Sampler settings are asserted before any survey is simulated,
#' and the first replicate that fails stops the sweep.
#' @param pins `numeric()` Simulation value of the parameters,
#' `NA` drawn from the prior.
#' @param fit_args `list` Sampler settings of each fit.
#' @param threads `integer(1)` Threads shared by the fits of a wave.
#' @param dots `list` Advanced setup of the check,
#' split between the simulation and the fits.
#' @param label `character(1)` Label of the progress bar.
#' @noRd
sweep_run_ <- function(hh_model, structure, sampler, houses, nrep, pins,
                       fit_args, threads, dots, reducer, progress, label) {
  adapt <- dots[names(dots) == "adapt_window"]
  do.call(resolve_setup_, c(fit_args, adapt))
  # fits read the wide survey only
  hh_model <- add_reporting_process(hh_model, long = NULL)
  # one simulation call per household count yields the whole cell,
  # the pins repeated one row per replicate
  sim_args <- c(list(draws = sweep_param_table_(pins, nrep)),
                dots[names(dots) == "time_step"])
  workers <- sweep_split_(nrep, threads)
  fit_args <- c(fit_args, adapt, threads = workers$inner, id_draws = ".none")
  waves <- split(seq_len(nrep), (seq_len(nrep) - 1L) %/% workers$outer)
  # the replicates at `n` households, reduced
  cell <- function(n, bar) {
    struct <- if (is.null(sampler)) structure else sampler(structure, n)
    surveys <- do.call(simulate_with_parameters_,
                       c(list(add_data_wide(hh_model, struct)), sim_args))
    truths <- lapply(seq_len(nrep),
                     \(rep) unlist(surveys$parameters[rep, , drop = FALSE]))
    # seeds drawn ahead of the fits: the reducers draw from the same stream
    # whatever the size of a wave
    seeds <- draw_seed_(nrep)
    unlist(lapply(waves, function(reps) {
      draws <- sweep_fit_wave_(hh_model, surveys$wide[reps], reps, n,
                               seeds[reps], fit_args, workers$outer, bar)
      Map(sweep_reduce_one_, draws, truths[reps], reps, seeds[reps],
          MoreArgs = list(houses = n, reducer = reducer))
    }), recursive = FALSE)
  }
  start <- proc.time()[["elapsed"]]
  ran <- with_progress_(progress, label, length(houses) * nrep, function(bar) {
    unlist(lapply(houses, cell, bar = bar), recursive = FALSE)
  })
  field <- \(k) lapply(ran, `[[`, k)
  list(estimates = sweep_stack_(field("estimates"), c("houses", "rep")),
       replicates = sweep_stack_(field("replicate"), c("houses", "rep")),
       warnings = sweep_warning_counts_(field("warnings")),
       secs = proc.time()[["elapsed"]] - start)
}

#' Run one wave: the surveys of `wides` built in order,
#' sampled `outer` at a time, housepi warnings held back.
#' * first survey failing -> stops every fit of the wave.
#' @param reps `integer()` Replicate of each survey.
#' @param houses `integer(1)` Household count of the surveys.
#' @param seeds `integer()` Engine seed of each survey.
#' @param bar Progress bar each fit fills by one unit, `NULL` for none.
#' @return Posterior draws of each survey.
#' @noRd
sweep_fit_wave_ <- function(hh_model, wides, reps, houses, seeds, fit_args,
                            outer, bar) {
  local_quiet_warnings_(TRUE)
  fits <- Map(sweep_prepare_one_, wides, seeds,
              MoreArgs = list(hh_model = hh_model, fit_args = fit_args))
  ran <- .engine$run_inference_list(lapply(fits, `[[`, "ptr"), outer, bar)
  if (!is.null(ran$error)) sweep_guard_(
    with_cpp_conditions(stop(ran$error)), reps[[ran$failed]], houses
  )
  Map(\(fit, i) wrap_inference_(fit$model, ran, i), fits, seq_along(fits))
}

#' Build the fit of one survey: its model and engine.
#' @noRd
sweep_prepare_one_ <- function(simulated, seed, hh_model, fit_args) {
  kept <- setdiff(names(simulated), get_latent(hh_model)[["name"]])
  model <- add_data_wide(hh_model, simulated[kept])
  ptr <- do.call(prepare_inference_, c(list(model), fit_args, seed = seed))
  list(model = model, ptr = ptr)
}

#' Fit one replicate, reduced against `truth`.
#' @noRd
sweep_reduce_one_ <- function(draws, truth, rep, seed, houses, reducer) {
  estimates <- cbind(houses = houses, rep = rep, reducer(draws, truth))
  replicate <- data.frame(houses = houses, rep = rep, seed = seed,
                          inference_secs = draws$inference_secs,
                          min_ess_bulk = min(estimates$ess_bulk),
                          flagged = any(estimates$flagged))
  list(estimates = estimates, warnings = draws$warnings, replicate = replicate)
}
