## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

# Replicate sweep of the checks that simulate household studies
# and fit them back: one simulation call per household count,
# one fit per replicate, each fit reduced to one row per parameter.

# Roxygen text shared by the checks, read at the use site as `r <name>`.

rd_sweep_nrep_ <- "`integer(1)` Number of replicates at each household count.
Each replicate is a household study simulated and fitted back."
rd_sweep_sampler_ <- "`function`
Sampler of the household structure the outbreaks are simulated on,
called as `structure_sampler(data, houses)`.
If `NULL`, the attached structure is reused verbatim and `houses` ignored."
rd_sweep_structure_ <- "Only transmission and observations are simulated.
Households, members and their covariates come from the attached data,
drawn once per household count by `structure_sampler`,
or reused whole without it, and shared by the replicates at that count."
rd_sweep_seed_ <- "Replicates are independent draws.
Call `set.seed()` before the check to replay a whole sweep,
at any `threads` up to `nrep`."
rd_sweep_reporting_ <- "Each simulated study goes through
the `dataset` reporting process of the model before it is fitted back,
see [add_reporting_process()].
Inference does not model it,
so a process losing information makes the check measure the misspecification."
rd_sweep_threads_ <- "`integer(1)` Number of threads.
Up to `threads` replicates are fitted at once.
Beyond `nrep`, the households of each fit share the threads left over."

#' Columns the outbreaks are simulated on: the attached data without
#' the observations and latent variables, which simulation regenerates,
#' and which a definition reading one of them must not see observed.
#' Status is kept as simulation receives it, drawn anew by transmission.
#' @noRd
sweep_base_structure_ <- function(hh_model) {
  data <- get_dataset(hh_model)
  if (!NROW(data)) stop_h(
    "No data is attached to the model.",
    paste("The check simulates its outbreaks on the household structure",
          "of the attached data: add one with {.fn add_dataset}.")
  )
  written <- c(get_observations(hh_model)[["name"]],
               get_latent(hh_model)[["name"]])
  data[setdiff(names(data), written)]
}

#' Simulation table: `values` repeated, one row per replicate.
#' * parameter without a column -> drawn from its prior.
#' @noRd
sweep_param_table_ <- function(values, nrep) {
  table <- data.frame(row.names = seq_len(nrep))
  table[names(values)] <- as.list(values)
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
sweep_diagnostics_ <- function(fit, variables) {
  diag <- fit$.__enclos_env__$private$diagnostics_
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

#' Household counts of a sweep, each with its number of replicates.
#' @noRd
sweep_cells_ <- function(replicates) {
  counts <- table(replicates[["houses"]])
  list(houses = as.integer(names(counts)), nrep = as.vector(counts))
}

#' Fingerprint of what the studies are simulated and fitted with:
#' * engine -> its module name;
#' * household structure the outbreaks are simulated on -> hashed;
#' * data every fit is conditioned on -> hashed, `NULL` for none;
#' * definitions as the engine imports them -> hashed per family;
#' * structure sampler and `dataset` reporting process -> kept whole;
#' * sampler settings and `time_step` -> one entry each.
#' @param settings `list` Sampler settings and `time_step` of the fits.
#' @noRd
sweep_fingerprint_ <- function(hh_model, structure, sampler, conditioned,
                               settings) {
  pvt <- hh_model$.__enclos_env__$private
  c(list(model = pvt$compiled_$module,
         household_data = rlang::hash(structure),
         conditioned_data =
           if (!is.null(conditioned)) rlang::hash(conditioned)),
    hash_payload_(pvt$payload_()),
    list(structure_sampler = sampler,
         reporting_process = pvt$reporting_process_$dataset),
    lapply(settings, as.numeric))
}

#' Data a study is fitted on: the households of `observed`,
#' then those of `simulated`, numbered past them.
#' * `observed` NULL -> `simulated` unchanged;
#' * numeric `hh_id` or `id` on both sides -> offset by the observed maximum,
#'   else character labels distinct from the observed ones;
#' * column on one side only -> `NA` on the other, of the same class.
#' @noRd
sweep_condition_ <- function(observed, simulated) {
  if (is.null(observed)) return(simulated)
  relabel <- function(sim, obs) {
    key <- match(sim, unique(sim))
    if (is.numeric(sim) && is.numeric(obs)) return(max(obs) + key)
    obs <- unique(as.character(obs))
    fresh <- make.unique(c(obs, paste0("sim_", seq_len(max(key)))))
    fresh[-seq_along(obs)][key]
  }
  for (key in intersect(c("hh_id", "id"), names(simulated))) {
    simulated[[key]] <- relabel(simulated[[key]], observed[[key]])
    if (!is.numeric(simulated[[key]]))
      observed[[key]] <- as.character(observed[[key]])
  }
  cols <- union(names(observed), names(simulated))
  fill <- function(d, other) {
    for (col in setdiff(cols, names(d)))
      d[[col]] <- rep(other[[col]][NA_integer_], NROW(d))
    d[cols]
  }
  out <- rbind(fill(observed, simulated), fill(simulated, observed))
  rownames(out) <- NULL
  out
}

#' Sweep of a recovery or calibration check,
#' held in the private field `sweep_` of the check.
#' @noRd
HousepiSweep <- R6::R6Class(
  classname = "HousepiSweep",

  ####  Private  ####

  private = list(
    # what the studies are simulated and fitted with, see sweep_fingerprint_()
    fingerprint_ = NULL,
    estimates_ = NULL,
    replicates_ = NULL,
    warnings_ = NULL,
    secs_ = NA_real_
  ),

  ####  Active  ####

  active = list(
    # `data.frame` One row per household count, replicate and parameter.
    estimates = function(value) {
      assert_read_only_(value, "estimates")
      private$estimates_
    },

    # `data.frame` One row per household count and replicate.
    replicates = function(value) {
      assert_read_only_(value, "replicates")
      private$replicates_
    },

    # `data.frame` Distinct warnings, with their counts.
    warnings = function(value) {
      assert_read_only_(value, "warnings")
      private$warnings_
    },

    # `numeric(1)` Wall-clock of the sweep.
    secs = function(value) {
      assert_read_only_(value, "secs")
      private$secs_
    },

    # `data.frame` Sweep cost per household count.
    overview = function(value) {
      assert_read_only_(value, "overview")
      sweep_overview_(private$replicates_)
    },

    # `integer()` Household counts swept.
    houses = function(value) {
      assert_read_only_(value, "houses")
      sweep_cells_(private$replicates_)$houses
    },

    # `integer()` Number of replicates per household count, order of `houses`.
    nrep = function(value) {
      assert_read_only_(value, "nrep")
      sweep_cells_(private$replicates_)$nrep
    }
  ),

  ####  Public  ####

  public = list(
    # Simulate `nrep` studies per household count, fit each one back
    # and reduce it, then report the warnings under the name of the check.
    # The outbreaks are simulated on the structure of the attached data.
    # The fitted model is built from `hh_model`,
    # never from the simulating one, holding the truth of the last replicate.
    # Fits run in waves of replicates sampled at once,
    # each wave reduced before the next one runs.
    # Structure, household counts and sampler settings are asserted
    # before any study is simulated,
    # and the first replicate that fails stops the sweep.
    # @param houses `integer()` Household counts, `NA` for the attached one.
    # @param values `data.frame` Simulation values, one row per replicate.
    # A parameter without a column is drawn from its prior.
    # @param fit_args `list` Sampler settings of each fit.
    # @param threads `integer(1)` Threads shared by the fits of a wave.
    # @param dots `list` Advanced settings of the check,
    # split between the simulation and the fits.
    # @param reducer `function(fit, truth)` Reduce the posterior draws `fit`
    # of one replicate against `truth`, the named values it was simulated at.
    # It returns one row per parameter, with columns `ess_bulk` and `flagged`.
    # @param label `character(1)` Label of the progress bar.
    # @param what `character(1)` Name of the check, heading its warnings.
    # @param conditioned `data.frame` Data every study is fitted alongside,
    # `NULL` for none.
    # @param one `logical(1)` Demand a single household count.
    initialize = function(hh_model, sampler, houses, nrep, values, fit_args,
                          threads, dots, reducer, progress, label, what,
                          conditioned = NULL, one = FALSE) {
      structure <- sweep_base_structure_(hh_model)
      houses <- sweep_houses_(structure, houses, is.null(sampler), one)
      do.call(resolve_settings_, fit_args)
      # fits read the dataset only
      hh_model <- add_reporting_process(hh_model, long_dataset = NULL)
      time_step <- dots[names(dots) == "time_step"]
      # a copy carries no settings: it reads the engine defaults
      advanced <- utils::modifyList(hh_model$settings["time_step"], time_step)
      private$fingerprint_ <- sweep_fingerprint_(
        hh_model, structure, sampler, conditioned, c(fit_args, advanced)
      )
      workers <- sweep_split_(nrep, threads)
      fit_args <- c(fit_args, list(threads = workers$inner, id_draws = NULL))
      waves <- split(seq_len(nrep), (seq_len(nrep) - 1L) %/% workers$outer)

      # the replicates at `n` households, reduced
      cell <- function(n, bar) {
        struct <- if (is.null(sampler)) structure else sampler(structure, n)
        sim_model <- add_dataset(hh_model, struct)
        studies <- do.call(sim_model$.__enclos_env__$private$simulate_,
                           c(list(draws = values), time_step))
        truths <- lapply(seq_len(nrep),
                         \(rep) unlist(studies$parameters[rep, , drop = FALSE]))
        # seeds drawn ahead of the fits: the reducers draw from the same stream
        # whatever the size of a wave
        seeds <- draw_seed_(nrep)
        reduced <- vector("list", nrep)
        for (reps in waves) {
          fits <- sweep_fit_wave_(
            hh_model, studies$dataset[reps], reps = reps, houses = n,
            seeds = seeds[reps], fit_args = fit_args, outer = workers$outer,
            bar = bar, conditioned = conditioned
          )
          reduced[reps] <- Map(sweep_reduce_one_, fits, truths[reps], reps,
                               seeds[reps],
                               MoreArgs = list(houses = n, reducer = reducer))
        }
        reduced
      }

      start <- proc.time()[["elapsed"]]
      ran <- with_progress_(progress, label, length(houses) * nrep, \(bar) {
        unlist(lapply(houses, cell, bar = bar), recursive = FALSE)
      })
      private$secs_ <- proc.time()[["elapsed"]] - start
      field <- \(k) lapply(ran, `[[`, k)
      keys <- c("houses", "rep")
      private$estimates_ <- sweep_stack_(field("estimates"), keys)
      private$replicates_ <- sweep_stack_(field("replicate"), keys)
      private$warnings_ <- sweep_warning_counts_(field("warnings"))
      sweep_report_warnings_(private$warnings_, what)
    },

    # Print `header`, the household counts and replicates, `rows`,
    # then the diagnostics, run times and warnings.
    # @param rows `character()` Rows of the check, glued in `.envir`.
    # Shared rows read this sweep as `sweep`.
    print_rows = function(header, rows, .envir = parent.frame()) {
      env <- list2env(list(sweep = self), parent = .envir)
      print_rows_(header, c(
        Households = "{sweep$houses}",
        Replicates = "{nrow(sweep$replicates)}",
        rows,
        `Flagged by diagnostics` = "{sum(sweep$replicates[['flagged']])}",
        `Run time` = "{format_secs_(sweep$secs)}",
        `Inference time per fit` =
          "{format_secs_(mean(sweep$replicates[['inference_secs']]))}"
      ), warnings = private$warnings_, .envir = env)
    }
  )
)

#' Combine checks of one class into a clone of the first one,
#' their sweeps pooled by sweep_pool_(), the fields of `pooled` concatenated.
#' * any other field differs -> error, see combine_r6_().
#' @param checks `list` Checks to combine, in order.
#' @param pooled `character()` Private fields concatenated across the checks.
#' @noRd
sweep_combine_ <- function(checks, pooled = character()) {
  comb <- combine_r6_(checks, c("sweep_", pooled), "checks")
  comb$pvt$sweep_ <- sweep_pool_(lapply(comb$privs, `[[`, "sweep_"))
  for (field in pooled)
    comb$pvt[[field]] <- unlist(lapply(comb$privs, `[[`, field))
  comb$out
}

#' Pool sweeps into a clone of the first one.
#' Replicates pool, numbered in sequence per household count,
#' warning counts and run times add up.
#' * fingerprints differ -> error, see combine_r6_();
#' * seed repeated at one household count -> error.
#' @param sweeps `list` Sweeps to pool, in order.
#' @noRd
sweep_pool_ <- function(sweeps) {
  pooled <- c("estimates_", "replicates_", "warnings_", "secs_")
  comb <- combine_r6_(sweeps, pooled, "checks")
  privs <- comb$privs
  # each check numbers 1..n per count, offset by earlier checks at that count
  shift <- function(d, i) {
    earlier <- unlist(lapply(privs[seq_len(i - 1L)], \(p) p$replicates_$houses))
    offset <- vapply(d[["houses"]], \(h) sum(earlier == h), integer(1L))
    d[["rep"]] <- d[["rep"]] + offset
    d
  }
  shifted <- \(field) Map(\(p, i) shift(p[[field]], i), privs, seq_along(privs))
  replicates <- sweep_stack_(shifted("replicates_"), c("houses", "rep"))
  repeated <- duplicated(replicates[c("houses", "seed")])
  if (any(repeated)) {
    # nolint next: object_usage_linter. read by the cli glue below
    houses <- unique(replicates[["houses"]][repeated])
    stop_h(
      "The checks repeat a replicate at {.val {houses}} households.",
      "Checks run after the same {.fn set.seed} call draw the same replicates."
    )
  }
  comb$pvt$estimates_ <- sweep_stack_(shifted("estimates_"), c("houses", "rep"))
  comb$pvt$replicates_ <- replicates
  comb$pvt$warnings_ <- sweep_warning_counts_(lapply(privs, `[[`, "warnings_"))
  comb$out
}

#' @title Combine checks of one model
#' @description
#' Pool the replicates of several recovery checks,
#' or of several calibration checks,
#' so a check grows without running its replicates again.
#' @param ... Objects of class [HousepiRecoveryCheck],
#' or of class [HousepiCalibrationCheck].
#' @return A check of the same class, holding the replicates of all checks.
#' @details
#' The checks must share the model, its definitions, the attached data,
#' the structure sampler, the reporting process and the sampler settings.
#' Recovery checks must also share the values the studies are simulated at.
#' Calibration checks must share their source, prior or posterior,
#' and posterior checks rank the draws of distinct runs of [run_inference()].
#' A sampler or a reporting process created inside another function
#' combines only with itself, within one R session.
#' Household counts may differ.
#' Replicates continue the numbering of their household count.
#' Checks run after the same `set.seed()` call draw the same replicates,
#' and do not combine.
#' @examples
#' hh_model <- flu_1918_frederick |>
#'   new_housepi_model("final-size") |>
#'   add_parameter(import = "template")
#' draws <- c(alpha = 0.15, beta = 0.2, gamma = 1)
#' recovery <- check_recovery(hh_model, draws, nrep = 2, houses = 30)
#' recovery_2 <- check_recovery(hh_model, draws, nrep = 2, houses = 30)
#' c(recovery, recovery_2)$nrep
#' @name combine_checks
#' @export
c.HousepiRecoveryCheck <- function(...) sweep_combine_(list(...))

#' @rdname combine_checks
#' @export
c.HousepiCalibrationCheck <- function(...) {
  checks <- list(...)
  seeds <- unlist(lapply(checks, \(x) x$.__enclos_env__$private$draws_seed_))
  if (anyDuplicated(seeds)) stop_h(
    "The checks rank draws of the same run.",
    "Run each posterior check on the draws of its own {.fn run_inference} call."
  )
  sweep_combine_(checks, "draws_seed_")
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

#' Run one wave: the studies of `datasets` built in order,
#' each fitted alongside `conditioned`,
#' sampled `outer` at a time, housepi warnings held back.
#' * first study failing -> stops every fit of the wave.
#' @param reps `integer()` Replicate of each study.
#' @param houses `integer(1)` Household count of the studies.
#' @param seeds `integer()` Engine seed of each study.
#' @param bar Progress bar each fit fills by one unit, `NULL` for none.
#' @param conditioned `data.frame` Data every study is fitted alongside,
#' `NULL` for none.
#' @return Posterior draws of each study.
#' @noRd
sweep_fit_wave_ <- function(hh_model, datasets, reps, houses, seeds, fit_args,
                            outer, bar, conditioned = NULL) {
  local_quiet_warnings_(TRUE)
  fits <- Map(sweep_prepare_one_, datasets, seeds,
              MoreArgs = list(hh_model = hh_model, conditioned = conditioned,
                              fit_args = fit_args))
  ran <- .engine$run_inference_list(lapply(fits, `[[`, "ptr"), outer, bar)
  if (!is.null(ran$error)) sweep_guard_(
    with_cpp_conditions(stop(ran$error)), reps[[ran$failed]], houses
  )
  Map(\(fit, i) wrap_inference_(fit$model, ran, i), fits, seq_along(fits))
}

#' Build the fit of one study, alongside `conditioned`: its model and engine.
#' @noRd
sweep_prepare_one_ <- function(simulated, seed, hh_model, conditioned,
                               fit_args) {
  kept <- setdiff(names(simulated), get_latent(hh_model)[["name"]])
  model <- add_dataset(hh_model, sweep_condition_(conditioned, simulated[kept]))
  ptr <- do.call(prepare_inference_, c(list(model), fit_args, seed = seed))
  list(model = model, ptr = ptr)
}

#' Reduce the fit of one replicate against `truth`.
#' @noRd
sweep_reduce_one_ <- function(fit, truth, rep, seed, houses, reducer) {
  estimates <- cbind(houses = houses, rep = rep, reducer(fit, truth))
  replicate <- data.frame(houses = houses, rep = rep, seed = seed,
                          inference_secs = fit$inference_secs,
                          min_ess_bulk = min(estimates$ess_bulk),
                          flagged = any(estimates$flagged))
  list(estimates = estimates, warnings = fit$warnings, replicate = replicate)
}
