#
# analysis.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title
#' Coerce traces to coda.
#' @description
#' Convert HouseholdModel traces to coda.
#'
#' @param trs
#' Traces as exported from HouseholdModel.
#'
#' @return
#' (`coda::mcmc` | `coda::mcmc.list`)\cr
#'
#' @noRd
coerce_traces_to_coda <- function(trs) {
  if (requireNamespace("coda", quietly = TRUE)) {
    if ((methods::is(trs, "mcmc") || methods::is(trs, "mcmc.list"))) {
      return(trs)

    } else if (is.matrix(trs) || (is.data.frame(trs))) {
      sweeps <- as.integer(rownames(df))
      return(coda::mcmc(data = trs,
                        start = sweeps[1],
                        end = sweeps[length(sweeps)],
                        thin = sweeps[2] - sweeps[1]))

    } else if (is.list(trs)) {
      sweeps <- as.integer(rownames(trs[[1]]))
      df_list <- lapply(trs, function(x) {
        coda::mcmc(
          data = x,
          start = sweeps[1],
          end = sweeps[length(sweeps)],
          thin = sweeps[2] - sweeps[1]
        )
      })
      return(coda::as.mcmc.list(df_list))
    }
  } else {
    stop.h(
      "Impossible to coerce to coda object.",
      "Package 'coda' not installed."
    )
  }
}

#' @title
#' Extract results from coda
#' @description
#' Collects relevant statistics from the the input mcmc traces.
#'
#' @param trs
#' Traces as exported from HouseholdModel.
#' @param credible_interval (`double(1)`)\cr
#' Probability for credible interval (default is 0.95).
#'
#' @return
#' (`data.frame`)
#'
#' @noRd
retrieve_traces_results <- function(trs,
                                    credible_interval = 0.95) {
  # change in case of multiple chains
  if (is.list(trs)) trs <- do.call(rbind, trs)

  cri_lower <- 0.5 - 0.5 * credible_interval
  cri_upper <- 0.5 + 0.5 * credible_interval
  out <- data.frame(
    mean = apply(trs, 2, mean, na.rm = TRUE),
    sd = apply(trs, 2, sd, na.rm = TRUE),
    low = apply(trs, 2, quantile, probs = cri_lower, na.rm = TRUE),
    up = apply(trs, 2, quantile, probs = cri_upper, na.rm = TRUE),
    check.names = FALSE
  )
  colnames(out)[3] <- sprintf("cri%.1f", cri_lower * 100)
  colnames(out)[4] <- sprintf("cri%.1f", cri_upper * 100)

  return(out)
}

#' @title
#' Apply on Homogeneous Frames
#' @description
#' Apply function `func` to corresponding elements
#' across a list of data frames that share the same dimensions.
#'
#' @param func (`closure()`)\cr
#' Function applied to corresponding elements across data frames in the list.
#' @param frame_list (`list()` | `data.frame`)\cr
#' List of data frames, all having the same dimensions.
#'
#' @return (`data.frame()`)\cr
#' Data frame of the same dimensions as the input frames.
#'
#' @noRd
apply_on_homogeneous_frames <- function(func, frame_list) {
  # perform element-wise operations
  output <- matrix(
    sapply(seq_along(frame_list[[1]]), function(x) {
      apply(sapply(frame_list, "[[", x), 1, func)
    }),
    nrow = nrow(frame_list[[1]]),
    byrow = FALSE
  )
  colnames(output) <- colnames(frame_list[[1]])
  return(as.data.frame(output))
}

#' @title
#' Retrieve Household Cases by Size
#' @description
#' Function to retrieve cases by size from a given data.frame.
#'
#' @param dataset (`data.frame`)\cr
#' Dataset containing the the individual 'id'
#' their 'status' and their 'house_id'.
#'
#' @return (`data.frame`)\cr
#'
#' @export
retrieve_cases_by_size <- function(dataset) {
  # check that dataset is consistent
  if (length(dataset) == 0) {
    stop_h("No valid dataset was provided.")
  } else {
    dataset <- as.data.frame(dataset)
  }

  # cast to 0\1 values
  dataset$status <- as.integer(dataset$status > 0)

  # aggregate data by house_id
  house_df <- stats::aggregate(status ~ house_id, data = dataset, FUN = sum)
  tmp_size <- table(dataset$house_id)
  house_df$size <- tmp_size[match(house_df$house_id, names(tmp_size))]

  # aggregate by size and cases
  cases_by_size <- as.data.frame(table(house_df$size, house_df$status))
  colnames(cases_by_size) <- c("size", "cases", "counts")

  # determine max size
  max_size <- max(as.integer(cases_by_size$size), na.rm = TRUE)

  # create a complete grid of size and cases
  grid <- expand.grid(size = 1:max_size, cases = 0:max_size)
  grid <- grid[grid$cases <= grid$size, ]

  # merge with cases_by_size and replace NAs with 0
  merged_data <- merge(grid, cases_by_size,
                       by = c("size", "cases"), all.x = TRUE)
  merged_data$counts[is.na(merged_data$counts)] <- 0

  # reshape to wide format
  result <- stats::reshape(merged_data,
                           idvar = "size",
                           timevar = "cases",
                           direction = "wide")

  # rename columns to match tidyverse output
  colnames(result) <- sub("counts\\.", "", colnames(result))

  # ensure rows are sorted by size
  output <- result[order(result$size), ]

  # drop row names
  rownames(output) <- NULL

  return(output)
}

############################################################################
############################################################################
############################################################################

#' @title
#' Housepi Output Expected Cases
#' @description
#' Output object from expected cases by size.
#'
#' @keywords internal
HousepiOutExpCases <- R6::R6Class(
  classname = "HousepiOutExpCases",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .obs_hh_case_table = NULL,
    .exp_hh_case_table = NULL,
    .q_lower = NULL,
    .q_upper = NULL,
    .traces = NULL, # temporary

    # @dscription
    # Assign a filtering on the case table.
    #
    # @param max_house_size (`integer(1)`)\cr
    # Maximum household size to consider.
    .filter_table = function(max_house_size = NA) {
      # load observed table
      obsdf <- private$.obs_hh_case_table

      # FIXME: what if size != 0 has all zeros!?
      mask_row <- rowSums(obsdf[colnames(obsdf) != "size"], na.rm = TRUE) != 0
      # implement max household size mask
      if (is.na(max_house_size)) {
        mask_col <- TRUE
      } else {
        col_max <- which(colnames(obsdf) == as.character(max_house_size))
        row_max <- which(rownames(obsdf) == as.character(max_house_size))

        mask_col <- seq_len(ncol(obsdf)) <= col_max
        mask_row <- seq_along(mask_row) <= row_max & mask_row
      }

      return(list(
        observed = obsdf[mask_row, mask_col],
        expected = private$.exp_hh_case_table[mask_row, mask_col],
        cri_lower = private$.q_lower[mask_row, mask_col],
        cri_upper = private$.q_upper[mask_row, mask_col]
      ))
    }
  ),

  ############################################################################
  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [HousepiOutExpCases][housepi::HousepiOutExpCases].
    #'
    #' @param house_model (`HouseholdModel`)\cr
    #' Object of class \[HouseholdModel][housepi::HouseholdModel]
    #' containing all model definitions after having run `$infer`.
    initialize = function(house_model,
                          credible_interval = 0.95,
                          progress = NA) {
      # extract observed household case table
      private$.obs_hh_case_table <- retrieve_cases_by_size(house_model$dataset)

      # extract traces
      traces <- house_model$traces

      if (is.matrix(traces) || (is.data.frame(traces))) {
        n_chains <- 1
        sample_size <- dim(traces)[1]
      } else {
        n_chains <- length(traces)
        sample_size <- dim(traces[[1]])[1]
      }

      # FLAG
      # FIXME: deep clone not working!
      # clone the model before proceeding further
      # hmc <- house_model$clone(deep = TRUE)
      # cleanup at end
      # withr::defer(rm(hmc))
      hmc <- house_model
      private$.traces <- traces # FIXME: temporary backup
      # FLAG-END

      # number of repetitions
      n_rep <- sample_size * n_chains

      # progress?
      progress <- isTRUE(progress) || 

      # generator seeds
      withr::with_seed(
        house_model$setup$seed,
        rng_seeds <- sample(1:max(1e8, n_rep), n_rep, replace = FALSE)
      )

      # extract traces columns
      traces_cols <- colnames(if (n_chains > 1) traces[[1]] else traces)
      mask_traces <- !(traces_cols %in% c("loglikelihood", "CHAIN__"))

      # extract constant parameters
      aux <- house_model$param_defs
      tmp <- aux[aux$constant, c("name", "value")]
      constants <- as.list(tmp$value)
      names(constants) <- as.list(tmp$name)

      # store simulated datasets
      if (isTRUE(progress)) {
        cli::cli_progress_bar("simulations", total = n_rep,
                              type = "tasks", clear = FALSE)
      }
      simulations_list <- list()

      for (c in seq_len(n_chains)) {
        for (s in seq_len(sample_size)) {
          r <- s + sample_size * (c - 1)
          traces_row <- if (n_chains > 1) traces[[c]][s, ] else traces[s, ]

          # associate parameters from posterior sample
          new_param_defs <- c(as.list(traces_row[mask_traces]), constants)
          hmc$load_param_list(from = new_param_defs, append = FALSE)
          # run simulation
          hmc$simulate(seed = rng_seeds[r], progress = FALSE)
          # append tables with cases by size
          simulations_list[[r]] <- retrieve_cases_by_size(hmc$dataset)

          if (isTRUE(progress)) cli::cli_progress_update()
        }
      }
      if (isTRUE(progress)) cli::cli_progress_done()

      # compute mean values
      private$.exp_hh_case_table <- apply_on_homogeneous_frames(
        function(x) {
          mean(x, na.rm = TRUE)
        }, simulations_list
      )

      # internal consistency check
      if (all(dim(private$.obs_hh_case_table) != dim(private$.exp_hh_case_table))) {
        stop_h("Internal inconsistency in observed and expected tables.")
      }

      # compute lower credible interval
      probs_lower <- 0.5 - 0.5 * credible_interval
      private$.q_lower = apply_on_homogeneous_frames(
        function(x) {
          stats::quantile(x, probs = probs_lower, na.rm = TRUE)
        }, simulations_list
      )

      # compute upper credible interval
      probs_upper <- 0.5 + 0.5 * credible_interval
      private$.q_upper = apply_on_homogeneous_frames(
        function(x) {
          stats::quantile(x, probs = 0.975, na.rm = TRUE)
        }, simulations_list
      )

      invisible(self)
    },

    #' @description
    #' Override the default printing behaviour.
    print = function() {
      cli::cli_alert_info(
        "Output of method `$expec_cases_by_size` of HouseholdModel."
      )
      invisible(self)
    },

    #' @dscription
    #' Perform Chi-squared test.
    #'
    #' @param max_house_size (`integer(1)`)\cr
    #' Maximum household size to consider.
    chisq_test = function(max_house_size = NA) {
      x <- private$.filter_table(max_house_size)
      res <- (x$observed - x$expected)^2 / x$expected
      return(list(
        chisq_value = sum(res, na.rm = TRUE),
        ddof = prod(dim(res)) - dim(res)[1] - sum(is.na(res)),
        p_value = stats::pchisq(chisq_value, df = ddof, lower.tail = FALSE)
      ))
    },

    #' @description
    #' Return a `flextable` combining the results of the check on household model
    #'
    #' @param max_house_size (`integer(1)`)\cr
    #' Maximum household size to consider.
    #'
    #' @return
    #' (`flextable` | `data.frame`)
    show_model_check_table = function(max_house_size = NA,
                                      show_chisq = FALSE,
                                      cri_format = "%.0f [%.0f, %.0f]") {
      # format counts to character
      format_obs <- function(x) {
        as.character(as.integer(x, na.rm = TRUE), na.rm = TRUE)
      }
      # format `mean` with its `lower` and `upper` values
      format_exp <- function(mean, lower, upper) {
        return(ifelse(is.na(mean) | is.na(lower) | is.na(upper),
                      NA, sprintf(cri_format, mean, lower, upper)))
      }

      # filter according to max_house_size
      cases <- private$.filter_table(max_house_size)

      # curate observed dataframe
      obsdf <- as.data.frame(apply(cases$observed, 2, format_obs))
      obsdf$distr <- "Observed"
      rownames(obsdf) <- seq(1, 2 * nrow(obsdf), 2)

      # curate expected dataframe
      expdf <- as.data.frame(
        mapply(format_exp, cases$expected, cases$cri_lower, cases$cri_upper)
      )
      expdf$size <- obsdf$size
      expdf$distr <- "Expected"
      rownames(expdf) <- seq(2, 2 * nrow(expdf), 2)

      # combine data frames
      combined <- rbind(obsdf, expdf)[as.character(1:(2 * nrow(obsdf))), ]

      # reset row names
      rownames(combined) <- NULL

      # reorder columns
      y <- colnames(combined)
      aux <- c("size", "distr")
      combined <- combined[, c(aux, y[!(y %in% aux)])]

      # format to extract labels with labelizor
      colnames(combined)[3:ncol(combined)] <- sapply(
        colnames(combined)[3:ncol(combined)], function(x) {
          return(paste("cases", x, sep="_"))
        }
      )

      if (requireNamespace("flextable", quietly = TRUE)) {
        ft <- flextable::flextable(combined)
        ft <- flextable::separate_header(ft)
        ft <- flextable::theme_vanilla(ft)
        ft <- flextable::align(ft, align = "center", part = "all")
        ft <- flextable::bold(ft, j = ~ size, bold = TRUE)
        ft <- flextable::merge_v(ft, j = ~ size)
        ft <- flextable::labelizor(
          ft, part = "all",
          labels = c(size = "Size",
                     distr = "Distribution",
                     cases = "Number of cases")
        )
        ft <- flextable::autofit(ft, add_w = 0, add_h=0)

        if (isTRUE(show_chisq)) {
          x <- self$chisq_test(max_house_size)
          footer <- sprintf("X-squared = %.1f, df = %.0f, p-value = %f",
                            x$chisq_value, x$ddof, x$p_value)
          ft <- flextable::add_footer_lines(ft, footer)
        }

        return(ft)

      } else {
        return(combined)
      }
    }
  )
)

#' @title
#' Print S3 method for housepi_model_check objects
#' @description
#' Provide a custom print method for objects of class `housepi_model_check`.
#'
#' @param x
#' An object of class `housepi_model_check`.
#' @param ...
#' Other arguments passed to print.
#'
#' @export
print.housepi_model_check <- function(x, ...) {
  cat("Household model check with Chi-squared test:\n")
  cat(sprintf("X-squared = %.1f, df = %s, p-value = %s",
              x$X_squared, x$df, x$p_value))
  other <- list(...)
  if (length(other) > 0) print(other)
}
