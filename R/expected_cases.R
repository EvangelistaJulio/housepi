#      __
#     / /_  ____  __  __ ____ __  ____  __
#    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#   / / / / /_/ / /_/ (__  )  __/ /_/ / /
#  /_/ /_/\____/\____/____/\___/ ____/_/
#                             /_/
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#
# expected_cases.R
#

#' @title
#' Measure SAR by household size
#' @description
#' Function to retrieve secondary attack rates from data.frame
#' with cases by size.
#'
#' @param dataset `data.frame`\cr
#' Dataset containing individual 'status' and 'house_id'.
#'
#' @return `data.frame`
#'
#' @export
measure_sar_by_house_size <- function(dataset) {
  # check that dataset is consistent
  if (length(dataset) == 0) {
    stop_h("No valid dataset was provided.")
  } else {
    dataset <- as.data.frame(dataset)
  }

  # cast to 0\1 values
  dataset$status <- as.integer(dataset$status > 0)

  # in tidyverse it would be
  # dataset |>
  #   group_by(house_id) |>
  #   summarise(size = n(), cases = sum(status)) |>
  #   filter(cases > 0) |>
  #   mutate(sar = (cases - 1) / (size - 1)) |>
  #   group_by(size) |>
  #   summarise(sar = mean(sar))

  agg1 <- stats::aggregate(
    status ~ house_id, data = dataset,
    FUN = function(x) c(size = length(x), cases = sum(x))
  )
  agg2 <- data.frame(
    house_id = agg1$house_id,
    size = agg1$status[, "size"],
    cases = agg1$status[, "cases"]
  )
  agg2 <- agg2[agg2$cases > 0, ]
  agg2$sar <- (agg2$cases - 1) / (agg2$size - 1)
  return(stats::aggregate(sar ~ size, data = agg2, FUN = mean))
}

#' @title
#' Measure Household Cases by Size
#' @description
#' Function to retrieve cases by size from a given data.frame.
#'
#' @param dataset `data.frame`\cr
#' Dataset containing individual 'status' and 'house_id'.
#'
#' @return `data.frame`
#'
#' @export
measure_cases_by_house_size <- function(dataset) {
  # check that dataset is consistent
  if (length(dataset) == 0) {
    stop_h("No valid dataset was provided.")
  } else {
    dataset <- as.data.frame(dataset)
  }

  # cast to 0\1 values
  dataset$status <- as.integer(dataset$status > 0)

  # in tidyverse it would be
  # dataset |>
  #   group_by(house_id) |>
  #   summarise(size = n(), cases = sum(status)) |>
  #   count(size, cases) |>
  #   complete(size = min(size):max(size), cases = 0:max(size)) |>
  #   mutate(n = if_else(is.na(n) & cases <= size, 0L, n)) |>
  #   pivot_wider(names_from = cases, values_from = n) |>
  #   as.data.frame()

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
#' Frame Apply
#' @description
#' Apply function `func` to corresponding elements
#' across a list of data.frame that share the same dimensions.
#'
#' @param func `closure()`\cr
#' Function applied to corresponding elements across data frames in the list.
#' @param frame_list `list()` | `data.frame`\cr
#' List of data frames, all having the same dimensions.
#'
#' @return `data.frame()`\cr
#' Data frame of the same dimensions as the input frames.
#' Columns and rownames of the first data.frame in the list.
#'
#' @noRd
frame_apply <- function(func, frame_list) {
  elem_n1 <- frame_list[[1]]

  # perform element-wise operations
  if (is.null(dim(elem_n1))) {
    # on scalar/vectors
    output <- sapply(seq_along(elem_n1), function(x) {
      func(sapply(frame_list, "[[", x))
    })
    names(output) <- names(elem_n1)
    if (is.list(elem_n1)) output <- as.list(output)

  } else {
    # on matrix/data.frames
    output <- matrix(
      sapply(seq_along(elem_n1), function(x) {
        apply(sapply(frame_list, "[[", x), 1, func)
      }),
      nrow = nrow(elem_n1),
      byrow = FALSE
    )
    colnames(output) <- colnames(elem_n1)
    rownames(output) <- rownames(elem_n1)
    if (is.data.frame(elem_n1)) output <- as.data.frame(output)
  }

  return(output)
}

#' @title
#' Frame Apply Mran
#' @description
#' Compute the mean across different dataframes.
#'
#' @param frame_list `list()` | `data.frame`\cr
#'
#' @return `data.frame()`
#' @noRd
frame_apply_mean <- function(frame_list) {
  func <- function(x) mean(x, na.rm = TRUE)
  frame_apply(func, frame_list)
}

#' @title
#' Frame Apply Credible Interval
#' @description
#' Compute the credible interval across different dataframes.
#'
#' @param frame_list `list()` | `data.frame`\cr
#' List of data frames, all having the same dimensions.
#' @param cri_level `numeric(1)`\cr
#' List of data frames, all having the same dimensions.
#' @param upper `logical(1)`\cr
#' If `TRUE`, computes the upper credible interval. Else the lower.
#'
#' @return `data.frame()`
#' @noRd
frame_apply_cri <- function(frame_list, cri_level, upper = TRUE) {
  func <- function(x) {
    probs <- 0.5 + ifelse(isTRUE(upper), 0.5, -0.5) * cri_level
    stats::quantile(x, probs = probs, na.rm = TRUE)
  }
  frame_apply(func, frame_list)
}

############################################################################
############################################################################
########################### Housepi Posterior Expec ########################
############################################################################
############################################################################

#' @title
#' Housepi Posterior Expectation
#' @description
#' Class to extract expecteted measures over simulations sampled from posterior.
#'
#' @keywords internal
HousepiPosteriorExpec <- R6::R6Class(
  classname = "HousepiPosteriorExpec",
  cloneable = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .expected = NULL,
    .observed = NULL,
    .cri_level = 0.95,

    # extract expected quantities from cases-by-size table list
    .extract_posterior_results = function(x) {
      cri <- private$.cri_level
      return(list(
        # compute mean values
        value = frame_apply_mean(x),
        # compute lower credible interval
        cri_lower = frame_apply_cri(x, cri, upper = FALSE),
        # compute upper credible interval
        cri_upper = frame_apply_cri(x, cri, upper = TRUE)
      ))
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field observed `list()`\cr
    #' Observed values of the measure of the posterior sample.
    observed = function(value) {
      if (missing(value)) return(private$.observed)
      stop_h("'observed' is read only.")
    },

    #' @field expected `list()`\cr
    #' Expected values of the measure of the posterior sample.
    expected = function(value) {
      if (missing(value)) return(private$.expected)
      stop_h("'expected' is read only.")
    }
  ),

  ############################################################################
  ################################### PUBLIC #################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [HousepiPosteriorExpec][housepi::HousepiPosteriorExpec].
    #'
    #' @param house_model `HouseholdModel`\cr
    #' Object of class [HouseholdModel][housepi::HouseholdModel]
    #' with complete definitions and traces already stored.
    #' @param FUN `closure()` | `list()` of `closure()`\cr
    #' Function to extract the desired measure from a dataset.
    #' @param progress `logical(1)`\cr
    #' Progress of simulations. Default is inherited from the model.
    initialize = function(house_model, FUN, progress = NA) {
      # function to be applied also in the case of multiple measures
      apply_func <- ifelse(
        is.list(FUN),
        function(x) lapply(FUN, function(f) f(x)),
        FUN
      )

      # initialize with results from observed dataset
      private$.observed <- apply_func(house_model$dataset)

      # deep clone model before proceeding
      # FIXME/UPDATE/NOTE: deep clone is not cloning traces e. g.!
      hmc <- house_model$clone(deep = TRUE)
      # cleanup at end
      withr::defer(rm(hmc))
      # extract traces without loglikelihood
      traces <- house_model$traces_no_loglike
      # stuck multiple chains
      if (is.list(traces)) traces <- do.call(rbind, traces)

      # credible interval from the model
      private$.cri_level <- hmc$cri_level

      # if progress not provided, import from model
      progress <- ifelse(is.na(progress), hmc$setup$progress, progress)
      # UPDATE: can't manually set the frequency here
      if (is.numeric(progress)) progress <- TRUE

      # number of repetitions
      n_rep <- dim(traces)[1]
      # generator seeds
      withr::with_seed(
        hmc$setup$seed,
        rng_seeds <- sample(1:max(.Machine$integer.max, n_rep),
                            n_rep, replace = FALSE)
      )

      # extract constant parameters
      aux <- hmc$param_defs
      tmp <- aux[aux$constant, c("name", "value")]
      constants <- as.list(tmp$value)
      names(constants) <- as.list(tmp$name)

      # define loop function working by row index
      loop_func <- function(r) {
        # associate parameters from posterior sample
        new_param_defs <- c(as.list(traces[r, ]), constants)
        hmc$load_param_list(from = new_param_defs, append = FALSE)
        # run simulation and extract cases by size
        apply_func(hmc$simulate(seed = rng_seeds[r], progress = FALSE))
      }

      # apply loop function
      if (isTRUE(progress)) {
        output <- lapply(
          cli::cli_progress_along(seq_len(n_rep), name = "simulations",
                                  type = "tasks", clear = FALSE),
          loop_func
        )
      } else {
        output <- lapply(seq_len(n_rep), loop_func)
      }

      if (is.list(FUN)) {
        output <- lapply(names(FUN), function(x) lapply(output, `[[`, x))
        names(output) <- names(FUN)
        private$.expected <- lapply(output, private$.extract_posterior_results)

      } else {
        private$.expected <- private$.extract_posterior_results(output)
      }

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info(
        "Output of method `$expec_posterior_measure` of HouseholdModel."
      )
      invisible(self)
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title
#' Combined Observed and Expected cases tables.
#' @description
#' Function to combine cases table.
#'
#' @param obsdf `data.frame`\cr
#' Observed cases table.
#' @param expdf `data.frame`\cr
#' Expected cases table.
#'
#' @return `data.frame`\cr
#'
#' @noRd
combine_obs_exp <- function(obsdf, expdf) {
  # 1. curate observed dataframe
  obsdf$distr <- "Observed"
  rownames(obsdf) <- seq(1, 2 * nrow(obsdf), 2)

  # 2. curate expected dataframe
  expdf$size <- obsdf$size
  expdf$distr <- "Expected"
  rownames(expdf) <- seq(2, 2 * nrow(expdf), 2)

  # 3. combine by row bind
  combined <- rbind(obsdf, expdf)[as.character(1:(2 * nrow(obsdf))), ]
  # reset row names
  rownames(combined) <- NULL
  # reorder columns
  y <- colnames(combined)
  aux <- c("size", "distr")
  combined <- combined[, c(aux, y[!(y %in% aux)])]

  # format columns to extract labels with labelizor
  colnames(combined)[3:ncol(combined)] <- sapply(
    colnames(combined)[3:ncol(combined)], function(x) {
      return(paste("cases", x, sep = "_"))
    }
  )
  return(combined)
}

############################################################################
############################################################################
############################ Expected Cases By Size ########################
############################################################################
############################################################################

#' @title
#' Housepi Expected Cases By SIze
#' @description
#' Output object of expected cases by size, method of household model.
#'
#' @keywords internal
ExpectedCasesBySize <- R6::R6Class(
  classname = "ExpectedCasesBySize",
  inherit = HousepiPosteriorExpec,
  cloneable = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    .max_house_size = NA,
    .counts = NULL,

    # @description
    # Number of household with non-zero cases per size.
    #
    # @return `data.frame`
    .count_non_zero_cases = function(dataset) {
      # in tidyverse it would be
      # dataset |>
      #   group_by(house_id) |>
      #   summarise(size = n(), cases = sum(status)) |>
      #   filter(cases > 0) |>
      #   group_by(size) |>
      #   summarise(count = n())

      dataset$status <- as.integer(dataset$status > 0)

      agg1 <- stats::aggregate(
        status ~ house_id, data = dataset,
        FUN = function(x) c(size = length(x), cases = sum(x))
      )
      agg2 <- data.frame(
        house_id = agg1$house_id,
        size = agg1$status[, "size"],
        cases = agg1$status[, "cases"]
      )
      agg2 <- agg2[agg2$cases > 0, ]
      return(stats::aggregate(cases ~ size, data = agg2, FUN = length))
    },

    # @description
    # Return a list with all filtered case tables.
    #
    # @return `list()`
    .filtered_tables = function() {
      # load observed table
      obsdf <- private$.observed$cases_by_size

      # remove rows without counts
      mask_row <- rowSums(obsdf[colnames(obsdf) != "size"], na.rm = TRUE) != 0

      # remove column with zero cases if none was observed
      mask_col <- (colSums(obsdf == 0, na.rm = TRUE) != nrow(obsdf)) | (colnames(obsdf) != "0")

      # implement max household size mask
      mhs <- private$.max_house_size
      if (!is.na(mhs)) {
        col_max <- which(colnames(obsdf) == as.character(mhs))
        row_max <- which(rownames(obsdf) == as.character(mhs))

        mask_col <- mask_col & (seq_len(ncol(obsdf)) <= col_max)
        mask_row <- mask_row & (seq_along(mask_row) <= row_max & mask_row)
      }

      expdf <- private$.expected$cases_by_size
      return(list(
        observed = obsdf[mask_row, mask_col],
        expected = expdf$value[mask_row, mask_col],
        cri_lower = expdf$cri_lower[mask_row, mask_col],
        cri_upper = expdf$cri_upper[mask_row, mask_col]
      ))
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field sar `data.frame`\cr
    #' Observed and expected secondary attack rates by household size.
    #' If a maximum house size is set, larger sizes are not returned.
    sar = function(value) {
      if (missing(value)) {
        # convert internal list to data.frame
        output <- data.frame(
          size = private$.observed$sar$size,
          observed = private$.observed$sar$sar,
          expected = private$.expected$sar$value$sar,
          cri_lower = private$.expected$sar$cri_lower$sar,
          cri_upper = private$.expected$sar$cri_upper$sar
        )

        # associate counts
        output$counts <- private$.counts$cases

        return(output)
      }
      stop_h("'sar' is read only.")
    },

    #' @field formatted_cases_table `data.frame`\cr
    #' Formatted observed and expected cases by household size.
    #' If a maximum house size is set, larger sizes are not returned.
    formatted_cases_table = function(value) {
      if (missing(value)) {
        # filter according to max_house_size
        cases <- private$.filtered_tables()

        # empty frame as a placeholder for Credible Interval
        na_df <- cases$observed
        na_df[] <- lapply(na_df, function(col) if (is.numeric(col)) NA else col)
        na_df$size <- cases$observed$size

        return(list(expected = combine_obs_exp(cases$observed, cases$expected),
                    cri_lower = combine_obs_exp(na_df, cases$cri_lower),
                    cri_upper = combine_obs_exp(na_df, cases$cri_upper)))
      }
      stop_h("'formatted_case_table' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @description
    #' New instance of class
    #' [ExpectedCasesBySize][housepi::ExpectedCasesBySize].
    #'
    #' @param house_model `HouseholdModel`\cr
    #' Object of class
    #' [HouseholdModel][housepi::HouseholdModel]
    #' containing all model definitions after having run `$infer`.
    #' @param progress `logical(1)`\cr
    #' Progress of simulations. Default is inherited from the model.
    initialize = function(house_model, progress = NA) {
      # ininitialize parent class HousepiPosteriorExpec
      super$initialize(
        house_model,
        FUN = c(cases_by_size = measure_cases_by_house_size,
                sar = measure_sar_by_house_size),
        progress = progress
      )
      private$.counts <- private$.count_non_zero_cases(house_model$dataset)

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info(
        "Output of method `$expec_cases_by_house_size` of HouseholdModel."
      )
      invisible(self)
    },

    #' @description
    #' Set maximum household size to consider in result outputs.
    #'
    #' @param n `integer(1)`\cr
    #' Maximum household size. If `NA`, restore observed maximum.
    set_max_house_size = function(n) {
      private$.max_house_size <- n
      invisible(self)
    },

    #' @description
    #' Restore maximum household size to observed maximum.
    drop_max_house_size = function() {
      self$set_max_house_size(NA)
      invisible(self)
    },

    #' @description
    #' Perform Chi-squared test and return list with results.
    #'
    #' @param lower_tail `logical(1)`\cr
    #' if TRUE (default), probabilities are P[X<=x], otherwise P[X>x].
    #'
    #' @return
    #' `list()`)
    chisq_test = function(lower_tail = FALSE) {
      x <- private$.filtered_tables()
      res <- (x$observed - x$expected)^2 / x$expected
      chisq_value <- sum(res, na.rm = TRUE)
      ddof <- prod(dim(res)) - dim(res)[1] - sum(is.na(res))
      return(list(
        chisq_value = chisq_value, ddof = ddof,
        p_value = stats::pchisq(chisq_value, df = ddof, lower.tail = lower_tail)
      ))
    }
  )
)
