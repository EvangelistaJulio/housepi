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

#' @title Measure SAR by household size
#' @description
#' Function to retrieve secondary attack rates from data.frame
#' with cases by household size.
#'
#' @param dataset `r rd_add_dataset_`
#'
#' @return `data.frame`
#' 
#' @details 
#' The secondary attack rate is computed for all household with size \eqn{n > 1} according to
#' \deqn{\frac{\mathbb{E[n_{+}\vert{n}] - 1}{n - 1}}}{(E[n_pos|n] - 1) / (n - 1)}
#'
#' @export
measure_sar_by_house_size <- function(dataset) {
  # check that dataset is consistent
  if (length(dataset) == 0) stop_h("No valid dataset was provided.")
  else dataset <- as.data.frame(dataset)

  dataset |>
    # cast to 0\1 values
    dplyr::mutate(
      status =  dplyr::if_else(is.na(status), 0, as.integer(status > 0))
    ) |>
    # aggregate data by house_id
    dplyr::group_by(house_id) |>
    dplyr::summarise(house_size =  dplyr::n(), cases = sum(status)) |>
    # compute sar according to definition
    dplyr::filter(cases > 0) |>
    dplyr::mutate(sar = (cases - 1) / (house_size - 1)) |>
    dplyr::group_by(house_size) |>
    dplyr::summarise(sar = mean(sar)) |> 
    as.data.frame()
}

#' @title Measure household cases by size
#' @description
#' Function to retrieve cases by household size from a given data.frame.
#'
#' @param dataset `r rd_add_dataset_`
#'
#' @return `data.frame`
#'
#' @export
measure_cases_by_house_size <- function(dataset) {
  # check that dataset is consistent
  if (length(dataset) == 0) stop_h("No valid dataset was provided.")

  housepi::survey |>
    # cast to 0\1 values
    dplyr::mutate(
      status =  dplyr::if_else(is.na(status), 0, as.integer(status > 0))
    ) |>
    # aggregate data by house_id
    dplyr::group_by(house_id) |>
    dplyr::summarise(house_size = dplyr::n(), cases = sum(status)) |>
    dplyr::count(house_size, cases) |>
    # create a complete grid of size and cases
    tidyr::complete(
      house_size = min(house_size):max(house_size), cases = 0:max(house_size)
    ) |>
    # merge NAs with 0
    dplyr::mutate(n = dplyr::if_else(is.na(n) & cases <= house_size, 0L, n)) |>
    # reshape to wide format
    tidyr::pivot_wider(names_from = cases, values_from = n) |>
    as.data.frame()
}

############################################################################
############################################################################
############################################################################

#' @title Frame apply
#' @description
#' Apply function `func` to corresponding elements
#' across a list of data.frame that share the same dimensions.
#'
#' @param frame_list `list` | `data.frame`\cr
#' List of data frames, all having the same dimensions.
#' @param func `closure`\cr
#' Function applied to corresponding elements across data frames in the list.
#'
#' @return `data.frame()`\cr
#' Data frame of the same dimensions as the input frames.
#' Columns and rownames of the first data.frame in the list.
#'
#' @keywords internal
#' @noRd
frame_apply <- function(frame_list, func) {
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

#' @title Frame apply mean
#' @description
#' Compute the mean across different dataframes.
#'
#' @param frame_list `list` | `data.frame`\cr
#'
#' @return `data.frame()`
#'
#' @keywords internal
#' @noRd
frame_apply_mean <- function(frame_list) {
  frame_apply(frame_list, function(x) {
    if (all(is.na(x))) return(NA)
    else if (all(is.numeric(x))) return(mean(x, na.rm = TRUE))
    else if (length(unique(x)) == 1) return(unique(x))
    else stop_h(
      "Undefined mean for non numeric, non unique values.",
      "Check output of your posterioir measure function."
    )
  })
}

#' @title Frame apply credible interval
#' @description
#' Compute the credible interval across different dataframes.
#'
#' @param frame_list `list` | `data.frame`\cr
#' List of data frames, all having the same dimensions.
#' @param cri_level `numeric(1)`\cr
#' List of data frames, all having the same dimensions.
#' @param upper `logical(1)`\cr
#' If `TRUE`, computes the upper credible interval. Else the lower.
#'
#' @return `data.frame()`
#'
#' @keywords internal
#' @noRd
frame_apply_cri <- function(frame_list, cri_level = 0.95, upper = TRUE) {
  p <- 0.5 + ifelse(isTRUE(upper), 0.5, -0.5) * cri_level
  frame_apply(frame_list, function(x) {
    if (all(is.na(x))) return(NA)
    else if (all(is.numeric(x))) return(stats::quantile(x, probs = p, na.rm = TRUE))
    else if (length(unique(x)) == 1) return(unique(x))
    else stop_h(
      "Undefined credible interval for non numeric, non unique values.",
      "Check output of your posterioir measure function."
    )
  })
}

############################################################################
############################################################################
########################### Housepi Posterior Expec ########################
############################################################################
############################################################################

#' @title Housepi Expectation
#' @description
#' Class to extract expecteted measures over simulations.
#'
#' @export HousepiExpectation
#' @docType class
HousepiExpectation <- R6::R6Class(
  classname = "HousepiExpectation",
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    .posterior = NULL,
    .observed = NULL,
    .cri_level = 0.95,

    # extract posterior quantities from cases-by-house-size table list
    .extract_posterior_results = function(x) {
      cri <- private$.cri_level
      return(list(
        # expected value
        expected = frame_apply_mean(x),
        # lower credible interval
        cri_lower = frame_apply_cri(x, cri, upper = FALSE),
        # upper credible interval
        cri_upper = frame_apply_cri(x, cri, upper = TRUE)
      ))
    },

    # Set evel for the credible interval.
    .set_cri_level = function(value) {
      if (is.na(value)) stop_h("Credible interval cannot be NA.")
      else if ((value <= 0.5) || (value >= 1)) stop_h(
        "Credible interval level must be a value between 0.5 and 1."
      )
      else private$.cri_level <- value
    }
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field observed `list`\cr
    #' Observed values of the measure of the posterior sample.
    observed = function(value) {
      if (missing(value)) return(private$.observed)
      stop_h("'observed' is read only.")
    },

    #' @field posterior `list`\cr
    #' Expected values of the measure of the posterior sample.
    posterior = function(value) {
      if (missing(value)) return(private$.posterior)
      stop_h("'posterior' is read only.")
    }
  ),

  ############################################################################
  ################################### PUBLIC #################################
  ############################################################################

  public = list(
    #' @noRd
    initialize = function(model, param_table, ...,
                          FUN = NULL, progress = NA, 
                          cri_level = 0.95) {
      if (is.null(FUN)) FUN <- list(...)
      if (is.null(FUN)) stop_h("No measure was provided.")

      # function to be applied also in the case of multiple measures
      apply_func <- ifelse(
        is.list(FUN),
        function(x) lapply(FUN, function(f) f(x)),
        FUN
      )

      private$.set_cri_level(cri_level)

      # extract draws without loglikelihood
      if (inherits(param_table, "HousepiPosteriorDraws")) {
        param_table <- param_table$param_draws |>
          posterior::as_draws_matrix() |>
          as.data.frame()

      } else if (posterior::is_draws(param_table)) {
        param_table <- param_table |>
          posterior::as_draws_matrix() |>
          as.data.frame()
          
      } else if (!is.data.frame(param_table)) stop_h(
        "Type of param_table not supported"
      )

      # initialize with results from observed dataset
      private$.observed <- apply_func(model$dataset)

      # deep clone model before proceeding
      # FIXME/ UPDATE: I might just need to restore...
      hmc <- model$clone(deep = TRUE)
      # cleanup at end
      withr::defer(rm(hmc))
      
      # if progress not provided, import from model
      progress <- ifelse(is.na(progress), hmc$setup$progress, progress)
      # UPDATE: can't manually set the frequency here
      if (is.numeric(progress)) progress <- TRUE

      # number of repetitions
      n_rep <- dim(param_table)[1]
      # generator seeds
      withr::with_seed(
        hmc$setup$seed,
        rng_seeds <- sample(1:max(.Machine$integer.max, n_rep),
                            n_rep, replace = FALSE)
      )

      # extract constant parameters
      aux <- hmc$parameters
      tmp <- aux[aux$const, c("name", "value")]
      constants <- as.list(tmp$value)
      names(constants) <- as.list(tmp$name)

      # define loop function working by row index
      # UPDATE: this should be improved with loops in C++ 
      loop_func <- function(r) {
        # associate parameters from table, plus constants
        hmc$load_parameter_list(
          from = c(as.list(param_table[r, ]), constants), 
          append = FALSE
        )
        # run simulation and extract cases by house size
        apply_func(hmc$simulate(seed = rng_seeds[r], progress = FALSE))
      }

      # apply loop function
      if (isTRUE(progress)) output <- lapply(
        cli::cli_progress_along(seq_len(n_rep), name = "simulations",
                                type = "tasks", clear = FALSE),
        loop_func
      )
      else output <- lapply(seq_len(n_rep), loop_func)

      if (is.list(FUN)) {
        output <- lapply(names(FUN), function(x) lapply(output, `[[`, x))
        names(output) <- names(FUN)
        private$.posterior <- lapply(output, private$.extract_posterior_results)

      } else {
        private$.posterior <- private$.extract_posterior_results(output)
      }

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info("Object of class `HousepiPosterior`.")
      invisible(self)
    }
  )
)

############################################################################
############################################################################
############################################################################

#' @title Combine observed and expected cases tables
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
#' @keywords internal
#' @noRd
combine_obs_exp <- function(obsdf, expdf) {
  # 1. curate observed dataframe
  obsdf$distr <- "Observed"
  rownames(obsdf) <- seq(1, 2 * nrow(obsdf), 2)

  # 2. curate expected dataframe
  expdf$house_size <- obsdf$house_size
  expdf$distr <- "Expected"
  rownames(expdf) <- seq(2, 2 * nrow(expdf), 2)

  # 3. combine by row bind
  combined <- rbind(obsdf, expdf)[as.character(1:(2 * nrow(obsdf))), ]
  # reset row names
  rownames(combined) <- NULL
  # reorder columns
  y <- colnames(combined)
  aux <- c("house_size", "distr")
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

# output naming convention
naming_convention_ <- function(x, cri_level) {
  qtag <- 100 * round(cri_level, 2)
  x <- gsub("expected", "mean", x) 
  x <- gsub("cri_lower", sprintf("q%.0f", 100 - qtag), x)
  x <- gsub("cri_upper", sprintf("q%.0f", qtag), x) 
  return(x)
}

#' @title Expected Cases By Size
#' @description
#' Output object of expected cases by household size, method of household model.
#'
#' @export ExpectedCasesBySize
#' @docType class
ExpectedCasesBySize <- R6::R6Class(
  classname = "ExpectedCasesBySize",
  inherit = HousepiExpectation,
  cloneable = TRUE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    .max_house_size = NA,
    .counts = NULL,

    # @description
    # Number of household with non-zero cases per household size.
    #
    # @return `data.frame`
    .count_non_zero_cases = function(dataset) {
      # UPDATE: in tidyverse it would be
      # dataset |>
      #   dplyr::group_by(house_id) |>
      #   dplyr::summarise(house_size = dplyr::n(), cases = sum(status)) |>
      #   dplyr::filter(cases > 0) |>
      #   dplyr::group_by(house_size) |>
      #   dplyr::summarise(count = dplyr::n())

      dataset$status <- as.integer(dataset$status > 0)

      agg1 <- stats::aggregate(
        status ~ house_id, data = dataset,
        FUN = function(x) c(house_size = length(x), cases = sum(x))
      )
      agg2 <- data.frame(
        house_id = agg1$house_id,
        house_size = agg1$status[, "house_size"],
        cases = agg1$status[, "cases"]
      )
      agg2 <- agg2[agg2$cases > 0, ]
      return(stats::aggregate(cases ~ house_size, data = agg2, FUN = length))
    },

    # @description
    # Return a list with all filtered case tables.
    #
    # @return `list`
    .filtered_tables = function() {
      # load observed table
      obsdf <- private$.observed$cases_by_size

      # remove rows without counts
      mask_row <- rowSums(obsdf[colnames(obsdf) != "house_size"], na.rm = TRUE) != 0

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

      expdf <- private$.posterior$cases_by_size
      return(list(
        observed = obsdf[mask_row, mask_col],
        expected = expdf$expected[mask_row, mask_col],
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
          house_size = private$.observed$sar$house_size,
          observed = private$.observed$sar$sar,
          expected = private$.posterior$sar$expected$sar,
          cri_lower = private$.posterior$sar$cri_lower$sar,
          cri_upper = private$.posterior$sar$cri_upper$sar
        )

        # associate counts
        output$counts <- private$.counts$cases
        colnames(output) <- naming_convention_(colnames(output), private$.cri_level)

        return(output)
      }
      stop_h("'sar' is read only.")
    },

    #' @field formatted_cases_table `data.frame`\cr
    #' Formatted observed and expected cases by household size.
    #' If a maximum household size is set, larger sizes are not returned.
    formatted_cases_table = function(value) {
      if (missing(value)) {
        # filter according to max_house_size
        cases <- private$.filtered_tables()

        # empty frame as a placeholder for Credible Interval
        na_df <- cases$observed
        na_df[] <- lapply(na_df, function(col) if (is.numeric(col)) NA else col)
        na_df$house_size <- cases$observed$house_size

        output <- list(expected = combine_obs_exp(cases$observed, cases$expected),
                       cri_lower = combine_obs_exp(na_df, cases$cri_lower),
                       cri_upper = combine_obs_exp(na_df, cases$cri_upper))
        names(output) <- naming_convention_(names(output), private$.cri_level)
        return(output)
      }
      stop_h("'formatted_case_table' is read only.")
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @noRd
    initialize = function(model, param_table,
                          progress = NA, cri_level = 0.95) {
      # ininitialize parent class HousepiExpectation
      super$initialize(
        model = model,
        param_table = param_table, 
        FUN = c(cases_by_size = measure_cases_by_house_size,
                sar = measure_sar_by_house_size),
        progress = progress,
        cri_level = cri_level
      )
      private$.counts <- private$.count_non_zero_cases(model$dataset)

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info("Object of class `ExpectedCasesBySize`.")
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
    #' if TRUE (default), probabilities are `P[X<=x]`, otherwise `P[X>x]`.
    #'
    #' @return
    #' `list`
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

#' @title Expected posterior measure
#' @description
#' Extract the posterior expectation of a given measure.
#' Perform simulations over given table of parameters.
#'
#' @param object `r rd_model_in_`
#' @param param_table `r rd_param_table_`
#' @param ...
#' Named arguments with respective functions.
#' If `FUN` is provided, arguments are ignored.
#' @param FUN `closure` | `list` of `closure`\cr
#' Function(s) to extract the desired measure from a dataset.
#' @param progress `r rd_progress_`
#'
#' @return
#' An object of class
#' [HousepiExpectation][housepi::HousepiExpectation].
#' 
#' @export 
expec_posterior_measure <- function(object, param_table, ...,
                                    FUN = NULL, progress = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  HousepiExpectation$new(object, param_table, ...,
                         FUN = FUN, progress = progress)
}

#' @title Expected cases by household size
#' @description
#' Extract the expected number of household cases per size.
#' It can be used to estimate the secondary attack rate.
#' Perform simulations over given table of parameters.
#'
#' @param object `r rd_model_in_`
#' @param param_table `r rd_param_table_`
#' @param progress `r rd_progress_`
#'
#' @return
#' An object of class
#' [ExpectedCasesBySize][housepi::ExpectedCasesBySize]
#' containing method to compare observed and expected cases.
#' 
#' @export
expec_cases_by_house_size <- function(object, param_table, progress = NA) {
  stopifnot(inherits(object, "HousepiModel"))
  ExpectedCasesBySize$new(object, param_table, 
                          progress = progress)
}
