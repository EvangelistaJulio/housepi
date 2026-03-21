## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' @importFrom rlang .data
NULL

#' @title Measure SAR by household size
#' @description
#' Retrieve secondary attack rates (SAR) from household data.
#' @param data_wide `r rd_add_data_wide_`
#' @return `data.frame`
#' @details
#' The secondary attack rate is computed for all household
#' with size \eqn{n > 1} according to
#' \eqn{\frac{\mathbb{E}[n_{+}\vert{n}] - 1}{n - 1}}
#' @export
measure_sar_by_house_size <- function(data_wide) {
  # check that data_wide is consistent
  if (length(data_wide) == 0) stop_h("No valid data_wide was provided.")
  else data_wide <- as.data.frame(data_wide)

  # cast to 0\1 values
  y <- dplyr::mutate(data_wide,
                     status = dplyr::if_else(is.na(.data$status), 0,
                                             as.integer(.data$status > 0)))
  # aggregate data by house_id
  y <- dplyr::group_by(y, .data$house_id)
  y <- dplyr::summarise(y, house_size =  dplyr::n(), cases = sum(.data$status))
  # compute sar according to definition
  y <- dplyr::filter(y, .data$cases > 0)
  y <- dplyr::mutate(y, sar = (.data$cases - 1) / (.data$house_size - 1))
  y <- dplyr::group_by(y, .data$house_size)
  y <- dplyr::summarise(y, sar = mean(.data$sar))
  as.data.frame(y)
}

#' @title Measure household cases by size
#' @description
#' Retrieve cases by household size from household data.
#' @param data_wide `r rd_add_data_wide_`
#' @return `data.frame`
#' @export
measure_cases_by_house_size <- function(data_wide) {
  # check that data_wide is consistent
  if (length(data_wide) == 0) stop_h("No valid data_wide was provided.")

  # cast to 0\1 values
  y <- dplyr::mutate(data_wide,
                     status = dplyr::if_else(is.na(.data$status), 0,
                                             as.integer(.data$status > 0)))
  # aggregate data by house_id
  y <- dplyr::group_by(y, .data$house_id)
  y <- dplyr::summarise(y, house_size = dplyr::n(), cases = sum(.data$status))
  y <- dplyr::count(y, .data$house_size, .data$cases)
  # create a complete grid of size and cases
  y <- tidyr::complete(
    y, house_size = min(.data$house_size):max(.data$house_size),
    cases = 0:max(.data$house_size)
  )
  # merge NAs with 0
  y <- dplyr::mutate(
    y, n = dplyr::if_else(is.na(.data$n) & .data$cases <= .data$house_size,
                          0L, .data$n)
  )
  # reshape to wide format
  y <- tidyr::pivot_wider(y, names_from = "cases", values_from = "n")
  as.data.frame(y)
}

############################################################################
############################################################################
######################### Housepi Posterior Prediction #####################
############################################################################
############################################################################

# Set level for the credible interval.
ppc_set_cri_level <- function(value) {
  if (is.na(value)) stop_h("Credible interval cannot be NA.")
  else if ((value <= 0.5) || (value >= 1)) stop_h(
    "Credible interval level must be a value between 0.5 and 1."
  )
  value
}

# store print messages at initialization
ppc_helper_print <- function(model_name, is_posterior, measure_names, cri) {
  sms <- c(i = sprintf("%s for model {.pkg %s}.",
                       ifelse(is_posterior, "Posterior predictions",
                              "Predictions"),
                       model_name),
           "Computed measures:")
  for (m in measure_names) sms <- c(sms, sprintf(" - {.emph %s}", m))
  c(sms, sprintf("Credible interval level set to %.0f%%.", 100 * cri))
}

#' Apply function `func` to corresponding elements
#' across a list of data.frame that share the same dimensions.
#' @param frame_list `list` of `data.frame`s. Homogeneous dataframes.
#' @param func `function`. To be applied to corresponding data elements.
#' @return A `data.frame` of the same dimensions as the input dataframes.
#' Columns and rownames of the first data.frame in the list.
#' @keywords internal
#' @noRd
ppc_frame_apply <- function(frame_list, func) {
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

  output
}

#' Compute the mean across different dataframes.
#' @param frame_list `list` of `data.frame`s.
#' @return `data.frame`
#' @keywords internal
#' @noRd
ppc_frame_apply_mean <- function(frame_list) {
  ppc_frame_apply(frame_list, function(x) {
    if (all(is.na(x))) NA_real_
    else if (all(is.numeric(x))) mean(x, na.rm = TRUE)
    else if (length(unique(x)) == 1) unique(x)
    else stop_h(
      "Undefined mean for non numeric, non unique values.",
      "Check output of your posterioir measure function."
    )
  })
}

#' Compute the credible interval across different dataframes.
#' @param frame_list `list` of `data.frame`s. Homogeneous dataframes.
#' @param cri_level `r rd_smpl_cri_level_`
#' @param upper `logical(1)`. Compute the upper credible interval.
#' If `FALSE`, compute the lower.
#' @return `data.frame`
#' @keywords internal
#' @noRd
ppc_frame_apply_cri <- function(frame_list, cri_level = 0.95, upper = TRUE) {
  p <- 0.5 + ifelse(isTRUE(upper), 0.5, -0.5) * cri_level
  ppc_frame_apply(frame_list, function(x) {
    if (all(is.na(x))) NA_real_
    else if (all(is.numeric(x))) stats::quantile(x, probs = p, na.rm = TRUE)
    else if (length(unique(x)) == 1) unique(x)
    else stop_h(
      "Undefined credible interval for non numeric, non unique values.",
      "Check output of your posterioir measure function."
    )
  })
}

# extract predicted quantities from cases-by-house-size table list
ppc_extract_predicted_results <- function(x, cri) {
  list(predicted = ppc_frame_apply_mean(x),
       cri_lower = ppc_frame_apply_cri(x, cri, upper = FALSE),
       cri_upper = ppc_frame_apply_cri(x, cri, upper = TRUE))
}

#' @title Prediction class
#' @description
#' Class to extract posterior predictive checks over measures on simulations.
#' @docType class
HousepiPrediction <- R6::R6Class(
  classname = "HousepiPrediction",
  cloneable = TRUE,
  lock_objects = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    cri_level_ = 0.95,
    stored_print_ = NULL
  ),

  ############################################################################
  ################################### PUBLIC #################################
  ############################################################################

  public = list(
    #' @description
    #' Given a model and a posterior sample of parameters,
    #' performs simulation and compute predicted values of a set of measures.
    #' @param hhmodel `r rd_model_in_`
    #' @param param_table `r rd_param_table_`
    #' @param obs_process `r rd_simul_observ_`
    #' @param obs_process_long `r rd_simul_obs_long_`
    #' @param ... `r rd_measure_dots_`
    #' @param measures_long `r rd_measure_long_`
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    #' @param cri_level `r rd_smpl_cri_level_`
    initialize = function(hhmodel,
                          param_table,
                          obs_process = NULL,
                          obs_process_long = NULL,
                          ...,
                          measures_long = NULL,
                          threads = 1,
                          progress = TRUE,
                          cri_level = 0.95) {
      # associate credible interval
      private$cri_level_ <- ppc_set_cri_level(cri_level)

      # load and classify measure functions (wide vs long)
      m_wide <- list(...)
      has_wide <- length(m_wide) > 0
      # FIXME this breaks if it's already list
      m_long <- list()
      if (!is.null(measures_long)) m_long <- list(measures_long)
      has_long <- length(m_long) > 0

      if (!has_wide && !has_long) stop_h("No measure was provided.")
      # FIXME: temporary forcing obs process long
      if (has_long && is.null(obs_process_long)) stop_h(
        "obs_process_long is required when long measures are provided."
      )

      # store (and test) results from observed data wide and/or long
      observations <- list()

      # wide-measure application
      apply_wide <- NULL
      if (has_wide) {
        observations <- c(
          observations, lapply(m_wide, function(func) func(hhmodel$data_wide))
        )
        apply_wide <- function(x) lapply(m_wide, function(func) func(x))
      }
      # long-measure application
      apply_long <- NULL
      if (has_long) {
        observations <- c(
          observations, lapply(m_long, function(func) func(hhmodel$data_long))
        )
        apply_long <- function(x) lapply(m_long, function(func) func(x))
      }

      # check measures
      if (any(lapply(observations, NROW) == 0)) stop_h(
        "One or more measures returned empty observation on data."
      )

      # run measures on posterior simulations
      output <- hhmodel$.__enclos_env__$private$run_measure_(
        param_table,
        apply_wide = apply_wide,
        apply_long = apply_long,
        obs_proc_wide = obs_process,
        obs_proc_long = obs_process_long,
        threads = threads,
        progress = progress
      )

      # parse results
      all_measure_names <- c(names(m_wide), names(m_long))

      parsed <- c()
      if (has_wide) parsed <- c(parsed,
        lapply(names(m_wide), function(x) lapply(output$wide, `[[`, x))
      )
      if (has_long) parsed <- c(parsed,
        lapply(names(m_long), function(x) lapply(output$long, `[[`, x))
      )
      names(parsed) <- all_measure_names # UPDATE: not necessary...?

      predictions <- lapply(parsed, ppc_extract_predicted_results,
                            cri = private$cri_level_)

      # change environment according to measures
      for (mn in all_measure_names) self[[mn]] <- list(
        observed = observations[[mn]],
        predicted = predictions[[mn]]$predicted,
        cri_lower = predictions[[mn]]$cri_lower,
        cri_upper = predictions[[mn]]$cri_upper
      )
      # restore locked environment
      lockEnvironment(self, bindings = TRUE)
      # store print message
      private$stored_print_ <- ppc_helper_print(
        hhmodel$name,
        inherits(param_table, "HousepiPosteriorDraws") ||
          posterior::is_draws(param_table),
        all_measure_names,
        private$cri_level_
      )

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_inform(private$stored_print_)
      invisible(self)
    }
  )
)

############################################################################
############################################################################
############################ Predicted Cases By Size #######################
############################################################################
############################################################################

# @description
# Number of household with non-zero cases per household size.
#
# @return `data.frame`
helper_count_non_zero_cases <- function(data_wide) {
  # cast to 0\1 values
  y <- dplyr::mutate(
    data_wide,
    status = dplyr::if_else(is.na(.data$status),
                            0, as.integer(.data$status > 0))
  )
  y <- dplyr::group_by(y, .data$house_id)
  y <- dplyr::summarise(y, house_size = dplyr::n(), cases = sum(.data$status))
  y <- dplyr::filter(y, .data$cases > 0)
  y <- dplyr::group_by(y, .data$house_size)
  y <- dplyr::summarise(y, count = dplyr::n())
  as.data.frame(y)
}

# @description
# Return a list with all filtered case tables.
#
# @return `list`
helper_filter_tables <- function(cases_table_by_size, max_house_size) {
  # load observed table to use as a type shape template
  obsdf <- cases_table_by_size$observed

  # remove rows without counts
  mask_row <- rowSums(obsdf[colnames(obsdf) != "house_size"], na.rm = TRUE) != 0

  # remove column with zero cases if none was observed
  mask_col <- colSums(obsdf == 0, na.rm = TRUE) != nrow(obsdf) | colnames(obsdf) != "0"

  # implement max household size mask
  if (!is.na(max_house_size)) {
    col_max <- which(colnames(obsdf) == max_house_size)
    row_max <- which(obsdf$house_size == max_house_size)

    mask_col <- mask_col & (seq_len(ncol(obsdf)) <= col_max)
    mask_row <- mask_row & (seq_along(mask_row) <= row_max & mask_row)
  }

  list(observed = obsdf[mask_row, mask_col],
       predicted = cases_table_by_size$predicted[mask_row, mask_col],
       cri_lower = cases_table_by_size$cri_lower[mask_row, mask_col],
       cri_upper = cases_table_by_size$cri_upper[mask_row, mask_col])
}

# Set naming convention of credible interval.
helper_get_name_convention <- function(x, cri) {
  qtag <- 100 * round(cri, 2)
  x <- gsub("predicted", "mean", x)
  x <- gsub("cri_lower", sprintf("q%.0f", 100 - qtag), x)
  x <- gsub("cri_upper", sprintf("q%.0f", qtag), x)
  x
}

# @description
# Function to combine cases table.
#
# @param obsdf A `data.frame` with observed cases.
# @param preddf A `data.frame` with predicted cases.
#
# @return `data.frame`
helper_combine_obs_pred <- function(obsdf, preddf) {
  # 1. curate observed dataframe
  obsdf$distr <- "Observed"
  rownames(obsdf) <- seq(1, 2 * nrow(obsdf), 2)

  # 2. curate predicted dataframe
  preddf$house_size <- obsdf$house_size
  preddf$distr <- "Predicted"
  rownames(preddf) <- seq(2, 2 * nrow(preddf), 2)

  # 3. combine by row bind
  combined <- rbind(obsdf, preddf)[as.character(1:(2 * nrow(obsdf))), ]
  # reset row names
  rownames(combined) <- NULL
  # reorder columns
  y <- colnames(combined)
  aux <- c("house_size", "distr")
  combined <- combined[, c(aux, y[!(y %in% aux)])]

  # format columns to extract labels with labelizor
  colnames(combined)[3:ncol(combined)] <- sapply(
    colnames(combined)[3:ncol(combined)],
    function(x) paste("cases", x, sep = "_")
  )

  combined
}

#' @title Predicted class with common measures of cases by size
#' @description
#' Predicted cases by household size, comparing observed and predicted values.
#' It estimates secondary attack rates and cases by size tables.
#' @docType class
PredictedCasesBySize <- R6::R6Class(
  classname = "PredictedCasesBySize",
  inherit = HousepiPrediction,
  cloneable = TRUE,
  lock_objects = FALSE,

  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################

  private = list(
    deep_clone = function(name, value) helper_clone(value),
    max_hh_size_ = NA_integer_,
    counts_ = NULL
  ),

  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################

  active = list(
    #' @field formatted_sec_attack_rate `data.frame`.
    #' Observed and predicted secondary attack rates by household size.
    #' If a maximum house size is set, larger sizes are not returned.
    formatted_sec_attack_rate = function(value) {
      if (!missing(value)) stop_h("'formatted_sec_attack_rate' is read only.")
      # convert internal list to data.frame
      output <- data.frame(
        house_size = self$sec_attack_rate$observed$house_size,
        observed = self$sec_attack_rate$observed$sar,
        predicted = self$sec_attack_rate$predicted$sar,
        cri_lower = self$sec_attack_rate$cri_lower$sar,
        cri_upper = self$sec_attack_rate$cri_upper$sar,
        counts = private$counts_$count
      )
      colnames(output) <- helper_get_name_convention(colnames(output),
                                                     private$cri_level_)
      # filter to maximum house size if applicable
      if (!is.na(private$max_hh_size_)) {
        output <- output[output$house_size <= private$max_hh_size_, ]
      }
      output
    },

    #' @field formatted_cases_table `data.frame`.
    #' Formatted observed and predicted cases by household size.
    #' If a maximum household size is set, larger sizes are not returned.
    formatted_cases_table = function(value) {
      if (!missing(value)) stop_h("'formatted_case_table' is read only.")
      # filter according to max_house_size
      cases <- helper_filter_tables(self$cases_table_by_size,
                                    private$max_hh_size_)

      # empty frame as a placeholder for Credible Interval
      na_df <- cases$observed
      na_df[] <- lapply(na_df, function(col) ifelse(is.numeric(col), NA, col))
      na_df$house_size <- cases$observed$house_size

      output <- list(
        predicted = helper_combine_obs_pred(cases$observed, cases$predicted),
        cri_lower = helper_combine_obs_pred(na_df, cases$cri_lower),
        cri_upper = helper_combine_obs_pred(na_df, cases$cri_upper)
      )
      names(output) <- helper_get_name_convention(names(output),
                                                  private$cri_level_)
      output
    }
  ),

  ############################################################################
  ################################## PUBLIC ##################################
  ############################################################################

  public = list(
    #' @param hhmodel `r rd_model_in_`
    #' @param param_table `r rd_param_table_`
    #' @param obs_process `r rd_simul_observ_`
    #' @param threads `r rd_thread_`
    #' @param progress `r rd_progress_`
    #' @param cri_level `r rd_smpl_cri_level_`
    initialize = function(hhmodel,
                          param_table,
                          obs_process = NULL,
                          threads = 1,
                          progress = TRUE,
                          cri_level = 0.95) {
      # ininitialize parent class HousepiPrediction
      super$initialize(hhmodel = hhmodel, param_table = param_table,
                       cases_table_by_size = measure_cases_by_house_size,
                       sec_attack_rate = measure_sar_by_house_size,
                       obs_process = obs_process,
                       threads = threads, progress = progress,
                       cri_level = cri_level)
      private$counts_ <- helper_count_non_zero_cases(hhmodel$data_wide)

      invisible(self)
    },

    #' @description Set maximum household size to consider in result outputs.
    #' @param n `integer(1)`. Maximum household size to consider.
    #' If `NA`, restore observed maximum.
    set_max_house_size = function(n) {
      private$max_hh_size_ <- n
      invisible(self)
    },

    #' @description Restore maximum household size to observed maximum.
    drop_max_house_size = function() {
      self$set_max_house_size(NA_integer_)
      invisible(self)
    },

    #' @description Perform Chi-squared test and return list with results.
    #' @param lower_tail
    #' `logical(1)`. Consider lower tail probabilities `P[X<=x]` (default).
    #' If `FALSE`, consider upper tail `P[X>x]`.
    #' @return `list`
    chisq_test = function(lower_tail = FALSE) {
      cases <- helper_filter_tables(self$cases_table_by_size,
                                    private$max_hh_size_)
      res <- (cases$observed - cases$predicted)^2 / cases$predicted
      x_val <- sum(res, na.rm = TRUE)
      ddof <- prod(dim(res)) - dim(res)[1] - sum(is.na(res))
      list(chisq_value = x_val, ddof = ddof,
           p_value = stats::pchisq(x_val, df = ddof, lower.tail = lower_tail))
    }
  )
)

#' @title Measures predicted by the model
#' @description
#' Extract the measures prediction according to given table of parameters.
#' Typically used to perform posterior predictive checks.
#' @param hhmodel `r rd_model_in_`
#' @param param_table `r rd_param_table_`
#' @param ... `r rd_measure_dots_`
#' @param obs_process `r rd_simul_observ_`
#' @param obs_process_long `r rd_simul_obs_long_`
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @return An object of class
#' [HousepiPrediction][housepi::HousepiPrediction].
#' @export
predicted_measures <- function(hhmodel,
                               param_table,
                               ...,
                               obs_process = NULL,
                               obs_process_long = NULL,
                               threads = 1,
                               progress = TRUE) {
  stopifnothousepi(hhmodel)
  HousepiPrediction$new(hhmodel$copy(), param_table, ...,
                        obs_process = obs_process,
                        obs_process_long = obs_process_long,
                        threads = threads,
                        progress = progress)
}

#' @title Common measures predicted by the model by household size
#' @description
#' Extract the predicted number of household cases per size.
#' It can be used to estimate the secondary attack rate (SAR).
#' @param hhmodel `r rd_model_in_`
#' @param param_table `r rd_param_table_`
#' @param obs_process `r rd_simul_observ_`
#' @param threads `r rd_thread_`
#' @param progress `r rd_progress_`
#' @return
#' An object of class
#' [PredictedCasesBySize][housepi::PredictedCasesBySize]
#' containing methods to compare observed and predicted cases.
#' @export
predicted_cases_by_house_size <- function(hhmodel,
                                          param_table,
                                          obs_process = NULL,
                                          threads = 1,
                                          progress = TRUE) {
  stopifnothousepi(hhmodel)
  PredictedCasesBySize$new(hhmodel$copy(), param_table,
                           obs_process = obs_process,
                           threads = threads,
                           progress = progress)
}

#' Format credible intervals into printable strings
#'
#' @description
#' Format elementwise mean and credible-interval bounds into a single string.
#' The function is vectorised over its arguments and will return a character
#' vector of the same length as the inputs.
#' @param mean `numeric()`. Mean values.
#' @param lower `numeric()`. Lower bounds of the credible interval.
#' @param upper `numeric()`. Upper bounds of the credible interval.
#' @param str_format `character(1)`.
#' Sprintf-style format with 3 '%s' placeholders.
#' Default is `"mean [lower, upper]"`.
#' @param as_integer `logical(1)`. Wheter values are formatted as integers.
#' Default is `FALSE`.
#' @return `character()`
#' @examples
#' format_credible_intervals(c(1.23, NA), c(0.9, NA), c(1.5, NA))
#' @export
format_credible_intervals <- function(mean, lower, upper,
                                      str_format = "%s [%s, %s]",
                                      as_integer = FALSE) {
  .fmt_vector <- function(mean, lower, upper) {
    if (length(mean) != length(lower) ||
          length(mean) != length(upper)) stop_h(
      "Mean, lower and upper values must be numeric of the same length."
    )

    if (length(gregexpr("%s", str_format, fixed = TRUE)[[1]]) != 3) stop_h(
      "Invalid choice for charachter 'str_format'.",
      "It must contain exactly three '%s' placeholders."
    )

    .fmt_num <- function(mn, low, up) {
      if (is.na(mn)) return(NA_character_)
      if (mn == low && mn == up) return(sprintf(str_format, mn, mn, mn))
      if (isTRUE(mn < low) || isTRUE(mn > up)) warning_h(
        "Found an inconsistent credible interval."
      )

      if (isTRUE(as_integer)) digits <- 0
      else if (mn == 0) digits <- max(0, ceiling(-log10(abs(low))) + 1)
      else digits <- max(0, ceiling(-log10(abs(mn))) + 1)

      mn_f <- formatC(mn, digits = digits, format = "f")
      low_f <- ifelse(is.na(low), NA_character_,
                      formatC(low, digits = digits, format = "f"))
      up_f <- ifelse(is.na(up), NA_character_,
                     formatC(up, digits = digits, format = "f"))

      sprintf(str_format, mn_f, low_f, up_f)
    }

    mapply(.fmt_num, mean, lower, upper, USE.NAMES = FALSE)
  }

  if (is.data.frame(mean) || is.matrix(mean)) {
    # check shapes
    if (!all(dim(mean) == dim(lower), dim(mean) == dim(upper))) stop_h(
      "Tabular inputs must have the same dimensions."
    )

    res <- as.data.frame(Map(.fmt_vector, as.data.frame(mean),
                             as.data.frame(lower), as.data.frame(upper)))
    names(res) <- names(mean)
    return(res)
  }
  .fmt_vector(mean, lower, upper)
}
