#
# analysis.R
#
# Copyright (C) 2025 CAMAGLIA Francesco - Institut Pasteur, Paris (France).
#

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
  # perform element-wise operations
  output <- matrix(
    sapply(seq_along(frame_list[[1]]), function(x) {
      apply(sapply(frame_list, "[[", x), 1, func)
    }),
    nrow = nrow(frame_list[[1]]),
    byrow = FALSE
  )
  colnames(output) <- colnames(frame_list[[1]])
  rownames(output) <- rownames(frame_list[[1]])
  return(as.data.frame(output))
}

#' @title
#' Retrieve Household Cases by Size
#' @description
#' Function to retrieve cases by size from a given data.frame.
#'
#' @param dataset `data.frame`\cr
#' Dataset containing the the individual 'id'
#' their 'status' and their 'house_id'.
#'
#' @return `data.frame`\cr
#'
#' @noRd
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

#' @title
#' Retrieve SAR from Cases by Size
#' @description
#' Function to retrieve secondary attack rates from data.frame
#' with cases by size.
#'
#' @param dataset `data.frame`\cr
#' Dataset containing the the individual 'id'
#' their 'status' and their 'house_id'.
#'
#' @return `data.frame`\cr
#'
#' @noRd
retrieve_sar_from_casesbysize <- function(cases_df) {
  # drop size column
  cases_df$size <- NULL
  # drop zero size column
  cases_df$`0` <- NULL
  cases_df <- cases_df[as.integer(row.names(cases_df)) > 1, ]

  # total n. cases (row names stay stored)
  n <- rowSums(cases_df, na.rm = TRUE)
  # normalization factors
  norm <- 1. / (as.vector(n) * (as.integer(names(n)) - 1))

  # extract vector of secondary cases
  sec_cases <- as.integer(colnames(cases_df)) - 1

  # substitute NAs with 0s
  cases_df[is.na(cases_df)] <- 0
  output <- data.frame(
    sar = as.vector(as.matrix(cases_df) %*% sec_cases) * norm,
    row.names = names(n)
  )

  return(output)
}

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
############################################################################

#' @title
#' Housepi Output Expected Cases
#' @description
#' Output object from expected cases by size.
#'
#' @keywords internal
ExpectedCases <- R6::R6Class(
  classname = "ExpectedCases",
  cloneable = TRUE,

  ############################################################################
  ############################################################################
  ################################## PRIVATE #################################
  ############################################################################
  ############################################################################

  private = list(
    .case_table = NULL,
    .sar = NULL,
    .cri_level = 0.95,
    .max_house_size = NA,

    # @description
    # Return a list with all filtered case tables.
    #
    # @return
    # `list()`)
    .filtered_tables = function() {
      # load observed table
      obsdf <- private$.case_table$observed

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

      return(list(
        observed = obsdf[mask_row, mask_col],
        expected = private$.case_table$expected[mask_row, mask_col],
        cri_lower = private$.case_table$cri_lower[mask_row, mask_col],
        cri_upper = private$.case_table$cri_upper[mask_row, mask_col]
      ))
    },

    # @description
    # Extract expected quantities from cases-by-size table list
    #
    # @param `list` of `data.frame`\cr
    # List of data.frame containing cases-by-size.
    .expected_from_casesbysize_list = function(y) {
      # extract list of secondary attack rates
      sar_list <- lapply(y, retrieve_sar_from_casesbysize)

      # compute mean values
      f_mean <- function(x) mean(x, na.rm = TRUE)
      private$.sar$expected <- frame_apply(f_mean, sar_list)
      private$.case_table$expected <- frame_apply(f_mean, y)

      # internal consistency check
      # note: normally it shouldn't be possible...
      if (any(dim(private$.case_table$observed) != dim(private$.case_table$expected))) {
        stop_h("Internal inconsistency in observed and expected tables.")
      }

      # compute lower credible interval
      f_cri_lower <- function(x) {
        stats::quantile(x, probs = 0.5 * (1 - private$.cri_level), na.rm = TRUE)
      }
      private$.sar$cri_lower <- frame_apply(f_cri_lower, sar_list)
      private$.case_table$cri_lower <- frame_apply(f_cri_lower, y)

      # compute upper credible interval
      f_cri_upper <- function(x) {
        stats::quantile(x, probs = 0.5 * (1 + private$.cri_level), na.rm = TRUE)
      }
      private$.sar$cri_upper <- frame_apply(f_cri_upper, sar_list)
      private$.case_table$cri_upper <- frame_apply(f_cri_upper, y)

      invisible(self)
    }
  ),

  ############################################################################
  ############################################################################
  ################################## ACTIVE ##################################
  ############################################################################
  ############################################################################

  active = list(
    #' @field sar `data.frame`\cr
    #' Observed and expected secondary attack rates by household size.
    #' If a maximum house size is set, larger sizes are not returned.
    sar = function(value) {
      if (missing(value)) {
        # convert internal list to data.frame
        output <- as.data.frame(private$.sar)
        colnames(output) <- names(private$.sar)

        output$house_size <- as.integer(row.names(output))
        if (!is.na(private$.max_house_size)) {
          output <- output[output$house_size <= private$.max_house_size, ]
        }

        n <- private$.case_table$observed[row.names(output), ]
        n$size <- NULL
        n$`0` <- NULL
        output$counts <- as.vector(rowSums(n, na.rm = TRUE))

        return(output)
      }
      stop_h("'sar' is read only.")
    },

    #' @field cases_table `data.frame`\cr
    #' Observed and expected cases by household size.
    #' If a maximum house size is set, larger sizes are not returned.
    cases_table = function(value) {
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
      stop_h("'cases_table' is read only.")
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
    #' [ExpectedCases][housepi::ExpectedCases].
    #'
    #' @param house_model `HouseholdModel`\cr
    #' Object of class
    #' [HouseholdModel][housepi::HouseholdModel]
    #' containing all model definitions after having run `$infer`.
    #' @param cri_level `double(1)`\cr
    #' Probability for credible interval (default is 0.95).
    #' @param progress `logical(1)`\cr
    #' Progress of simulations. Default is inherited from the model.
    initialize = function(house_model,
                          cri_level = 0.95,
                          progress = NA) {
      # check credible interval
      if ((cri_level <= 0) || (cri_level >= 1)) {
        stop_h(
          "Credible interval level must be a value between 0 and 1."
        )
      }
      private$.cri_level <- cri_level

      # initialize with observed case table
      obs <- retrieve_cases_by_size(house_model$dataset)
      private$.sar <- list(observed = retrieve_sar_from_casesbysize(obs))
      private$.case_table <- list(observed = obs)

      # extract traces without loglikelihood
      traces <- house_model$traces_no_loglike
      # stuck multiple chains
      if (is.list(traces)) traces <- do.call(rbind, traces)

      # number of repetitions
      n_rep <- dim(traces)[1]

      # deep clone model before proceeding
      # FIXME/UPDATE/NOTE: deep clone is not cloning traces e. g.!
      hmc <- house_model$clone(deep = TRUE)
      # cleanup at end
      withr::defer(rm(hmc))

      if (is.na(progress)) {
        progress <- isTRUE(hmc$setup$progress)
      } else {
        progress <- isTRUE(progress)
      }

      # generator seeds
      withr::with_seed(
        hmc$setup$seed,
        rng_seeds <- sample(1:max(1e8, n_rep), n_rep, replace = FALSE)
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
        return(retrieve_cases_by_size(
          hmc$simulate(seed = rng_seeds[r], progress = FALSE)
        ))
      }

      if (isTRUE(progress)) {
        results <- lapply(
          cli::cli_progress_along(
            seq_len(n_rep), name = "simulations",
            type = "tasks", clear = FALSE
          ), loop_func
        )
      } else {
        results <- lapply(seq_len(n_rep), loop_func)
      }
      # extract expected quantities from cases-by-size table list
      private$.expected_from_casesbysize_list(results)

      invisible(self)
    },

    #' @description
    #' Override the default printing behavior.
    print = function() {
      cli::cli_alert_info(
        "Output of method `$expec_cases_by_size` of HouseholdModel."
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
    #' @param lower.tail `logical(1)`\cr
    #' if TRUE (default), probabilities are P[X<=x], otherwise P[X>x].
    #'
    #' @return
    #' `list()`)
    chisq_test = function(lower.tail = FALSE) {
      x <- private$.filtered_tables()
      res <- (x$observed - x$expected)^2 / x$expected
      chisq_value <- sum(res, na.rm = TRUE)
      ddof <- prod(dim(res)) - dim(res)[1] - sum(is.na(res))
      return(list(
        chisq_value = chisq_value, ddof = ddof,
        p_value = stats::pchisq(chisq_value, df = ddof, lower.tail = lower.tail)
      ))
    },

    #' @description
    #' Plot observed vs expected secondary attack rates.
    #'
    #' @return
    #' `ggplot` object)
    plot_sar = function() {
      if (require_suggested_library("ggplot2")) {
        col_dict <- c("Expected" = "coral", "Observed" = "black")
        g <- ggplot2::ggplot(self$sar, ggplot2::aes(x = house_size)) +
          # expected SAR
          ggplot2::geom_point(ggplot2::aes(y = expected, color = "Expected")) +
          ggplot2::geom_errorbar(ggplot2::aes(
            y = expected, ymin = cri_lower, ymax = cri_upper, color = "Expected"
          ), width = 0.2) +
          # observed SAR
          ggplot2::geom_point(ggplot2::aes(y = observed, color = "Observed")) +
          # text with observed counts
          ggplot2::geom_text(ggplot2::aes(
            x = house_size + 0.25, y = observed, label = sprintf("[%s]", counts)
          ),
          color = col_dict["Observed"]) +
          ggplot2::labs(x = "Household Size",
                        y = "Secondary Attack Rate",
                        color = "Distribution") +
          ggplot2::scale_color_manual(values = col_dict)
        return(g)
      }
    },

    #' @description
    #' Return a `flextable` combining observed and expected cases by size.
    #'
    #' @param show_chisq `\logical(1)`\cr
    #' Whether to show footer with Chi-squared test.
    #' @param cri_format `character(1)`\cr
    #' Formatting string for the credible interval passed to `sprintf`.
    #' It expects 3 double values (mean, CRI lower and CRI upper).
    #'
    #' @return
    #' `flextable`
    show_cases_table = function(show_chisq = FALSE,
                                cri_format = "%.0f [%.0f, %.0f]") {
      # load cases-table
      cases_df <- self$cases_table
      combined <- cases_df$expected

      # select expected cases
      mask_exp <- cases_df$expected$distr == "Expected"
      mask_cols <- setdiff(names(cases_df$expected), c("size", "distr"))
      # format `mean` with its `lower` and `upper` values
      format_exp <- function(mean, lower, upper) {
        return(ifelse(is.na(mean) | is.na(lower) | is.na(upper),
                      NA, sprintf(cri_format, mean, lower, upper)))
      }
      combined[mask_exp, mask_cols] <- as.data.frame(mapply(
        format_exp,
        cases_df$expected[mask_exp, mask_cols],
        cases_df$cri_lower[mask_exp, mask_cols],
        cases_df$cri_upper[mask_exp, mask_cols]
      ))

      if (require_suggested_library("flextable")) {
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
        ft <- flextable::autofit(ft, add_w = 0, add_h = 0)
        ft <- flextable::colformat_num(ft, big.mark = "",
                                       na_str = "", nan_str = "")

        if (isTRUE(show_chisq)) {
          x <- self$chisq_test()
          footer <- sprintf("X-squared = %.1f, df = %.0f, p-value = %f",
                            x$chisq_value, x$ddof, x$p_value)
          ft <- flextable::add_footer_lines(ft, footer)
        }

        return(ft)
      }
    }
  )
)
