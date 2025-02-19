#
# analysis.R
#
# Copyright (C) 2025 Francesco Camaglia - Institut Pasteur, Paris (France).
#

#' @title
#' Extract results from coda
#' @description
#' Collects relevant statistics from the the input mcmc traces.
#' 
#' @param trs
#' An object of class `mcmc` or `mcmc.list`.
#' 
#' @return
#' (`data.frame`)
#' 
#' @export
retrieve_traces_results <- function(trs) {
  if (!(methods::is(trs, "mcmc") || methods::is(trs, "mcmc.list"))) {
    stop.h(
      "Input should be a `coda::mcmc` or a `coda::mcmc.list` object."
    )
  }

  # note: coda is changing format for single columns...
  smr <- summary(trs)
  if (is.null(dim(smr$statistics))) {
    output <- data.frame(
      "mean" = smr$statistics["Mean"],
      "sd" = smr$statistics["SD"],
      "q2.5%" = smr$quantiles["2.5%"],
      "q97.5%" = smr$quantiles["97.5%"], 
      check.names = FALSE
    )
    row.names(output) <- colnames(trs)
  } else {
    output <- data.frame(
      "mean" = smr$statistics[, "Mean"],
      "sd" = smr$statistics[, "SD"],
      "q2.5%" = smr$quantiles[, "2.5%"],
      "q97.5%" = smr$quantiles[, "97.5%"], 
      check.names = FALSE
    )
    row.names(output) <- row.names(smr$statistics)
  }

  return(output)
}

#' @title
#' Apply on Homogeneous Frames
#' @description
#' Apply function `FUN` to corresponding elements 
#' across a list of data frames that share the same dimensions.
#'
#' @param FUN (`closure()`)\cr
#' Function to be applied element-wise across corresponding positions in the data frames.
#' @param frame_list (`list()` | `data.frame`)\cr
#' List of data frames, all having the same dimensions.
#'
#' @return (`data.frame()`)\cr
#' Data frame of the same dimensions as the input frames.
#' 
apply_on_homogeneous_frames <- function (FUN, frame_list) {
  # perform element-wise operations
  output <- matrix(
    sapply(seq_along(frame_list[[1]]), function(x) {
      apply(sapply(frame_list, "[[", x), 1, FUN)
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
#' Dataset containing the the individual 'id', their 'status' and their 'house_id'.
#' 
#' @return (`data.frame`)\cr
#' 
#' @export
retrieve_household_cases_by_size <- function(dataset) {
  # check that dataset is consistent
  if (length(dataset) == 0) {
    stop.h("No valid dataset was provided.")
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
  merged_data <- merge(grid, cases_by_size, by = c("size", "cases"), all.x = TRUE)
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
#' Retrieve Expected Cases by Size
#' @description
#' For a given `house_model` simulate different household datasets
#' according to parameters sampled from the model posterior distribution.
#' 
#' @param house_model (`HouseholdModelClass`)\cr
#' Object of class \[HouseholdModelClass][housepi::HouseholdModelClass]
#' containing all model definitions.
#' @param dataset (`character(1)` | `data.frame`)\cr
#' "Path/to/file.csv" or `data.frame` containing the dataset.
#' @param size (`integer(1)`)\cr
#' Sample size of a single chain according to posterior distribution.
#' @param n_burnin (`integer(1)`)\cr
#' Number of steps to discard at the beginning of the sample.
#' @param thin (`integer(1)`)\cr
#' Number of steps to discard between each sample.
#' @param n_adapt (`integer(1)`)\cr
#' Number of steps with finite adaptation before sampling.
#' @param n_chains (`integer(1)`)\cr
#' Number of independent MCMC chains.
#' @param n_threads (`integer(1)`)\cr
#' Number of threads for the parallel computations.
#' @param rand_init (`logical(1)`)\cr
#' If `TRUE` the initial value of the parameters is drawn from prior.
#' @param seed (`integer(1)`)\cr
#' Seed of the random number generator.
#' If `NA` uses the internal clock (default).
#' @param progress (`logical(1)`)\cr
#' If `FALSE` suppress prints about progress state.
#' @param max_house_size (`integer(1`)\cr
#' Maximum household size to consider.
#' 
#' @export
retrieve_expected_cases_by_size <- function(house_model,
                                            dataset = NULL,
                                            size = 100,
                                            n_burnin = 0,
                                            thin = 1,
                                            n_adapt = 0,
                                            n_chains = 1,
                                            n_threads = 1,
                                            rand_init = TRUE,
                                            seed = NA,
                                            progress = TRUE,
                                            max_house_size = NA) {
  # FIXME deep clone not working
  # clone the model before proceeding further
  # hmc <- house_model$clone(deep = TRUE) 
  # cleanup at end
  #withr::defer(rm(hmc))
  hmc <- house_model
  
  # generator seeds
  if (is.na(seed)) seed <- sample(1:max(1e8), 1)
  withr::with_seed(
    seed,
    rng_seeds <- sample(1:max(1e8, size), 
                        size + 1, replace = FALSE)
  )
  
  # associate dataset through 'infer()' and then perform sampling
  hmc$infer(dataset,
    size = size, n_burnin = n_burnin, thin = thin, 
    n_adapt = n_adapt, n_chains = n_chains,
    n_threads = n_threads, rand_init = rand_init, 
    progress = progress, seed = rng_seeds[size + 1])
  
  # extract traces
  traces <- hmc$traces_no_loglike
  
  # extract constant parameters
  params <- hmc$param_defs
  tmp <- params[params$constant, c("name", "value")]
  constants <- as.list(tmp$value)
  names(constants) <- as.list(tmp$name)
  
  # store simulated datasets
  if (isTRUE(progress)) {
    cli::cli_progress_bar("simulations", total = n_repetitions, 
                          type = "tasks", clear = FALSE)
  }
  simulations_list <- list()
  
  for (r in 1:n_repetitions) {
    # extract parameters from posterior sample
    hmc$load_param_list(c(as.list(traces[r,]), constants), append = FALSE)
    
    # extract cases by size tables
    # UPDATE: remove `datasete` once fixed pre-loading issue
    hmc$simulate(dataset, seed = rng_seeds[r], progress = FALSE)
    simulations_list[[r]] <- retrieve_household_cases_by_size(hmc$dataset)
    
    if (isTRUE(progress)) cli::cli_progress_update()
  }
  if (isTRUE(progress)) cli::cli_progress_done()
  
  # extract observed household case table
  obs_hh_case_table <- retrieve_household_cases_by_size(dataset)
  
  # keep only observed house sizes and discard according to maximum
  a <- obs_hh_case_table
  m <- apply(a[colnames(a) != "size"], 1, 
                function(x) { !all(x == 0, na.rm = TRUE) })
  if (!is.na(max_house_size)) m <- m & (a$size <= max_house_size)
  obs_hh_case_table <- a[m, ]
  
  # compute mean values
  exp_hh_case_table <- apply_on_homogeneous_frames(
    function(x) { mean(x, na.rm = TRUE) }, 
    simulations_list
  )[m, ]
  
  # consistency check
  if (all(dim(obs_hh_case_table) != dim(exp_hh_case_table))) {
    stop.h(
      "Internal inconsistency in observed and expected tables."
    )
  }
  
  # perform Pearson chi-squared test
  res <- (obs_hh_case_table - exp_hh_case_table)^2 / exp_hh_case_table
  chisq_value <- sum(res, na.rm = TRUE)
  ddof <- prod(dim(res)) - dim(res)[1] - sum(is.na(res))
  p_value <- stats::pchisq(chisq_value, df = ddof, lower.tail = FALSE)
  
  # create a S3 object for the output
  output <- structure(
    list(
      "traces" = traces,
      "observed_cases_by_size" = obs_hh_case_table,
      "expected_cases_by_size" = exp_hh_case_table,
      "Q2.5" = apply_on_homogeneous_frames(
        function(x) {stats::quantile(x, probs = 0.025, na.rm = TRUE)}, 
        simulations_list)[m, ],
      "Q97.5" = apply_on_homogeneous_frames(
        function(x) {stats::quantile(x, probs = 0.975, na.rm = TRUE)},
        simulations_list)[m, ],
      "X_squared" = chisq_value,
      "df" = ddof,
      "p_value" = p_value,
      "seed" = seed
    ),
    class = "housepi_model_check"
  )
  
  return(output)
}

#' @title
#' Print S3 method for housepi_model_check objects
#' @description
#' Provide a custom print method for objects of class `housepi_model_check`.
#' 
#' @param x 
#' An object of class `housepi_model_check`.
#' 
#' @export
print.housepi_model_check <- function(x, ...) {
  cat("Household model check with Chi-squared test:\n")
  cat(sprintf("X-squared = %.1f, df = %s, p-value = %s", 
              x$X_squared, x$df, x$p_value))
}

#' @title
#' Show model check table
#' @description
#' Return a flextable combining the results of the check on household model
#' 
#' @param x 
#' An object of class `housepi_model_check`.
#' 
#' @return 
#' (`flextable` | `data.frame`)
#' 
#' @export
show_model_check_table <- function(x) {
  
  if (!methods::is(x, "housepi_model_check")) {
    stop.h(
      "Input should be the result of function 'check_household_model'."
    )
  }
  
  # curate observed dataframe
  A <- as.data.frame(apply(x$observed_cases_by_size, 2, function(x) {
    as.character(as.integer(x, na.rm = T), na.rm = T)}))
  A$distr <- "Observed"
  rownames(A) <- seq(1, 2*nrow(A), 2)
  
  # curate expected dataframe

  # format `mean` with its `lower` and `upper` values
  format_expression <- function(mean, lower, upper) {
    return(ifelse(is.na(mean) | is.na(lower) | is.na(upper),
                  NA, sprintf("%.1f[%.0f-%.0f]", mean, lower, upper)))
  }
  B <- as.data.frame(mapply(format_expression, 
                            x$expected_cases_by_size, x$Q2.5, x$Q97.5))
  B$size <- A$size
  B$distr <- "Expected"
  rownames(B) <- seq(2, 2*nrow(B), 2)
  
  # combine data frames
  combined <- rbind(A, B)[as.character(1:(2*nrow(A))), ]
  
  # reset row names
  rownames(combined) <- NULL
  
  # reorder columns
  y <- colnames(combined)
  msk <- c("size", "distr")
  combined <- combined[, c(msk, y[!(y %in% msk)])]
  colnames(combined)[1] <- "Household size"
  colnames(combined)[2] <- "Distribution"

  if (requireNamespace("flextable", quietly = TRUE)) {
    # format output with flextable
    ft <- flextable::flextable(combined)
    ft <- flextable::bold(ft, ~ `Household size`, bold = TRUE)
    ft <- flextable::theme_vanilla(ft)
    ft <- flextable::merge_v(ft, j = ~ `Household size`)
    ft <- flextable::add_header_row(ft, values = c("", "Number of cases"), colwidths = c(2, ncol(combined) - 2))
    ft <- flextable::set_table_properties(ft, align = "center", layout = "autofit")
    ft <- flextable::add_footer_lines(ft, sprintf("X-squared = %.1f, df = %s, p-value = %s", x$X_squared, x$df, x$p_value))

    return(ft)

  } else {
    return(combined)
  }
}
