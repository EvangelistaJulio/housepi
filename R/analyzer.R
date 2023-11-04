#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))


# options
MIN_MCMC_VAR <- 1e-10
UNIROOT_TOL <- 1e-6
UNIROOT_MAXITER <- 1e3
PERC_BURN_IN <- 0.10 

# Define a function to parse command-line arguments
parse_command_line_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  output = list()
  # Check if there are at least two arguments
  if (length(args) == 0) {
    cat("Usage: Rscript analyzer.R <results_file> [<setup_file>] ][<params_file>]\n")
    quit("no", 1)
  } else {
    output[["results_file"]] <- args[1]

    if (length(args) > 1 ) {
      output[["setup_file"]] <- args[2]
    } else {
      output[["setup_file"]] <- NULL
    }
    if (length(args) > 2 ) {
      output[["params_file"]] <- args[3]
    } else {
      output[["params_file"]] <- NULL
    }
  }

  return(output)
}

# Main script
main <- function() {

  # Path to files
  args <- parse_command_line_args()
  file <- args[["results_file"]]
  setup_file <- args[["setup_file"]]
  params_file <- args[["params_file"]]

  out_path <- file.path(dirname(file), "analysis")
  dir.create(out_path, showWarnings = FALSE)
  out_file <- file.path(out_path, tools::file_path_sans_ext(basename(file)))
  out_file <- paste0(out_file, "-a.", tools::file_ext(file))

  # Load setup 
  if (is.null(setup_file)) {
    delim <- ","
    num_adapt_iter <- 0
    burn_in <- 0
  } else if (file.exists(setup_file)) {
    setup <- read_csv(setup_file, col_names = FALSE, show_col_types = FALSE)
    setup <- setup %>% column_to_rownames(., var = "X1")
    delim <- as.character(setup["delim", ])
    num_iter <- as.integer(setup["num_iter", ])
    num_adapt_iter <- as.integer(setup["num_adapt_iter", ])
    burn_in <- floor(PERC_BURN_IN * (num_iter - num_adapt_iter))
  } else {
    cat("Setup file ", setup_file, " not found.\n")
    quit("no", 1)
  }

  # Load results 
  if (file.exists(file)) {
    results <- read_delim(file, delim = delim, show_col_types = FALSE)
  } else {
    cat("File ", file, " not found.\n")
    quit("no", 1)
  }

  # Load parameters 
  if (is.null(params_file)) {
    params <- tibble()
  } else if (file.exists(params_file)) {
    params <- read_delim(params_file, delim = delim, col_names = FALSE, show_col_types = FALSE)
    params <- params %>% column_to_rownames(., var = "X1")
  } else {
    cat("Parameters file ", params_file, " not found.\n")
    quit("no", 1)
  }

  # Define working data frame
  # extract adaptive steps and burn_in
  df <- results %>% filter(.data[["iter"]] > num_adapt_iter + burn_in)
  thinning = df[["iter"]][2] - df[["iter"]][1]

  sample_size <- length(df[["iter"]])
  # find columns of interest 
  cols <- colnames(df)
  cols <- cols[cols != "iter"]
  cols <- cols[!grepl("__a", cols)]

  # compute mean and dev.std from the MCMC chains
  # UPDATE : this can be more compact in tidyverse piping
  output_df <- data.frame()
  for (p in cols) {
    signal <- as.vector(df[[p]])

    if (var(signal) < MIN_MCMC_VAR) {
      # pathological case: signal is constant
      key_data <- data.frame(key = p, mean = mean(signal), sd = 0, period = sample_size * thinning, ESS = 1)
    } else {
      # compute autocorrelation function
      acorr <- acf(signal - mean(signal), lag.max=floor(length(signal) / 2), plot=FALSE)[["acf"]]
      # linear inerpolation
      interp_func <- approxfun(0:(length(acorr) - 1), acorr, method="linear") 

      # find first root
      x <- seq(0, length(acorr) - 1, 1)
      y <- interp_func(x)
      end <- x[which(y < 0)[1]]

      if (end < length(acorr) - 1) {
        result <- uniroot(interp_func, c(0, end), tol=UNIROOT_TOL, maxiter=UNIROOT_MAXITER, check.conv=TRUE)
        period <- as.integer(round(result[["root"]]))

      } else {
        # pathological case
        # FIXME : what about the (unphysical) case where acorr converges asymptotically to 0?
        period <- length(acorr) - 1
      }

      # mask for effective sample
      mask <- seq(1, sample_size, period)
      # effective sample size
      ess <- length(mask)
      # compute statistics on effective sample
      mu <- mean(signal[mask])
      sigma <- sd(signal[mask])

      # statistically relevant num. of digits
      r <- max(0, ceiling(1 - log10(sigma)))

      key_data <- data.frame(key = p, mean = round(mu, r), sd = round(sigma, r), period = period * thinning, ESS = ess )
      output_df <- rbind(output_df, key_data)
    }
  }

  # save computations
  write_delim(output_df, out_file, delim=delim )

}

# Execute
main()
