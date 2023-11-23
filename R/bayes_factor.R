#!/usr/bin/env Rscript

# LIBRARIES
suppressPackageStartupMessages(library(tidyverse))

# Define a function to parse command-line arguments
parse_command_line_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  output <- list()

  # Check if there are at least two arguments
  if (length(args) != 1) {
    cat("Usage: Rscript bayes_factor.R <working_directory>\n")
    quit("no", 1)
  } else {
    output[["working_directory"]] <- args[1]
  }

  return(output)
}

# MAIN
main <- function() {
  args <- parse_command_line_args()
  working_directory <- args[["working_directory"]]
  out_file <- file.path(working_directory, "bayes_factor.csv")

  # default path for
  thermo_ladder_path <- file.path(working_directory, "thermo-integr/analysis")
  thermo_integr <- data.frame(rung = integer(0), ll_mean = numeric(0), ll_sd = numeric(0))

  # load the list of files with the MCMC of the thermodynamic ladder
  thermo_ladder_files <- list.files(thermo_ladder_path)
  num_rungs <- length(thermo_ladder_files)
  if (num_rungs == 0) {
    cat("The folder ", thermo_ladder_path,  " is empty.\n")
    quit("no", 1)
  }

  # load loglike statistics precedently computed
  for (file_name in thermo_ladder_files) {
    # Regular expression to extract rung from the filename
    r <- as.integer(gsub("infer-(\\d+)-a.csv", "\\1", file_name))

    # Add the rung filename to the dictionary
    df <- read_csv(file.path(thermo_ladder_path, file_name), show_col_types = FALSE)
    df <- df %>% column_to_rownames(., var = "key")

    # append the results
    key_res <- data.frame(rung = r, ll_mean = df["loglike", "mean"], ll_sd = df["loglike", "sd"] )
    thermo_integr <- rbind(thermo_integr, key_res)
  }

  # Reorder according to rung
  thermo_integr <- thermo_integr %>%
    arrange(.data[["rung"]]) %>%
    column_to_rownames(., var = "rung")

  # Power Posterior estimation of marginal likelihood
  ll_mean <- thermo_integr[["ll_mean"]]
  ll_var <- thermo_integr[["ll_sd"]] ^ 2

  # UPDATE : temporary choice of the rung - inverse temperature relation
  inver_temp <- (seq(0, num_rungs - 1, 1) / (num_rungs - 1))^5

  # trapezoidal estimation of the integral
  m_0 <- 1:num_rungs-1
  m_1 <- 2:num_rungs
  d_it <- inver_temp[m_1] - inver_temp[m_0]
  m_ll_mean <- (ll_mean[m_1] + ll_mean[m_0])
  d_ll_var <- (ll_var[m_1] - ll_var[m_0])
  marg_ll <- sum(d_it * m_ll_mean / 2 - (d_it^2) * d_ll_var / 12)

  # save results
  out_df <- data.frame(what = character(0), val = numeric(0))
  out_df <- rbind(out_df, data.frame(what = "num_rungs", val = num_rungs))
  out_df <- rbind(out_df, data.frame(what = "log_marg_like", val = marg_ll))

  write_csv(out_df, out_file, col_names = FALSE)
}

# execute
main()