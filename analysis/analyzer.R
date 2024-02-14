#!/usr/bin/env Rscript

# clean environment
rm(list = ls())

# LIBRARIES
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales)) # for hue_pal
suppressPackageStartupMessages(library(ggcorrplot)) # for ggcorrplo
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(optparse))
options(warn = -1)

# FORMAT : options
delim <- ","
dict_delim <- "/"

# ANALYSIS : options
perc_burn_in <- 0.10
min_mcmc_var <- 1e-10
uniroot_tol <- 1e-6
uniroot_maxiter <- 1e3
acorr_n_period <- 5
num_hist_bins <- 22

# FIGURE : options
clr_stat <- "black"
clr_prior <- "grey"
clr_target <- "red"
clr_heatmap <- c("#297db8", "white", "#e55b11")
# Define a custom theme and set it as the default theme for all plots
custom_theme <- theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA)
  )
theme_set(custom_theme)

# Define a function to parse command-line arguments
parse_command_line_args <- function() {
  if (length(commandArgs(trailingOnly = TRUE)) == 0) {
    cat("Usage: Rscript analyzer.R -i <file_to_analyze> [options]\n")
    quit("no", 1)
  } else {
    option_list <- list(
      make_option(c("-i", "--input_file"), type = "character", help = "Path/to/file to analyze"),
      make_option(c("-o", "--output_file"), type = "character", default = NULL, help = "Path/to/file where to store analysis"),
      make_option(c("-s", "--setup_file"), type = "character", default = NULL, help = "Path/to/file with setup"),
      make_option(c("-p", "--params_file"), type = "character", default = NULL, help = "Path/to/file with parameters"),
      make_option(c("-f", "--figure_file"), type = "character", default = NULL, help = "Path/to/file with figure"),
      make_option(c("--save_fig"), action = "store_true", default = FALSE, help = "Save figure to default"),
      make_option(c("--is_target"), action = "store_true", default = FALSE, help = "Parameters are considered target")
    )

    # Parse the command-line arguments
    return(parse_args(OptionParser(option_list = option_list)))
  }
}

# MAIN
main <- function() {
  # Load options
  args <- parse_command_line_args()
  input_file <- args[["input_file"]]
  output_file <- args[["output_file"]]
  setup_file <- args[["setup_file"]]
  params_file <- args[["params_file"]]
  figure_file <- args[["figure_file"]]
  save_fig <- args[["save_fig"]]
  is_target <- args[["is_target"]]

  # Load setup
  if (is.null(setup_file)) {
    burn_in <- 0
  } else if (file.exists(setup_file)) {
    setup <- read_csv(setup_file, col_names = FALSE, show_col_types = FALSE)
    setup <- setup %>% column_to_rownames(., var = "X1")
    burn_in <- floor(perc_burn_in * as.integer(setup["num_iter", ]))
  } else {
    cat("Setup file ", setup_file, " not found.\n")
    quit("no", 1)
  }

  # Load input file
  if (file.exists(input_file)) {
    input_df <- read_delim(input_file, delim = delim, show_col_types = FALSE)
  } else {
    cat("File ", input_file, " not found.\n")
    quit("no", 1)
  }

  # Load output
  if (is.null(output_file)) {
    # default
    out_path <- file.path(dirname(input_file), "analysis")
    dir.create(out_path, showWarnings = FALSE)
    output_file <- file.path(
      out_path,
      tools::file_path_sans_ext(basename(input_file))
    )
    output_file <- paste0(output_file, "-a.", tools::file_ext(input_file))
  } else if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file))
  }
  output_file_corr <- paste0(tools::file_path_sans_ext(output_file), "-corr.", tools::file_ext(output_file))

  # Load figure
  if (is.null(figure_file)) {
    if (save_fig == TRUE) {
      # default
      pdf_path <- file.path(dirname(input_file), "figures")
      dir.create(pdf_path, showWarnings = FALSE)
      figure_file <- file.path(
        pdf_path,
        tools::file_path_sans_ext(basename(input_file))
      )
      figure_file <- paste0(figure_file, "-fig.pdf")
    }
  } else {
    save_fig <- TRUE
    if (!dir.exists(dirname(figure_file))) {
      dir.create(dirname(figure_file))
    }
  }
  if (save_fig == TRUE) {
    figure_file_corr <- paste0(tools::file_path_sans_ext(figure_file), "-corr.", tools::file_ext(figure_file))
  }

  # Load parameters
  if (is.null(params_file)) {
    params <- tibble()
    is_target <- FALSE
  } else if (file.exists(params_file)) {
    params <- read_delim(
      params_file,
      delim = delim,
      col_names = FALSE,
      show_col_types = FALSE
    )
    params <- params %>% column_to_rownames(., var = "X1")
  } else {
    cat("Parameters file ", params_file, " not found.\n")
    quit("no", 1)
  }

  # display a short summary
  cat("Performing analysis on file : ", input_file, "\n")
  cat("Results are saved on file : ", output_file, "\n")
  if (save_fig == TRUE) {
    cat("Figure is saved on file : ", figure_file, "\n")
  }

  # define columns of interest
  selec_cols <- input_df %>%
    select(-iter) %>%
    select(-contains("__a")) %>%
    colnames()

  # initialize output
  output_df <- data.frame()
  if (save_fig == TRUE) {
    plot_list <- list()
    colors <- hue_pal()(length(selec_cols))
  }

  # define working data frame and remove burn_in
  df <- input_df %>%
    filter(iter > burn_in)

  # Correlation Matrix
  corr_mtx <- df %>%
    select(all_of(selec_cols)) %>%
    select(-loglike) %>%
    cor()

  # save computations
  write_delim(
    data.frame(corr_mtx),
    output_file_corr,
    delim = delim
  )

  if (save_fig == TRUE) {
    corr_mtx_plot <- ggcorrplot(
      corr_mtx,
      hc.order = TRUE,
      outline.col = "white",
      colors = clr_heatmap
    )

    ggsave(
      figure_file_corr,
      plot = corr_mtx_plot,
      device = "pdf",
      dpi = 300,
      width = length(selec_cols),
      height = length(selec_cols)
    )
  }


  sample_size <- length(df[["iter"]])
  # extract thinning
  thinning <- df[["iter"]][2] - df[["iter"]][1]

  # compute mean and dev.std for selected MCMC chains
  for (i in seq(1, length(selec_cols), 1)) {
    p <- selec_cols[i]
    signal <- as.vector(df[[p]])

    if (var(signal) < min_mcmc_var) {
      # pathological case: signal is constant
      # UPDATE : find a way to not block the code
      cat("The variance of `", p, "` is lower than minimum MCMC variance ", min_mcmc_var, "\n")
      quit("no", 1)
    }

    # UPDATE : apprently this can be more compact in tidyverse piping
    # compute autocorrelation function

    # se non funziona è colpa di GitHub Copilot
    #acorr <- acf(
    #  signal - mean(signal),
    #  lag.max = floor(0.5 * length(signal)),
    #  plot = FALSE
    #)[["acf"]]

    acorr <- signal %>%
      subtract(mean(signal)) %>%
      acf(lag.max = floor(0.5 * length(signal)), plot = FALSE) %>%
      extract("acf")

    x <- seq(0, length(acorr) - 1, 1)
    # linear inerpolation
    interp_func <- approxfun(x, acorr, method = "linear")

    # find first root
    end <- x[which(acorr < 0)[1]]
    if (end < length(acorr) - 1) {
      x_0 <- uniroot(
        interp_func,
        c(0, end),
        tol = uniroot_tol,
        maxiter = uniroot_maxiter,
        check.conv = TRUE
      )
      period <- as.integer(round(x_0[["root"]]))
    } else {
      # pathological case : acorr converges asymptotically to 0
      period <- length(acorr) - 1
    }

    # mask for effective sample
    mask <- seq(1, sample_size, period)
    # effective sample size
    ess <- length(mask)
    # compute statistics on effective sample
    mu <- mean(signal[mask])
    sigma <- sd(signal[mask])

    # append results to the output dataframe
    key_data <- data.frame(
      key = p,
      mean = mu,
      sd = sigma,
      period = period * thinning,
      ESS = ess
    )
    output_df <- rbind(output_df, key_data)

    # plot summary of the analysis
    if (save_fig == TRUE) {
      # options
      min_val <- min(df[[p]])
      max_val <- max(df[[p]])
      binwidth <- (max_val - min_val) / num_hist_bins

      # parameters
      show_prior <- (p %in% rownames(params))
      if (show_prior == TRUE) {
        # load prior
        prior_value <- params[p, "X2"]
        prior_name <- params[p, "X4"]
        prior_args_string <- params[p, "X5"]
        prior_args <- list()
        aux <- unlist(strsplit(prior_args_string, dict_delim))
        for (k in seq(1, length(aux), 2)) {
          prior_args[[aux[k]]] <- as.numeric(aux[k + 1])
        }
        # compute prior pdf
        prior_x <- seq(min_val, max_val, binwidth)
        if (prior_name == "uniform") {
          prior_y <- dunif(
            prior_x,
            min = prior_args[["low"]],
            max = prior_args[["up"]]
          )
        } else if (prior_name == "normal") {
          prior_y <- dnorm(
            prior_x,
            mean = prior_args[["mean"]],
            sd = prior_args[["sd"]]
          )
        } else if (prior_name == "lognorm") {
          prior_y <- dlnorm(
            prior_x,
            meanlog = prior_args[["meanlog"]],
            sdlog = prior_args[["sdlog"]]
          )
        } else if (prior_name == "gamma") {
          prior_y <- dgamma(
            prior_x,
            shape = prior_args[["shape"]],
            scale = prior_args[["scale"]]
          )
        } else if (prior_name == "exp") {
          prior_y <- dexp(
            prior_x,
            rate = prior_args[["rate"]]
          )
        }
        prior_df <- tibble(xaxis = prior_x, yaxis = prior_y)
      }

      # 1 : plot Monte Carlo Markov Chain
      mcmc_plot <- ggplot(df, aes(x = iter, y = .data[[p]])) +
        geom_line(color = colors[i]) +
        geom_hline(yintercept = mu, color = clr_stat) +
        geom_hline(yintercept = mu - sigma, linetype = "dashed", color = clr_stat) +
        geom_hline(yintercept = mu + sigma, linetype = "dashed", color = clr_stat)
      if ((is_target == TRUE) && (show_prior == TRUE)) {
        mcmc_plot <- mcmc_plot +
          geom_hline(yintercept = prior_value, color = clr_target)
      }
      mcmc_plot <- mcmc_plot +
        ylim(min_val, max_val) +
        labs(x = ifelse(i == length(selec_cols), "iterations", ""), y = p)

      # 2 : plot histogram of values
      hist_plot <- ggplot()
      if (p %in% rownames(params)) {
        hist_plot <- hist_plot +
          geom_line(data = prior_df, mapping = aes(x = xaxis, y = yaxis), color = clr_prior)
      }
      hist_plot <- hist_plot +
        geom_freqpoly(data = df, mapping = aes(x = .data[[p]], y = after_stat(density)), binwidth = binwidth, color = colors[i]) +
        geom_vline(xintercept = mu, color = clr_stat) +
        geom_vline(xintercept = mu - sigma, linetype = "dashed", color = clr_stat) +
        geom_vline(xintercept = mu + sigma, linetype = "dashed", color = clr_stat)
      if ((is_target == TRUE) && (show_prior == TRUE)) {
        hist_plot <- hist_plot +
          geom_vline(xintercept = prior_value, color = clr_target)
      }
      hist_plot <- hist_plot +
        xlim(min_val, max_val) +
        labs(x = p, y = ifelse(i == length(selec_cols), "density", "")) +
        coord_flip()

      # 3 : plot corrleation function
      x <- thinning * seq(0, acorr_n_period * period, 1)
      acorr_df <- data.frame(iter = x, acorr = acorr[seq(1, length(x), 1)])
      corr_plot <- ggplot(acorr_df, aes(x = iter, y = acorr)) +
        geom_hline(yintercept = 0, color = clr_stat) +
        geom_line(color = colors[i]) +
        geom_vline(xintercept = period * thinning, linetype = "dashed", color = clr_stat)
      corr_plot <- corr_plot +
        labs(x = ifelse(i == length(selec_cols), "iterations", ""), y = "auto-corr")

      # append plots to the plot list
      plot_list[[p]] <- ggplotGrob(mcmc_plot)
      plot_list[[paste0(p, "-hist")]] <- ggplotGrob(hist_plot)
      plot_list[[paste0(p, "-corr")]] <- ggplotGrob(corr_plot)
    }
  }

  # save computations
  write_delim(
    output_df,
    output_file,
    delim = delim
  )

  # save figure
  if (save_fig == TRUE) {
    arranged_plots <- plot_grid(
      plotlist = plot_list,
      ncol = 3,
      nrow = length(selec_cols),
      rel_widths = c(2, 0.75, 1.25),
      align = "v"
    )
    ggsave(
      figure_file,
      plot = arranged_plots,
      device = "pdf",
      dpi = 300,
      width = 11,
      height = 1.5 * length(selec_cols)
    )
  }
}

# EXECUTE
main()
