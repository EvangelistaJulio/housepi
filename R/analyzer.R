#!/usr/bin/env Rscript

# LIBRARIES
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(optparse))
options(warn=-1)

# ANALYSIS : options
PERC_BURN_IN <- 0.10 
MIN_MCMC_VAR <- 1e-10
UNIROOT_TOL <- 1e-6
UNIROOT_MAXITER <- 1e3
ACORR_N_PERIOD <- 5
NUM_BINS <- 22

# FIGURE : options
COL_STAT <- "black"
COL_PRIOR <- "grey"
# Define a custom theme and set it as the default theme for all plots
custom_theme <- theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = "black"),
    panel.border = element_rect(color = "black", fill = NA),
    #text = element_text(size=16,  family="Comic Sans MS")
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
      make_option(c("--save_fig"), action = "store_true", default = FALSE, help = "Save figure to default")
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

  # Load setup 
  if (is.null(setup_file)) {
    delim <- ","
    burn_in <- 0
  } else if (file.exists(setup_file)) {
    setup <- read_csv(setup_file, col_names = FALSE, show_col_types = FALSE)
    setup <- setup %>% column_to_rownames(., var = "X1")
    delim <- as.character(setup["delim", ])
    num_iter <- as.integer(setup["num_iter", ])
    burn_in <- floor(PERC_BURN_IN * num_iter)
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
  if (is.null(setup_file)) {
    # default
    out_path <- file.path(dirname(input_file), "analysis")
    dir.create(out_path, showWarnings = FALSE)
    output_file <- file.path(out_path, tools::file_path_sans_ext(basename(input_file)))
    output_file <- paste0(output_file, "-a.", tools::file_ext(input_file))
  } else if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file))
  }

  # Load figure 
  if (is.null(figure_file)) {
    if (save_fig == TRUE) {
      # default 
      pdf_path <- file.path(dirname(input_file), "figures")
      dir.create(pdf_path, showWarnings = FALSE)
      figure_file <- file.path(out_path, tools::file_path_sans_ext(basename(input_file)))
      figure_file <- paste0(figure_file, "-fig.pdf")
    }
  } else {
    save_fig = TRUE
    if (!dir.exists(dirname(figure_file))) {
      dir.create(dirname(figure_file))
    }
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

  # display a short summary
  cat("Performing analysis on file : ", input_file, "\n")
  cat("Results are saved on file : ", output_file, "\n")
  if (save_fig == TRUE) {
    cat("Figure is saved on file : ", figure_file, "\n")
  }

  # define working data frame and remove burn_in
  df <- input_df %>% filter(.data[["iter"]] > burn_in)
  sample_size <- length(df[["iter"]])
  # extract thinning
  thinning = df[["iter"]][2] - df[["iter"]][1]

  # define columns of interest 
  cols <- colnames(df)
  cols <- cols[cols != "iter"]
  cols <- cols[!grepl("__a", cols)]

  # initialize output
  output_df <- data.frame()
  if (save_fig == TRUE)
  {
    plot_list <- list()
    colors <- hue_pal()(length(cols))
  }

  # compute mean and dev.std for selected MCMC chains
  for (i in 1:length(cols)) {
    p <- cols[i]
    signal <- as.vector(df[[p]])

    if (var(signal) < MIN_MCMC_VAR) {
      # pathological case: signal is constant // UPDATE : find a way to not block the code
      cat("The variance of `", p, "` is lower than minimum MCMC variance ", MIN_MCMC_VAR, "\n")
      quit("no", 1)
    } 

    # UPDATE : apprently this can be more compact in tidyverse piping
    # compute autocorrelation function
    acorr <- acf(signal - mean(signal), lag.max=floor(length(signal) / 2), plot=FALSE)[["acf"]]
    x <- seq(0, length(acorr)-1, 1)
    # linear inerpolation
    interp_func <- approxfun(x, acorr, method="linear") 

    # find first root
    end <- x[which(acorr < 0)[1]]
    if (end < length(acorr)-1) {
      x_0 <- uniroot(interp_func, c(0, end), tol=UNIROOT_TOL, maxiter=UNIROOT_MAXITER, check.conv=TRUE)
      period <- as.integer(round(x_0[["root"]]))
    } else {
      # pathological case : acorr converges asymptotically to 0
      period <- length(acorr)-1
    }

    # mask for effective sample
    mask <- seq(1, sample_size, period)
    # effective sample size
    ess <- length(mask)
    # compute statistics on effective sample
    mu <- mean(signal[mask])
    sigma <- sd(signal[mask])

    key_data <- data.frame(key = p, mean = mu, sd = sigma, period = period * thinning, ESS = ess )
    output_df <- rbind(output_df, key_data)

    if (save_fig == TRUE)
    {
      # 1 : plot Monte Carlo Markov Chain
      min_y <- min(df[[p]])
      max_y <- max(df[[p]])
      mcmc_plot <- ggplot(df, aes(x = iter, y = .data[[p]])) + 
        geom_line(color = colors[i]) +
        geom_hline(yintercept = mu, color = COL_STAT) +
        geom_hline(yintercept = mu - sigma, linetype = "dashed", color = COL_STAT) +
        geom_hline(yintercept = mu + sigma, linetype = "dashed", color = COL_STAT) +
        ylim(min_y, max_y) + 
        labs(x = "iterations", y = p)
      plot_list[[p]] <- ggplotGrob(mcmc_plot)

      # 2 : plot histogram of values
      binwidth <- (max_y - min_y) / NUM_BINS
      hist_plot <- ggplot()
      if (p %in% rownames(params)) {
        # load prior
        prior_name <- params[p, "X4"]
        prior_args_string <- params[p, "X5"]
        prior_args <- list()
        aux <- unlist(strsplit(prior_args_string, "/"))
        for (k in seq(1, length(aux), 2)) {
          prior_args[[aux[k]]] <- as.numeric(aux[k+1])
        }
        # compute prior pdf
        X <- seq(min_y, max_y, binwidth)
        if (prior_name == "uniform") {
          Y <- dunif(X, min = prior_args[["low"]], max = prior_args[["up"]])
        } else if (prior_name == "normal") {
          Y <- dnorm(X, mean = prior_args[["mean"]], sd = prior_args[["sd"]])
        } else if (prior_name == "lognorm") {
          Y <- dlnorm(X, meanlog = prior_args[["mean"]], sdlog = prior_args[["sd"]])
        } else if (prior_name == "gamma") {
          Y <- dgamma(X, shape = prior_args[["shape"]], scale = prior_args[["scale"]])
        } else if (prior_name == "exp") {
          Y <- dexp(X, rate = prior_args[["rate"]])
        }

        hist_plot <- hist_plot +
          geom_line(data = tibble(xaxis=X, yaxis=Y), mapping = aes(x=xaxis, y=yaxis), color=COL_PRIOR)
      }
      hist_plot <- hist_plot  + 
        geom_freqpoly(data = df, mapping = aes(x = .data[[p]], y = after_stat(density)), binwidth = binwidth, color = colors[i]) +
        geom_vline(xintercept = mu, color = COL_STAT) +
        geom_vline(xintercept = mu - sigma, linetype = "dashed", color = COL_STAT) +
        geom_vline(xintercept = mu + sigma, linetype = "dashed", color = COL_STAT) +
        xlim(min_y, max_y) + 
        labs(x = p, y = "density") +
        coord_flip()
      plot_list[[paste0(p, "-hist")]] <- ggplotGrob(hist_plot)

      # 3 : plot corrleation function
      x <- thinning * seq(0, 5 * period,1 )
      acorr_df <- data.frame(iter = x, acorr = acorr[1:length(x)])
      corr_plot <- ggplot(acorr_df, aes(x = iter, y = acorr)) + 
        geom_hline(yintercept = 0, color = COL_STAT) +
        geom_line(color = colors[i]) +
        geom_vline(xintercept = period * thinning, linetype = "dashed", color = COL_STAT) + 
        labs(x = "iterations", y = "auto-corr")
      plot_list[[paste0(p, "-corr")]] <- ggplotGrob(corr_plot)
    }
  }
  
  # save computations
  write_delim(output_df, output_file, delim=delim )

  # save figure 
  if (save_fig == TRUE)
  {
    # FIXME : how do I avoid default opening of the figure?
    arranged_plots <- grid.arrange(grobs = plot_list, ncol = 3, widths = c(2, 0.75, 1.25))
    ggsave(figure_file, plot = arranged_plots, device = "pdf", dpi = 300, width = 8, height = 1.5 * length(cols))
    dev.off()
  }
}

# EXECUTE
main()
