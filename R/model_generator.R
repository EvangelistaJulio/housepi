
#################################
#######  Global Variables #######
#################################

UNDEFINED <- "None"

###################################
####### Auxiliary Functions #######
###################################

# Get the system directory of a model.
get_usr_tmp_file <- function(hemi, name) {
  usr_dir <- system.file(package = attr(class(hemi), "package"))
  output <- file.path(usr_dir, "tmp", paste0(name, ".csv"))
  dir.create(dirname(output),
             recursive = TRUE,
             showWarnings = FALSE)
  return(output)
}

# Create temporary csv file if input is a table (wrapper).
# FIXME : this should be improved allowing dataframe I/O within Rcpp
cpp_csv_commuter <- function(hemi, input, name = UNDEFINED) {
  if (is.character(input)) {
    output <- input
  } else if (is.data.frame(input) || is.matrix(input) || is.list(input)) {
    output <- get_usr_tmp_file(hemi, name)
    write.csv(input,
              file = output,
              row.names = FALSE,
              quote = FALSE)
  } else {
    stop(cli::format_error(paste0("Unrecognized input for ", name)))
  }
  return(output)
}

# Default model list message.
show_default_names <- function() {
  models_dir <- file.path(system.file(package = "housepi"), "models")
  cli::cli_alert_info("Default model names:")
  ulid <- cli::cli_ul()
  for (name in list.files(models_dir)) {
    cli::cli_li(name)
  }
  cli::cli_end(ulid)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Copy a default model to a given directory.
#' @param which default model name.
#' @param path path where to copy implementation and header.
#' @export
template_model <- function(which = UNDEFINED, path = UNDEFINED) {
  # system path to housepi
  models_dir <- file.path(system.file(package = "housepi"), "models")

  # check if default exists
  if (!is.character(which)) {
    stop(cli::format_error(
      c("Provided model is not a default model name.",
      i = "To see the available default models, run: template_model()")))
    
  } else if (which == UNDEFINED) {
    show_default_names()

  } else {
    if (which %in% list.files(models_dir)) {
      # create directory
      if (path == UNDEFINED) {
        new_dir <- file.path(getwd(), paste0("user_", which))
      } else {
        new_dir <- path
      }

      dir.create(new_dir, recursive = TRUE, showWarnings = TRUE)
      if(!dir.exists(new_dir)) {
        stop(cli::format_error(
          c("Impossibile to create a temporary directory for the model.",
          i = "Check your user permissions.")))
      }
      
    } else {
      stop(cli::format_error(
        c("Provided model is not a default model name.",
        i = "To see the available default models, run: template_model()")))
    }

    # copy files
    for (f in c("model.cpp", "model.h")) {
      from_file <- file.path(models_dir, which, f)
      to_file <- file.path(new_dir, f)
      if(!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.",
          i = "Check your user permissions.")))
      }
    }
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Generate a model
#' 
#' @description
#' `generate_model` generates a model from a default or from a user script.
#' 
#' @param model a string.\itemize{
#' \item default model name.
#' \item path/to/folder with implementation model.cpp and header model.h.
#' }
#' @param quiet hide compilation commands.
#' @param force generate a model even if a cached model already exist.
#' @returns Return the household model object.
#' @export
generate_model <- function(model, quiet = TRUE, force = FALSE) {

  # system path to housepi
  housepi_dir <- system.file(package = "housepi")
  
  # temporary per-session directory where to store model environment
  tmp_env_path <- tempdir()

  if (dir.exists(model)) {
    # ----- Custom user model ----- #
    # implementation file
    cpp_file <- file.path(model, "model.cpp")
    if (!file.exists(cpp_file)) {
      stop(cli::format_error(paste0("A file `", cpp_file, "` does not exist.")))
    }

    # header file
    hpp_file <- file.path(model, "model.h")
    if (!file.exists(hpp_file)) {
      stop(cli::format_error(paste0("A file `", hpp_file, "` does not exist.")))
    }

    # temporary folder
    new_name <- sprintf("housepi%x", as.integer(file.info(cpp_file)$mtime))
    new_dir <- file.path(tmp_env_path, basename(model), new_name)

  } else if (model %in% list.files(file.path(housepi_dir, "models"))) {
    # ----- Default model ----- #
    # implementation file
    cpp_file <- file.path(housepi_dir, "models", model, "model.cpp")

    # header file
    hpp_file <- file.path(housepi_dir, "models", model, "model.h")

    # temporary folder
    new_name <- paste0("housepi", gsub("_", "", model))
    new_dir <- file.path(tmp_env_path, new_name)

  } else {
    stop(cli::format_error(
      c(paste0("Unrecognized model: `", model, 
      "` is nor a default name nor a path to existing folder."),
      i = "To see the available default models, run: template_model()")))
  }

  # ----- Create model package ----- #

  # Returns the path/to/DLL file supposedly generated in this session
  executable <- file.path(new_dir, "src", paste0(new_name, .Platform$dynlib.ext))

  # option to force the removal of the package directory
  if ((force == TRUE) && (file.exists(executable))) {
    unlink(new_dir, recursive = TRUE)
  }

  # deal with existing version
  if (file.exists(executable)) {
    # use cached version
    cli::cli_alert_info("Using cached {.pkg {pkgload::pkg_name(new_dir)}}")

  } else {
    # copy the package to the new directory
    dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)
    if(!dir.exists(new_dir)) {
      stop(cli::format_error(
        c("Impossibile to create a temporary directory for the model.",
        i = "Check your user permissions.")))
    }

    # copy files
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(housepi_dir, "template", f)
      to_file <- file.path(new_dir, f)
      if(!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.",
          i = "Check your user permissions.")))
      }
    }

    # copy folders
    for (fldr in c("src", "R", "man")) {
      from_path <- file.path(housepi_dir, "template", fldr)
      if(!file.copy(from_path, new_dir, recursive = TRUE, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.",
          i = "Check your user permissions.")))
      }
    }

    # update with new name
    src_files <- c("man/housepi12345-package.Rd",
                   "man/Modelhousepi12345.Rd",
                   "src/export_bindings.cpp",
                   "R/zzz.R",
                   "NAMESPACE",
                   "DESCRIPTION")
    for (f in src_files) {
      fpath <- file.path(new_dir, f)
      content <- readLines(fpath, warn = FALSE)
      modified <- gsub("housepi12345", new_name, content)
      writeLines(modified, fpath)
      file.rename(fpath, gsub("housepi12345", new_name, fpath))
    }

    # copy implementation
    if(!file.copy(cpp_file, file.path(new_dir, "src", "model.cpp"), overwrite = TRUE)) {
      stop(cli::format_error(c("Failed to setup implementation file.",
                               i = "Check your user permissions.")))
    }

    # copy header
    if(!file.copy(hpp_file, file.path(new_dir, "src", "model.h"), overwrite = TRUE)) {
      stop(cli::format_error(c("Failed to setup header file.",
                               i = "Check your user permissions.")))
    }

    # ----- compile package ----- #
    # FIXME : how to enhance pre-compilation?
    pkgbuild::compile_dll(new_dir, quiet = quiet, debug = FALSE)
  }

  # load locally the package
  devtools::load_all(new_dir, recompile = FALSE, quiet = FALSE)

  # extract the class definition and provide new object
  return(methods::new(get(paste0("Model", new_name))))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Build model to run inference or simulation.
#' @param hemi household model object.
#' @param setup file path for the setup configurations.
#' @param params file path for the model parameters.
#' @param data file path for the input data.
#' @param augm file path for augmented data.
#' @param extra file path for any extra configurations.
#' @param verbose boolean flag to suppress automatic prints.
#' @export
build_model <- function(hemi,
                        setup = UNDEFINED,
                        params = UNDEFINED,
                        data = UNDEFINED,
                        augm = UNDEFINED,
                        extra = UNDEFINED,
                        verbose = FALSE) {
  hemi$build(cpp_csv_commuter(hemi, setup, name = "setup"), 
             cpp_csv_commuter(hemi, params, name = "params"),
             cpp_csv_commuter(hemi, data, name = "data"), 
             cpp_csv_commuter(hemi, augm, name = "augm"), 
             cpp_csv_commuter(hemi, extra, name = "extra"), 
             verbose)
}

#' Run the inference.
#' @param hemi household model object.
#' @param progbar boolean flag to suppress progress bar.
#' @returns Return a coda mcmc.
#' @export
run_inference <- function(hemi,
                          progbar = TRUE) {
  tmp_out_file <- get_usr_tmp_file(hemi, "infer")
  hemi$infer(tmp_out_file, progbar)
  df <- utils::read.csv(tmp_out_file, header = TRUE)
  unlink(tmp_out_file)
  
  # extract setup information for coda
  setup <- hemi$setup()
  
  # Check if multiple chains
  if ("CHAIN__" %in% colnames(df)) {
    df_list <- lapply(split(df, df[["CHAIN__"]]), function(x) {
      coda::mcmc(data = x[, !names(x) %in% "CHAIN__"],
                 start = setup[["burn_in"]] + setup[["thinning"]],
                 end = setup[["num_iter"]],
                 thin = setup[["thinning"]])
    })
    output <- coda::as.mcmc.list(df_list)
    
  } else {
    output <- coda::mcmc(data = df,
                         start = setup[["burn_in"]] + setup[["thinning"]],
                         end = setup[["num_iter"]],
                         thin = setup[["thinning"]])
  }
  
  return(output)
}