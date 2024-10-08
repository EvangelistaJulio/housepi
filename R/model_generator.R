#
# model_generator.R
#
# Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
#

#' Copy a default model to a given directory.
#' @param which default model name.
#' @param path path where to copy implementation and header.
#' @export
template_model <- function(which = NULL, path = NULL) {
  # system path to housepi
  models_dir <- file.path(system.file(package = "housepi"), "models")
  
  if (is.null(which)) {
    cli::cli_alert_info("Default model names:")
    ulid <- cli::cli_ul()
    for (name in list.files(models_dir)) {
      cli::cli_li(name)
    }
    cli::cli_end(ulid)

  } else if (is.character(which)){
    if (which %in% list.files(models_dir)) {
      # create directory
      if (is.null(path)) {
        new_dir <- file.path(getwd(), paste0("user_", which))
      } else {
        new_dir <- path
      }

      dir.create(new_dir, recursive = TRUE, showWarnings = TRUE)
      if (!dir.exists(new_dir)) {
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
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.",
          i = "Check your user permissions.")))
      }
    }
  } else {
    stop(cli::format_error(
      c("Provided model is not a default model name.",
      i = "To see the available default models, run: template_model()")))
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Generate a model
#' @description
#' `generate_model` generates a model from a default or from a user script.
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
      c(paste0("`", model, "` is not an existing folder nor a default model."),
      i = "To see available default models, run: template_model()")))
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
    if (!dir.exists(new_dir)) {
      stop(cli::format_error(
        c("Impossibile to create a temporary directory for the model.",
        i = "Check your user permissions.")))
    }

    # copy files
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(housepi_dir, "template", f)
      to_file <- file.path(new_dir, f)
      if (!file.copy(from_file, to_file, overwrite = TRUE)) {
        stop(cli::format_error(
          c("Failed to setup internal files.",
          i = "Check your user permissions.")))
      }
    }

    # copy folders
    for (fldr in c("src", "R", "man")) {
      from_path <- file.path(housepi_dir, "template", fldr)
      if (!file.copy(from_path, new_dir, recursive = TRUE, overwrite = TRUE)) {
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
    if (!file.copy(cpp_file, file.path(new_dir, "src", "model.cpp"), overwrite = TRUE)) {
      stop(cli::format_error(c("Failed to setup implementation file.",
                               i = "Check your user permissions.")))
    }

    # copy header
    if (!file.copy(hpp_file, file.path(new_dir, "src", "model.h"), overwrite = TRUE)) {
      stop(cli::format_error(c("Failed to setup header file.",
                               i = "Check your user permissions.")))
    }
    
    # copy templated parameters
    templ_pars_file <- file.path(housepi_dir, "models", model, "templated_params")
    if (!file.copy(templ_pars_file, file.path(new_dir, "src", "templated_params"), overwrite = TRUE)) {
      cli::cli_alert_warning(cli::format_warning(c("Missing templated parameters.")))
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
#' @param setup list or file path with the setup definitions
#' @param params dataframe or file path with parameters definitions.
#' @param data dataframe or file path with input data.
#' @param augm dataframe or file with augmented variable definitions.
#' @param verbose boolean flag for automatic prints.
#' @export
build_model <- function(hemi,
                        setup = NULL,
                        params = NULL,
                        data = NULL,
                        augm = NULL,
                        verbose = FALSE) {
  store_scipen <- options()$scipen
  options(scipen = 500)

  #===== 1) SETUP =====#
  if (is.null(setup)) {
    hemi$read_setup("<NA>")

  } else if (is.character(setup)) {
    hemi$read_setup(setup)

  } else if ((is.list(setup)) && (!is.data.frame(setup))) {
    hemi$import_setup(lapply(setup, as.character))

  } else {
    stop(cli::format_error(c("Wrong type for setup.",
    i = "It should be list or character (path/to/file).")))
  }

  #===== 2) PARAMETERS =====#
  if (is.null(params)) {
    hemi$read_params("<NA>")

  } else if (is.character(params)) {
    hemi$read_params(params)

  } else if (is.data.frame(params)) {
    hemi$import_params(data.frame(lapply(params, as.character)))

  } else {
    stop(cli::format_error(c("Wrong type for params",
    i = "It should be datafram or character (path/to/file).")))
  }
  
  #===== 3) DATA =====#
  if (is.null(data)) {
    hemi$read_data("<NA>")

  } else if (is.character(data)) {
    hemi$read_data(data)

  } else if (is.data.frame(data)) {
    hemi$import_data(data.frame(lapply(data, as.character)))

  } else {
    stop(cli::format_error(c("Wrong type for data",
    i = "It should be datafram or character (path/to/file).")))
  }

  #===== 4) AUGMENTED VARIABLES =====#
  if (is.null(augm)) {
    hemi$read_augm("<NA>")

  } else if (is.character(augm)) {
    hemi$read_augm(augm)

  } else if (is.data.frame(augm)) {
    hemi$import_augm(data.frame(lapply(augm, as.character)))

  } else {
    stop(cli::format_error(c("Wrong type for augm",
    i = "It should be datafram or character (path/to/file).")))
  }

  options(scipen = store_scipen)
}

#' Run the inference.
#' @param hemi household model object.
#' @param progbar boolean flag to suppress progress bar.
#' @returns Return a coda mcmc.
#' @export
run_inference <- function(hemi, progbar = TRUE) {

  hemi$infer(progbar)
  df <- hemi$traces()
  
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

#' Run the simulation.
#' @param hemi household model object.
#' @param progbar boolean flag to suppress progress bar.
#' @returns Returns a data.frame.
#' @export
run_simulation <- function(hemi, progbar = TRUE) {

  hemi$simulate(progbar)

  return(hemi$data())
}