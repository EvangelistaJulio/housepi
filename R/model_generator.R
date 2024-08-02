# ----- Auxiliary Functions ----- #

# get the system directory of a model
get_usr_tmp_file <- function(x, name) {
  usr_dir <- system.file(package = attr(class(x), "package"))
  output <- file.path(usr_dir, "tmp", paste0(name, ".csv"))
  dir.create(dirname(output),
             recursive = TRUE,
             showWarnings = FALSE)
  
  return(output)
}

# create temporary csv file if input is a table (wrapper)
# note : can this be improved allowing dataframe I/O within Rcpp ?
cpp_csv_commuter <- function(x, input, name = "None") {
  if (is.character(input)) {
    output <- input
    
  } else if (is.data.frame(input) ||
             is.matrix(input) || is.list(input)) {
    output <- get_usr_tmp_file(x, name)
    write.csv(input,
              file = output,
              row.names = FALSE,
              quote = FALSE)
    
  } else {
    stop(paste0("Unrecognized input for ", name))
  }
  
  return(output)
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Copy a default model to a given directory.
#' @param which default model name.
#' @param path path where to copy implementation and header.
#' @export
template_model <- function(which = "None", path = "None") {
  # system path to housepi
  models_dir <- file.path(system.file(package = "housepi"), "models")
  
  # check if default exists
  if ((!is.character(path)) || (which == "None")) {
    defaults <- "Default models:"
    for (name in list.files(models_dir)) {
      defaults <- paste0(defaults, "\n* ", name)
    }
    
    if (!is.character(path)) {
      oss <- paste0("Unknown default model `", which, "`.\n", defaults)
    } else {
      stop(oss)
    }
    
  } else {
    if (which %in% list.files(models_dir)) {
      # create directory
      if (path == "None") {
        new_dir <- file.path(getwd(), paste0("user_", which))
        
      } else {
        new_dir <- path
      }
      
      dir.create(new_dir, recursive = TRUE, showWarnings = TRUE)
      
    } else {
      stop("The provided path is not valuable.")
    }
    
    # copy files
    for (f in c("model.cpp", "model.h")) {
      from_file <- file.path(models_dir, which, f)
      to_file <- file.path(new_dir, f)
      file.copy(from_file, to_file, overwrite = TRUE)
    }
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Generate a model from a cpp script.
#' @param modelpath folder with implementation model.cpp and header model.h.
#' @param quiet quiet the compilation messages (default).
#' @param force generate a model even if a cached model already exist.
#' @returns Return the household model object.
#' @export
generate_model <- function(modelpath, 
                           quiet = TRUE, 
                           force = FALSE) {
  # system path to housepi
  housepi_dir <- system.file(package = "housepi")
  
  # ----- Check that files exist ----- #
  
  if (dir.exists(modelpath)) {
    # implementation file
    cpp_file <- file.path(modelpath, "model.cpp")
    if (!file.exists(cpp_file)) {
      stop("A file ", cpp_file, " does not exist.")
    }
    
    # header file
    hpp_file <- file.path(modelpath, "model.h")
    if (!file.exists(hpp_file)) {
      stop("A file ", hpp_file, " does not exist.")
    }
    
    # temporary folder
    new_name <- sprintf("housepi%x", as.integer(file.info(cpp_file)$mtime))
    new_dir <- file.path(tempdir(), basename(modelpath), new_name)
    
  } else if (modelpath %in% list.files(file.path(housepi_dir, "models"))) {
    # implementation file
    cpp_file <- file.path(housepi_dir, "models", modelpath, "model.cpp")
    
    # header file
    hpp_file <- file.path(housepi_dir, "models", modelpath, "model.h")
    
    # temporary folder
    new_name <- paste0("housepi", gsub("_", "", modelpath))
    new_dir <- file.path(tempdir(), new_name)
    
  } else {
    stop(paste0("Unrecognized model path: ", modelpath))
  }
  
  # ----- Create model package ----- #
  
  # executable name
  executable <- file.path(new_dir, "src", paste0(new_name, ".so"))
  
  # option to force the removal of the package directory
  if ((force == TRUE) || (file.exists(executable))) {
    unlink(new_dir, recursive = TRUE)
  }
  
  if (file.exists(executable)) {
    # use cached version
    warning(paste0("Using cached version of ", new_name))
    
  } else {
    # copy the package to the new directory
    dir.create(new_dir, recursive = TRUE, showWarnings = FALSE)
    
    # copy files
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(housepi_dir, "template", f)
      to_file <- file.path(new_dir, f)
      file.copy(from_file, to_file, overwrite = TRUE)
    }
    
    # copy folders
    for (fldr in c("src", "R")) {
      from_path <- file.path(housepi_dir, "template", fldr)
      file.copy(from_path,
                new_dir,
                recursive = TRUE,
                overwrite = TRUE)
    }
    
    # update with new name
    src_files <- c("NAMESPACE",
                   "DESCRIPTION",
                   "src/export_bindings.cpp",
                   "R/zzz.R")
    for (f in src_files) {
      content <- readLines(file.path(new_dir, f), warn = FALSE)
      modified <- gsub("housepi12345", new_name, content)
      writeLines(modified, file.path(new_dir, f))
    }
    
    # copy implementation
    file.copy(cpp_file, file.path(new_dir, "src/model.cpp"), overwrite = TRUE)
    
    # copy header
    file.copy(hpp_file, file.path(new_dir, "src/model.h"), overwrite = TRUE)
    
    # ----- compile package ----- #
    # FIXME : how to enhance pre-compilation?
    # FIXME : annoying Warning related to NAMESPACE
    pkgbuild::compile_dll(new_dir, quiet = quiet, debug = FALSE)
  }
  
  # FIXME : should I run devtools::document(new_dir) ?
  devtools::load_all(new_dir)
  
  # extract the class definition and provide new object
  return(methods::new(get(paste0("Model", new_name))))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Build model to run inference or simulation.
#' @param hemi household model object.
#' @param params file path for the model parameters.
#' @param data file path for the input data.
#' @param augm file path for augmented data.
#' @param setup file path for the setup configurations.
#' @param extra file path for any extra configurations.
#' @param verbose boolean flag to suppress automatic prints.
#' @export
build_model <- function(hemi,
                        params = "None",
                        data = "None",
                        augm = "None",
                        setup = "None",
                        extra = "None",
                        verbose = FALSE) {
  hemi$build(cpp_csv_commuter(hemi, params, name = "params"),
             cpp_csv_commuter(hemi, data, name = "data"), 
             cpp_csv_commuter(hemi, augm, name = "augm"), 
             cpp_csv_commuter(hemi, setup, name = "setup"), 
             cpp_csv_commuter(hemi, extra, name = "extra"), 
             verbose)
}

#' Run the inference.
#' @param hemi household model object.
#' @param progbar boolean flag to suppress progress bar.
#' @param verbose boolean flag to suppress automatic prints.
#' @returns Return a coda mcmc.
#' @export
run_inference <- function(hemi,
                          progbar = TRUE,
                          verbose = FALSE) {
  tmp_out_file <- get_usr_tmp_file(hemi, "infer")
  hemi$infer(tmp_out_file, verbose, progbar)
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

#' Run the simulation.
#' @param hemi household model object.
#' @returns data frame.
#' @export
run_simulation <- function(hemi) {
  tmp_out_file <- get_usr_tmp_file(hemi, "simulate")
  hemi$simulate(tmp_out_file, FALSE)
  data <- utils::read.csv(tmp_out_file, header = TRUE)
  unlink(tmp_out_file)
  
  return(data)
}