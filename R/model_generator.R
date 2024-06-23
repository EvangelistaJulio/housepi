#' Generate a model from a cpp script.
#' 
#' @param model file containing the implementation file .cpp. The same folder must contain the header file .h.
#' @param quiet quiet the compilation messages (default).
#' @param force generate a model even if a cached model already exist
#' 
#' @returns Return an user model object.
generate_model <- function(model, quiet=TRUE, force=FALSE) {
  
  # temporary folder
  new_name <- sprintf("%x", as.integer(file.info(model)$mtime))
  new_name <- paste0("hemirt", new_name)
  
  new_dir <- file.path(tempdir(), basename(dirname(model)), new_name)
  
  executable <- file.path(new_dir, "src", paste0(new_name, ".so"))
  if ((force == TRUE) | (!file.exists(executable))) {
    dir.create(new_dir, recursive=TRUE, showWarnings = FALSE)
    
    # copy the package to the new directory
    src_dir <- file.path(system.file(package = "hemirt"), "template")
    src_files <- c("NAMESPACE", "DESCRIPTION", "src", "R")
    
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(src_dir, f)
      to_file <- file.path(new_dir, f)
      file.copy(from_file, to_file, overwrite = TRUE)
    }
    
    for (fldr in c("src", "R")) {
      from_path <- file.path(src_dir, fldr)
      file.copy(from_path, new_dir, recursive=TRUE, overwrite = TRUE)
    }
    
    # update with new name 
    src_files <- c("NAMESPACE", "DESCRIPTION", 
                   "src/export_bindings.cpp", "R/zzz.R") 
    for (f in src_files) {
      content <- readLines(file.path(new_dir, f), warn = FALSE)
      modified <- gsub("hemirt12345", new_name, content)
      writeLines(modified, file.path(new_dir, f))
    }
    
    # copy user header and implementation
    noext <- tools::file_path_sans_ext(model)
    
    # copy header
    if (file.exists(paste0(noext, ".h"))) {
      file.copy(
        paste0(noext, ".h"), file.path(new_dir, "src/model.h"),
        recursive = FALSE, overwrite = TRUE)
    } else if (file.exists(paste0(noext, ".hpp"))) {
      file.copy(
        paste0(noext, ".hpp"), file.path(new_dir, "src/model.h"),
        recursive = FALSE, overwrite = TRUE)
    } else {
      stop("A file ", paste0(noext, ".h"), " (or .hpp) does not exist.")
    }
      
    # copy implementation
    c_file <- paste0(noext, ".cpp")
    if (file.exists(c_file)) {
      file.copy(
        c_file, file.path(new_dir, "src/model.cpp"),
        recursive = FALSE, overwrite = TRUE)
    } else {
      stop("A file ", cpp_file, " does not exist.")
    }
    
    # FIXME : how to sue pre-compilation?
    # compile 
    pkgbuild::compile_dll(new_dir, 
                          quiet = quiet,
                          debug = FALSE)
  } else {
    warning(paste0("Using cached version of ", new_name))
  }
  
  #devtools::document(new_dir)
  devtools::load_all(new_dir)
  
  # extract the class definition and provide new object
  return(new(get(paste0("Model_", new_name))))
}

#' Build user model to run inference or simulation
#'
#' @param x user model object
#' @param params file path for the model parameters
#' @param data file path for the input data
#' @param augm file path for augmented data
#' @param setup file path for the setup configurations
#' @param extra file path for any extra configurations
#' @param is_infer boolean flag to indicate whether to perform inference
build_model <- function (
    x,
    params = "None",
    data = "None", 
    augm = "None", 
    setup = "None", 
    extra = "None",
    is_infer = TRUE,
    verbose = TRUE
) {
  
  if (params == "None") {
    stop("Paramter file can't be undefined.")
  }
  
  if (data == "None") {
    stop("Paramter file can't be undefined.")
  }
  
  x$build(is_infer, 
          params, data, augm,
          setup, extra, verbose)
}

#' Run the inference
#'
#' @param x user model object
#' 
#' @returns Return a coda mcmc.
run_hemirt_infer <- function(x) {
  temp_output_file <- tempfile(fileext = "csv")
  x$inference(temp_output_file)
  mcmc_data <- read.csv(temp_output_file, header = TRUE)
  unlink(temp_output_file)
  
  return(coda::as.mcmc(mcmc_data))
}


#' Run the simulation
#'
#' @param x user model object
#' 
#' @returns Return a data frame
run_hemirt_simulate <- function(x) {
  temp_output_file <- tempfile(fileext = "csv")
  x$simulation(temp_output_file)
  data <- read.csv(temp_output_file, header = TRUE)
  unlink(temp_output_file)
  
  return(data)
}