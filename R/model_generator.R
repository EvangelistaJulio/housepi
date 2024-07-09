#' Generate a model from a cpp script.
#' 
#' @param modelpath folder with implementation model.cpp and header model.h.
#' @param quiet quiet the compilation messages (default).
#' @param force generate a model even if a cached model already exist
#' 
#' @returns Return the user model object.
#' @export
generate <- function(modelpath, quiet=TRUE, force=FALSE) {
  # system path to hemirt
  hemirt_dir <- system.file(package = "hemirt")

  # ----- Check that files exist ----- #

  if(dir.exists(modelpath)) {

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
    new_name <- sprintf("hemirt%x", as.integer(file.info(cpp_file)$mtime))
    new_dir <- file.path(tempdir(), basename(modelpath), new_name)
  }
  else if(modelpath %in% list.files(file.path(hemirt_dir, "models"))){

    # implementation file
    cpp_file <- file.path(hemirt_dir, "models", modelpath, "model.cpp")

    # header file
    hpp_file <- file.path(hemirt_dir, "models", modelpath, "model.h")

    # temporary folder
    new_name <- paste0("hemirt", gsub("_", "", modelpath))
    new_dir <- file.path(tempdir(), new_name)

  }
  else {
    stop(paste0("Unrecognized model path: ", modelpath))
  }

  # ----- Create model package ----- #

  # executable name
  executable <- file.path(new_dir, "src", paste0(new_name, ".so"))

  # option to force the removal of the package directory
  if ((force == TRUE) | (file.exists(executable))) {
    unlink(new_dir, recursive = TRUE)
  }

  if (file.exists(executable)) {
    # use cached version
    warning(paste0("Using cached version of ", new_name)) 

  } else {
    # copy the package to the new directory
    dir.create(new_dir, recursive=TRUE, showWarnings = FALSE)
    src_files <- c("NAMESPACE", "DESCRIPTION", "src", "R")
    
    # copy files
    for (f in c("NAMESPACE", "DESCRIPTION")) {
      from_file <- file.path(hemirt_dir, "template", f)
      to_file <- file.path(new_dir, f)
      file.copy(from_file, to_file, overwrite = TRUE)
    }
    
    # copy folders
    for (fldr in c("src", "R")) {
      from_path <- file.path(hemirt_dir, "template", fldr)
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
        
    # copy implementation
    file.copy(cpp_file, file.path(new_dir, "src/model.cpp"), 
    recursive = FALSE, overwrite = TRUE)

    # copy header
    file.copy(hpp_file, file.path(new_dir, "src/model.h"),
    recursive = FALSE, overwrite = TRUE)
          
    # ----- compile package ----- #
    # FIXME : how to enhance pre-compilation?
    # FIXME : annoying Warning related to NAMESPACE
    pkgbuild::compile_dll(new_dir, 
                          quiet = quiet,
                          debug = FALSE)
  }
  
  #devtools::document(new_dir)
  devtools::load_all(new_dir)
  
  # extract the class definition and provide new object
  return(methods::new(get(paste0("Model_", new_name))))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#' Build user model to run inference or simulation
#'
#' @param x user model object
#' @param params file path for the model parameters
#' @param data file path for the input data
#' @param augm file path for augmented data
#' @param setup file path for the setup configurations
#' @param extra file path for any extra configurations
#' @param is_infer boolean flag to indicate whether to perform inference
#' @param verbose boolean flag to suppress automatic prints
#' @export
build <- function (
    x,
    params = "None",
    data = "None", 
    augm = "None", 
    setup = "None", 
    extra = "None",
    is_infer = TRUE,
    verbose = FALSE
) {
    
  x$build(params, data, augm,
          setup, extra, 
          is_infer, verbose)
}

#' Run the inference
#'
#' @param x user model object
#' @param verbose boolean flag to suppress automatic prints
#' @param progbar boolean flag to suppress progress bar
#' 
#' @returns Return a coda mcmc.
#' @export
infer <- function(x, verbose = FALSE, progbar = TRUE) {
  temp_output_file <- tempfile(fileext = "csv")
  x$inference(temp_output_file, verbose, progbar)
  mcmc_data <- utils::read.csv(temp_output_file, header = TRUE)
  unlink(temp_output_file)
  
  return(coda::as.mcmc(mcmc_data))
}

#' Run the simulation
#'
#' @param x user model object
#' 
#' @returns Return a data frame
#' @export
simulate <- function(x) {
  temp_output_file <- tempfile(fileext = "csv")
  x$simulation(temp_output_file, FALSE)
  data <- utils::read.csv(temp_output_file, header = TRUE)
  unlink(temp_output_file)
  
  return(data)
}

#' Compute the model loglikelihood
#'
#' @param x user model object
#' @param verbose boolean flag to suppress automatic prints
#' 
#' @returns Return a double
#' @export
compute_loglike <- function(x, verbose = FALSE) {
  return(x$loglikelihood(verbose))
}