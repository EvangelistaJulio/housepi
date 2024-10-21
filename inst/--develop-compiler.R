model <- "chain_binom"
quiet <- TRUE
debug <- FALSE

#----

# system path to housepi
housepi_dir <- system.file(package = "housepi")
template_dir <- file.path(housepi_dir, "template")

# FLAG
template_dir <- getwd()
# FLAG-END

#----- template model -----#
source_files <- file.path(housepi_dir,
                          "models",
                          model,
                          c("model.cpp", "model.h", "templated_params"))

# temporary folder
env_name <- paste0("housepi", gsub("_", "", model))
# FIXME :
# env_name <- sprintf("housepi%x", as.integer(file.info(cpp_file)$mtime))
# temporary session directory where to store model environment
exec_dir <- file.path(tempdir(), env_name)

# ----- Create model package ----- #

# Returns the path/to/DLL file supposedly generated in this session
src_dir <- file.path(exec_dir, "src")
executable <- file.path(src_dir,
                        paste0("housepienvironment", .Platform$dynlib.ext))

# FLAG
unlink(exec_dir, recursive = TRUE)
# FLAG-END

# copy the template folder
fs::dir_copy(template_dir, exec_dir)

# copy source files in the new directory
if (!all(file.copy(source_files, src_dir, overwrite = TRUE))) {
  stop(cli::format_error("Failed to setup model implementation files."))
}

# <=====
  
# Navigate to src directory
old_dir <- setwd(src_dir)
on.exit(setwd(old_dir), add = TRUE)

# load shared library
dyn.load(executable)

# expose internally the module
m <- Rcpp::Module("model_housepi_rcpp", PACKAGE = "housepienvironment")

#----- 

obj <- m$ModelHousepi$new()
