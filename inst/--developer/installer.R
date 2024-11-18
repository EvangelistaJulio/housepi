mtime <- function(x) {
  x <- x[file.exists(x)]
  if (length(x) == 0) {
    return(NULL)
  }
  max(file.info(x)$mtime)
}

# List all source files in the package
sources <- function(path = ".") {
  srcdir <- file.path(path, "src")
  src <- c(
    dir(srcdir, "\\.(c.*|f|rs)$", recursive = TRUE, full.names = TRUE),
    dir(srcdir, "^Cargo\\.toml$", recursive = TRUE, full.names = TRUE)
  )
  extra <- desc::desc_get("Config/build/extra-sources", path)
  
  if (!is.na(extra)) {
    glb <- trimws(strsplit(extra, ",", fixed = TRUE)[[1]])
    xs <- file.path(path, globs(path, glb))
    xfls <- unlist(lapply(
      xs,
      dir,
      recursive = TRUE,
      full.names = TRUE
    ))
    src <- c(src, xs, xfls)
  }
  src
}

# List all header files in the package
headers <- function(path = ".") {
  incldir <- file.path(path, "inst", "include")
  srcdir <- file.path(path, "src")
  
  c(
    dir(srcdir, "^Makevars.*$", recursive = TRUE, full.names = TRUE),
    dir(srcdir, "\\.h.*$", recursive = TRUE, full.names = TRUE),
    dir(incldir, "\\.h.*$", recursive = TRUE, full.names = TRUE)
  )
}

dll_path <- function(path = ".") {
  name <- paste0(path, .Platform$dynlib.ext)
  file.path(path, "src", name)
}

needs_compile <- function(path = ".") {
  source <- mtime(c(sources(path), headers(path)))
  # no source files, so doesn't need compile
  if (is.null(source)) {
    return(FALSE)
  }
  
  dll <- mtime(dll_path(path))
  # no dll, so needs compile
  if (is.null(dll)) {
    return(TRUE)
  }
  
  source > dll
}

#=======================#
#===== COMPILE_DLL =====#
#=======================#

path = "tester"
force = FALSE
quiet = FALSE
debug = TRUE

# Determine if a full compile is needed
if (!needs_compile(path) && !isTRUE(force)) {
  message("No changes detected; skipping recompilation.")
  return(invisible())
}

pkgbuild::check_build_tools(quiet = TRUE)

# Rebuild RcppExports if needed
unlink(file.path(path, c("R/RcppExports.R", "src/RcppExports.cpp")))
Rcpp::compileAttributes(path, verbose = !quiet)

# Identify all source files and corresponding object files
src_path <- file.path(path, "src")
cpp_files <- list.files(src_path, pattern = "\\.cpp$", full.names = TRUE)
o_files <- gsub("\\.cpp$", ".o", cpp_files)

# Determine which files have changed
mask <- file.exists(o_files)
modified_files <- cpp_files[!mask]
if (any(mask))
{
  modified_files <- c(
    modified_files,
    cpp_files[mtime(cpp_files[mask]) > mtime(o_files[mask])]
  )
}

# Compile only the modified source files
if (length(modified_files) > 0) {
  for (file in modified_files) {
    callr::rcmd_safe("SHLIB", c(file), env = c(PKG_LIBS = "-shared"))
  }
} else {
  message("No modified files detected; skipping source recompilation.")
}

# Generate the DLL file by linking all object files
dll_name <- paste0(pkgbuild::pkg_name(path), .Platform$dynlib.ext)
dll_target <- file.path(src_path, dll_name)

o_files <- list.files(src_path, pattern = "\\.o$", full.names = TRUE)

if (length(o_files) == 0) {
  stop("No object files found to link into a shared library.")
}

message("Linking object files into DLL: ", dll_name)
callr::rcmd_safe(
  "SHLIB",
  c(paste(o_files, collapse = " "), "-o", dll_target),
  env = c(PKG_LIBS = "-shared")
)

message("DLL successfully generated: ", dll_target)
