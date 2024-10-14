compile_dll <- function(path = ".",
                        force = FALSE,
                        quiet = FALSE,
                        debug = TRUE) {
  path <- pkg_path(path)
  
  if (!needs_compile(path) && !isTRUE(force)) {
    return(invisible())
  }
  
  pkgbuild::check_build_tools(quiet = TRUE)
  # FIXME : is this necessary
  unlink(file.path(path, c("R/RcppExports.R", "src/RcppExports.cpp")))
  Rcpp::compileAttributes(path, verbose = !quiet)
  
  # Mock install the package to generate the DLL
  # FIXME : this requires to be recoded
  xflags <- should_add_compiler_flags()
  
  if (!quiet) {
    cli::cli_alert_info(c(
      "Compiling {.pkg {pkg_name(path)}}",
      if (debug && xflags) " (debug build)"
    ))
  }
  
  install_dir <- tempfile("devtools_install_")
  dir.create(install_dir)
  
  build <- function() {
    install_min(
      path,
      dest = install_dir,
      components = "libs",
      args = if (needs_clean(path)) "--preclean",
      quiet = quiet
    )
    
    invisible(dll_path(file.path(install_dir, pkg_name(path))))
  }
  
  if (xflags) {
    withr_with_makevars(compiler_flags(debug), {
      if (debug) {
        withr_with_envvar(c(DEBUG = "true"), build())
      } else {
        build()
      }
    })
  } else {
    build()
  }
}

clean_dll <- function(path = ".") {
  path <- pkg_path(path)
  
  # Clean out the /src/ directory and children:
  # * individual object files
  # * overall package definition file
  # * symbols.rds (added when run inside R CMD check)
  pattern <- sprintf(
    "\\.(o|sl|so|dylib|a|dll)$|(%s\\.def)$|^symbols.rds$",
    pkg_name(path)
  )
  files <- dir(
    file.path(path, "src"),
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE
  )
  unlink(files)
  
  invisible(files)
}

# Returns the full path and name of the DLL file
dll_path <- function(path = ".") {
  name <- paste(pkg_name(path), .Platform$dynlib.ext, sep = "")
  file.path(path, "src", name)
}

mtime <- function(x) {
  x <- x[file.exists(x)]
  if (length(x) == 0) {
    return(NULL)
  }
  max(file.info(x)$mtime)
}

globs <- function(path = ".", x) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(path)
  Sys.glob(x)
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

# Does the package need a clean compile?
# (i.e. is there a header or Makevars newer than the dll)
needs_clean <- function(path = ".") {
  headers <- mtime(headers(path))
  # no headers, so never needs clean compile
  if (is.null(headers)) {
    return(FALSE)
  }
  
  dll <- mtime(dll_path(path))
  # no dll, so needs compile
  if (is.null(dll)) {
    return(TRUE)
  }
  
  headers > dll
}

install_min <- function(path = ".", dest, components = NULL, args = NULL, quiet = FALSE) {
  stopifnot(is.character(dest), length(dest) == 1, file.exists(dest))
  
  poss <- c("R", "data", "help", "demo", "inst", "docs", "exec", "libs")
  if (!is.null(components)) {
    components <- match.arg(components, poss, several.ok = TRUE)
  }
  no <- setdiff(poss, components)
  no_args <- paste0("--no-", no)
  
  rcmd_build_tools(
    "INSTALL",
    c(
      path,
      paste("--library=", dest, sep = ""),
      no_args,
      "--no-multiarch",
      "--no-test-load",
      args
    ),
    fail_on_status = TRUE,
    quiet = quiet
  )
  
  invisible(file.path(dest, pkg_name(path)))
}