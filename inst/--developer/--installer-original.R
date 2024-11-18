
compile_dll <- function(path = ".",
                        force = FALSE,
                        quiet = FALSE,
                        debug = TRUE) {
  path <- pkgbuild::pkg_path(path)
  
  if (!needs_compile(path) && !isTRUE(force)) {
    return(invisible())
  }
  
  pkgbuild::check_build_tools(quiet = TRUE)
  # FIXME : is this necessary?
  unlink(file.path(path, c("R/RcppExports.R", "src/RcppExports.cpp")))
  Rcpp::compileAttributes(path, verbose = !quiet)
  
  # Mock install the package to generate the DLL
  # FIXME : this requires to be recoded
  xflags <- pkgbuild::should_add_compiler_flags()
  
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
    
    invisible(dll_path(file.path(install_dir, pkgbuild::pkg_name(path))))
  }
  
  if (xflags) {
    pkgbuild::withr_with_makevars(compiler_flags(debug), {
      if (debug) {
        pkgbuild::withr_with_envvar(c(DEBUG = "true"), build())
      } else {
        build()
      }
    })
  } else {
    build()
  }
}

clean_dll <- function(path = ".") {
  path <- pkgbuild::pkg_path(path)
  
  # Clean out the /src/ directory and children:
  # * individual object files
  # * overall package definition file
  # * symbols.rds (added when run inside R CMD check)
  pattern <- sprintf(
    "\\.(o|sl|so|dylib|a|dll)$|(%s\\.def)$|^symbols.rds$",
    pkgbuild::pkg_name(path)
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
  name <- paste(pkgbuild::pkg_name(path), .Platform$dynlib.ext, sep = "")
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
  
  invisible(file.path(dest, pkgbuild::pkg_name(path)))
}

# Call R CMD 'command' with build tools active
#
# This is a wrapper around `callr::rcmd_safe()` that checks
# that you have build tools available, and on Windows, automatically sets
# the path to include Rtools.
#
# @param ... Parameters passed on to `rcmd_safe`.
# @param env Additional environment variables to set. The defaults from
#   `callr::rcmd_safe_env()` are always set.
# @inheritParams with_build_tools
# @inheritParams build
# @export
# @examples
# # These env vars are always set
# callr::rcmd_safe_env()
#
# if (has_build_tools()) {
#   rcmd_build_tools("CONFIG", "CC")$stdout
#   rcmd_build_tools("CC", "--version")$stdout
# }
rcmd_build_tools <- function(..., env = character(), required = TRUE, quiet = FALSE) {
  env <- c(callr::rcmd_safe_env(), env)

  if (!quiet) {
    cli::cat_rule(paste0("R CMD ", ..1), col = "cyan")
  }

  warn_for_potential_errors()

  callback <- if (cli::is_dynamic_tty()) {
    block_callback(quiet)
  } else {
    simple_callback(quiet)
  }

  res <- with_build_tools(
    {
      withCallingHandlers(
        callr::rcmd_safe(...,
          env = env, spinner = FALSE, show = FALSE,
          echo = FALSE, block_callback = callback, stderr = "2>&1"
        ),
        error = function(e) {
          if (!quiet) e$echo <- TRUE
          asNamespace("callr")$err$throw(e)
        }
      )
    },
    required = required
  )

  msg_for_long_paths(res)

  invisible(res)
}

msg_for_long_paths <- function(output) {
  if (is_windows() &&
    any(grepl("over-long path length", output$stdout))) {
    message(
      "\nIt seems that this package contains files with very long paths.\n",
      "This is not supported on most Windows versions. Please contact the\n",
      "package authors and tell them about this. See this GitHub issue\n",
      "for more details: https://github.com/r-lib/remotes/issues/84\n"
    )
  }
}

warn_for_potential_errors <- function() {
  if (is_windows() && grepl(" ", R.home()) &&
    getRversion() <= "3.4.2") {
    warning(
      immediate. = TRUE,
      "\n!!! Building will probably fail!\n",
      "This version of R has trouble with building packages if\n",
      "the R HOME directory (currently '", R.home(), "')\n",
      "has space characters. Possible workarounds include:\n",
      "- installing R to the C: drive,\n",
      "- installing it into a path without a space, or\n",
      "- creating a drive letter for R HOME via the `subst` windows command, and\n",
      "  starting R from the new drive.\n",
      "See also https://github.com/r-lib/remotes/issues/98\n"
    )
  }
}