## Custom script to build website
##

## locate template Rd files
##
model_dirs <- list.dirs("inst/models", recursive = FALSE)
rd_files <- file.path(model_dirs, "template_info.Rd")
# << assert
stopifnot(all(file.exists(rd_files)))
stopifnot(length(rd_files) == length(model_dirs))

## copy template Rd files into man/
##
copied_files <- file.path(
  "man", paste0("_hhmodel_", basename(dirname(rd_files)), ".Rd")
)
# << assert
stopifnot(all(file.copy(rd_files, copied_files, overwrite = TRUE)))
# add \keyword{internal} to copied files
for (f in copied_files) {
  lines <- readLines(f)
  if (!any(grepl("^\\\\keyword\\{internal\\}", lines))) {
    lines <- c(lines, "\\keyword{internal}")
    writeLines(lines, f)
  }
}

## build pkgdown site
##
pkgdown::build_site(install = FALSE)

## cleanup
##
file.remove(copied_files)
