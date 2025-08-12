# locate template Rd files
model_dirs <- list.dirs("inst/models", recursive = FALSE)
rd_files <- file.path(model_dirs, "template_info.Rd")
rd_files <- rd_files[file.exists(rd_files)]

stopifnot(length(rd_files) > 0)

# copy Rd files into man/
copied_files <- file.path(
  "man", paste0("_template_", basename(dirname(rd_files)), ".Rd")
)
file.copy(rd_files, copied_files, overwrite = TRUE)

# parse Rd titles for the index
template_info <- lapply(seq_along(rd_files), function(i) {
  rd <- tools::parse_Rd(rd_files[i])
  title_idx <- which(sapply(rd, attr, "Rd_tag") == "\\title")
  title_txt <- if (length(title_idx)) paste(unlist(rd[[title_idx]]), collapse = "") else ""
  list(
    folder = basename(dirname(rd_files[i])),
    topic  = paste0("_template_", basename(dirname(rd_files))),
    title  = if (nzchar(title_txt)) title_txt else basename(dirname(rd_files[i]))
  )
})

# write templates/index.md
lines <- c(
  "---",
  'title: "Templates"',
  "---",
  "",
  "# Templates",
  "",
  "The following templates are available:"
)

for (info in template_info) {
  lines <- c(
    lines,
    sprintf("- **%s** — [%s](../reference/%s.html)",
            info$folder, info$title, info$topic)
  )
}

dir.create("pkgdown/templates", recursive = TRUE, showWarnings = FALSE)
writeLines(lines, "pkgdown/templates/index.md")

# build pkgdown site
pkgdown::build_site(install = F)

# cleanup
file.remove(copied_files)

