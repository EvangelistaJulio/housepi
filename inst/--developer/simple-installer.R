#----------------------------#
#----- 1st installation -----#
#----------------------------#

pkgbuild::compile_dll("tester", debug = FALSE)
#dyn.load("tester/src/dummyjukebox.so")
devtools::load_all("tester", reset = TRUE, 
                   recompile = FALSE, export_all = FALSE)
x <- JukeBox$new()
x$pezzali()
x$collins()

# note: it compiles only changed cpp files!

#-------------------------------#
#----- custom installation -----#
#-------------------------------#

unlink("product", recursive = TRUE)
dir.create("product", recursive = TRUE)

# copy single package files
for (f in c("NAMESPACE", "DESCRIPTION")) {
  from_file <- file.path("tester", f)
  to_file <- file.path("product", f)
  if (!file.copy(from_file, to_file, overwrite = TRUE)) {
    stop(paste("Failed:", f))
  }
}

# copy folders from template
for (fldr in c("src", "R")) {
  from_path <- file.path("tester", fldr)
  if (!file.copy(from_path,
                 "product",
                 recursive = TRUE,
                 overwrite = TRUE)) {
    stop(paste("Failed:", fldr))
  }
}

pkgbuild::compile_dll("product", debug = FALSE)

#dyn.load("tester/src/dummyjukebox.so")
devtools::load_all("tester", reset = TRUE, 
                   recompile = FALSE, export_all = FALSE)


