setwd("inst/template/")
withr::with_envvar(
  housepi:::get_cpp_environment_flags(),
  devtools::document()
)
setwd("../../")
