library(housepi)

official <- "/Users/fcamagli/Documents/SANDBOX/housepi/inst/models/augmented_infection_time/housepiaugmentedinfectiontime.so"

choice <- "augmented_infection_time"

for (discr in c(FALSE, TRUE, FALSE, TRUE)) {
  print(paste0("MODEL DISCR=", discr))
  hhgen <- housepi:::HousepiGenerator$new(
    choice, discr = discr, quiet = TRUE)
  private <- hhmodel$.__enclos_env__$private
  print(private$.dll)
  print(private$.last_user_flags)
  
  m <- housepi:::get_rcpp_module_from_dll(official)
  y <- m$HousepiRcppCompiledFlags$new()
  print(paste0("SO=", y$get_flags()))
  print("")
}


