devtools::load_all()

# check module cache is empty
env <- environment(get_rcpp_module_from_dll)
names(env$module_cache)

# load base DLL
base_dll <- "/Users/fcamagli/Documents/SANDBOX/housepi/inst/models/augmented_infection_time/housepiaugmentedinfectiontime.so"
m <- get_rcpp_module_from_dll(base_dll)

# check available modules
env <- environment(get_rcpp_module_from_dll)
names(env$module_cache)

# return base Model
env$module_cache[[base_dll]]$HousepiRcppModel$new()

# compile with DISCRETE flags and return store Model internally
hhgen <- HousepiGenerator$new("augmented_infection_time", discr = TRUE, quiet = TRUE)

# check available modules
env <- environment(get_rcpp_module_from_dll)
names(env$module_cache)

# return base Model again
env$module_cache[[base_dll]]$HousepiRcppModel$new()
