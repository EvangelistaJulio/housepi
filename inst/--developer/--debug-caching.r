
devtools::load_all()

# check module cache is empty
env <- environment(get_rcpp_module_from_dll)
names(env$module_cache)

choice_1 <- "augmented_infection_time"
# load base DLL
dll_1 <- file.path(template_models_dir(), template_dll_file(choice_1))
# register module and return Model 1
m_1 <- get_rcpp_module_from_dll(dll_1)$HousepiRcppModel$new()
# method I know as different return
m_1$export_infec_profile()

choice_2 <- "chain_binomial"
# load base DLL
dll_2 <- file.path(template_models_dir(), template_dll_file(choice_2))
# register module and return Model 0
m_2 <- get_rcpp_module_from_dll(dll_2)$HousepiRcppModel$new()
m_2$export_infec_profile()

# check available modules
env <- environment(get_rcpp_module_from_dll)
names(env$module_cache)

# return base Model again
m_1 <- get_rcpp_module_from_dll(dll_1)$HousepiRcppModel$new()
m_1$export_infec_profile()

# double check
m_3 <- env$module_cache[[dll_1]]$HousepiRcppModel$new()
m_3$export_infec_profile()
