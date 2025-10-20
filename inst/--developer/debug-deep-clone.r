source("inst/--developer/debug.R")

hhsimul <- new_housepi_model("blank")
hhsimul$add_data_wide(dplyr::mutate(housepi::survey, status = NA))
head(hhsimul$data_wide)
hhsimul$.__enclos_env__$private$.rcpp

x <- hhsimul$clone(deep = TRUE)
# FIXME this should not be the same!
x$.__enclos_env__$private$.rcpp

x$add_data_wide(housepi::survey)
# FIXME this should not change!
head(hhsimul$data_wide)

#NOTE : this is different
# coso <- new_housepi_model("blank")
# coso$.__enclos_env__$private$.rcpp

# dll <- "/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/housepi/template/src/housepienvironment.so"
# m <- housepi:::get_rcpp_module_from_dll("housepi_model_core", dll)
# m$ModelHousepi$new()
