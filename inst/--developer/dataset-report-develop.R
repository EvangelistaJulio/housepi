
hhmodel <- housepi::HouseholdModel$new("chain-binomial")
hhmodel$load_param_template()
hhmodel$loglikelihood(housepi::survey)

hhmodel$.__enclos_env__$private$.model$.__enclos_env__$private$.rcpp$export_data_report()

hhmodel$.__enclos_env__$private$.model$.__enclos_env__$private$.rcpp$export_house_size_report()
