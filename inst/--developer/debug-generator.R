library(housepi)

model_dir <- "../debugger"
template_household_model("chain-binomial", model_dir)
hhmodel <- HouseholdModel$new(model_dir,
                              as_template = T)
hhmodel$traces_results
