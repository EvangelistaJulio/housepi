# source loads
source("R/utils.R")
source("R/template.R")
source("R/augmented.R")
source("R/parameter.R")
source("R/infectivity.R")
source("R/group_rate.R")
source("R/expected_cases.R")
source("R/model_generator.R")
source("R/model_wrapper.R")
source("R/model_class.R")
source("R/sampler_class.R")

hhmodel <- HouseholdModel$new("chain-binomial-heterogeneous")

# load parameter definitions  from template
hhmodel$load_param_template()

dataset <- housepi::survey
dataset$agegroup <- ifelse(dataset$age <= 15, "child", "adult")

# import dataset from active field
hhmodel$dataset <- dataset

# compact definition from list
hhmodel$load_group_rate_list(
  infectivity = list(by = "agegroup"),
  susceptibility = list(by = c("agegroup", "sex"), independent = TRUE)
)

private <- hhmodel$.__enclos_env__$private

# THIS is WRONG
# 1: community is missing
# 2: susceptibiltiy has multiple 
hhmodel$group_rate_defs

private$.export_to_cpp_standalone()
