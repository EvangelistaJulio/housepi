
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

hhmodel <- HouseholdModel$new("chain-binomial-heterogeneous")
dataset <- housepi::survey
dataset$child <- as.integer(dataset$age < 15)
hhmodel$dataset <- dataset
hhmodel$load_param_template()

hhmodel$add_group_rate("susceptibility", by = c("sex", "child"), independent = TRUE)
hhmodel$.__enclos_env__$private$.build_group_rate(NULL)
hhmodel$.__enclos_env__$private$.group_rate$logic_tables
hhmodel$param_defs

# FIXME: successive calls to .build_group_rate won't compile...
# FIXME: old parameters are not dropped
hhmodel$add_group_rate("susceptibility", by = c("sex", "child"), independent = FALSE)
hhmodel$.__enclos_env__$private$.build_group_rate(NULL)
hhmodel$.__enclos_env__$private$.group_rate$logic_tables
hhmodel$param_defs


#hhmodel$loglikelihood(seed = 1)

