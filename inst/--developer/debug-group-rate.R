
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

hhmodel$load_group_rate_list(
  susceptibility = list(by = c("sex", "child"), independent = TRUE)
)
hhmodel$group_rate_tables
hhmodel$param_defs


hhmodel$loglikelihood(seed = 1)

hhmodel_2 <- hhmodel$clone(deep = TRUE)
hhmodel_2$add_group_rate("susceptibility", by = c("sex", "child"), independent = FALSE)
hhmodel_2$loglikelihood(seed = 1)

# FIXME: current caching system works, but...
# it overwrites everythin in tempdir()/housepichainbinomialheterogeneousgrouprates
# meaning that old model needs to be recompiled...
hhmodel$loglikelihood(seed = 1)
