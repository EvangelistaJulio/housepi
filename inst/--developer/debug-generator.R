
source("R/utils.R")
source("R/template.R")
source("R/augmented.R")
source("R/parameter.R")
source("R/infectivity.R")
source("R/group_rate.R")
source("R/expected_cases.R")
source("R/model_generator.R")
source("R/model_class.R")

# create custom package
model_choice <- "../debugger"
#template_household_model("chain-binomial", model_choice)

# # FIXME : dyn.load(private$.dll) is WRONG
# hhmodel <- HousepiGenerator$new(model_choice, as_template = TRUE)
#
# hhmodel$.__enclos_env__$private$.need_compile
# hhmodel$cache_dir
# hhmodel$.__enclos_env__$private$.dll
#
# hhmodel$generate()
#
# hhmodel$.__enclos_env__$private$.need_compile
# hhmodel$cache_dir
# hhmodel$.__enclos_env__$private$.dll
#
# # second round
# hhmodel <- HousepiGenerator$new(model_choice, as_template = TRUE)
#
# hhmodel$.__enclos_env__$private$.need_compile
# hhmodel$cache_dir
# hhmodel$.__enclos_env__$private$.dll
#
# hhmodel$generate()

# load or compile custom package
hhmodel <- HouseholdModel$new(model_choice,
                              as_template = TRUE, force = FALSE)

# test correct version is loaded
hhmodel$load_param_template()
hhmodel$loglikelihood(housepi::survey)

