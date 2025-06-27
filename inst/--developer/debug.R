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

# CHECK:
# I shall not pollute the registrations
getDLLRegisteredRoutines("housepi")$`.Call`
getDLLRegisteredRoutines("housepienvironment")$`.Call`
getLoadedDLLs()