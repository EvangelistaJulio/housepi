# source loads
source("R/utils.R")
source("R/template.R")
source("R/augmented.R")
source("R/parameter.R")
source("R/infectivity.R")
source("R/group_rates.R")
source("R/expected_cases.R")
source("R/model_generator.R")
source("R/model_class.R")
source("R/sampler_class.R")
source("R/expose_class.R")
source("R/expose_internal.R")
source("R/posterior.R")

## CHECK:
##----------
# getDLLRegisteredRoutines("housepi")$`.Call`
# getDLLRegisteredRoutines("housepienvironment")$`.Call`
# getLoadedDLLs()
