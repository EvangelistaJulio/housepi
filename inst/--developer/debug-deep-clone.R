library(testthat)

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

# load or compile custom package
hhmodel <- HouseholdModel$new("augm-onset-single-intro")

# FIXME early clone breaks it
hhclone <- hhmodel$clone(deep = TRUE)

hhmodel$load_param_template()
hhmodel$load_augm_template()
hhmodel$loglikelihood(housepi::survey)

hhclone <- hhmodel$clone(deep = TRUE)

# definitions
expect_equal(hhmodel$param_defs, hhclone$param_defs)
expect_equal(hhmodel$augm_defs, hhclone$augm_defs)
expect_equal(hhmodel$setup, hhclone$setup)

# dataset
expect_equal(hhmodel$dataset, hhclone$dataset)

expect_equal(hhmodel$show_summary(), hhclone$show_summary())
