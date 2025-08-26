#library(testthat)

options(error = NULL)
source("inst/--developer/debug.R")

# debug deep clone of housepi wrapper class
hhwrap <- HousepiWrapper$new("chain-binomial")
hhclone <- hhwrap$clone(deep = TRUE)
hhwrap$generate()
hhclone <- hhwrap$clone(deep = TRUE)

# devug deep clone of final hosuepi model class
hhmodel <- HouseholdModel$new("chain-binomial")
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
