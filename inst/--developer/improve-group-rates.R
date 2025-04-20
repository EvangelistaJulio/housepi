
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

library(tidyverse)
# FIXME: I want to be able to assign group rates to this kind og grouping
hhmodel$dataset <- dataset %>% as_tibble() %>%
  mutate(agegroup = cut(age, breaks = c(0, 12, 25, 50, 100),
                        right = FALSE, include.lowest = TRUE))
