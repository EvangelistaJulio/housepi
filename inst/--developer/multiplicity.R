
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
hhmodel$load_param_template()

dataset <- housepi::survey
dataset$agegroup <- as.integer(dataset$age <= 15)
hhmodel$dataset <- dataset

hhmodel$add_group_rate("susceptibility", by = c("sex", "agegroup"),
                       independent = FALSE, reference = list(sex = 2, agegroup = 0))
hhmodel$group_rate_tables

# multiplicity engine #

# dataset required by the implementation
store <- hhmodel$.__enclos_env__$private$.require_dataset()
# dataset required by group rate
for (this_rate in hhmodel$.__enclos_env__$private$.group_rate$logic_tables) {
  for (this_table in this_rate) {
    store <- c(store, names(this_table))
  }
}
store <- setdiff(store, c(".counts", ".param.name"))


library(tidyverse)

# household multiplicities table
multipl_hh_table <- dataset %>%
  select(house_id, status, !!store) %>%
  nest(data = -house_id) %>%
  mutate(
    # sort each group's individuals to make order irrelevant
    sorted_data = map(data, ~ arrange(.x, across(everything()))),
    # convert each group's data to a string representation
    hash = map_chr(sorted_data, ~ paste(apply(.x, 1, paste, collapse = "-"), collapse = "|"))
  ) %>%
  group_by(hash) %>%
  summarise(
    class_id = first(house_id),
    multiplicity = n(),
    full_class = list(house_id),
    .groups = "drop"
  )

output <- data.frame(house_id = unique(dataset$house_id))
rownames(output) <- output$house_id
output$multiplicity <- 0
output$class_id <-NA

output[ multipl_hh_table$class_id, ]$multiplicity <- multipl_hh_table$multiplicity

for (i in seq_len(nrow(multipl_hh_table))) {
  output[multipl_hh_table$full_class[i][[1]], ]$class_id <- multipl_hh_table$class_id[i]
}
output
