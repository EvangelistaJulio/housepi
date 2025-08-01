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
source("R/expose.R")

original <- new_housepi_model("chain-binomial-heterogeneous") |>
  load_parameter_template()
original$dataset <- housepi::survey
 
unlist(original$.__enclos_env__$private$.group_rates$defs$susceptibility)

# FIXME cloning is not working
cloned <- original$clone(deep = TRUE)
cloned$add_group_rate("susceptibility", by = "sex")
unlist(original$.__enclos_env__$private$.group_rates$defs$susceptibility)

# Check if .group_rates itself is cloned
original_gr <- original$.__enclos_env__$private$.group_rates
cloned_gr <- cloned$.__enclos_env__$private$.group_rates

# Check memory addresses
cat("Original .group_rates address:", data.table::address(original_gr), "\n")
cat("Cloned .group_rates address:", data.table::address(cloned_gr), "\n")
cat("Original defs address:", data.table::address(original_gr$defs), "\n")
cat("Cloned defs address:", data.table::address(cloned_gr$defs), "\n")

# # Are they the same R6 object?
# cat("Same .group_rates object?", identical(original_gr, cloned_gr), "\n")
# 
# # Check the defs field specifically
# cat("Same defs object?", identical(original_gr$defs, cloned_gr$defs), "\n")
