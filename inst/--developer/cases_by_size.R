library(housepi)

OUTPUT="vignettes/debug"

dataset = housepi::survey
size = 500
n_burnin = 0
n_adapt = 2500
thin = 5
seed = 82754390

################################################################################
################################################################################
################################################################################

cb_params <-  list(
  alpha = list(prior = "unif", args = list(min = 0, max = 5),
               rule = "lnorm", rule_scale = 1),
  beta = list(prior = "unif", args = list(min = 0, max = 5),
              rule = "lnorm", rule_scale = 1)
)

hmc1 <- HouseholdModelClass$new("chain-binomial")
hmc1$load_param_list(cb_params)

res1 <- retrieve_expected_cases_by_size(
  hmc1,
  dataset = dataset,
  size = size,
  n_burnin = n_burnin,
  thin = thin, 
  n_adapt = n_adapt,
  seed = seed
)

t1 <- res1$traces
pdf(file.path(OUTPUT, "homogeneous-traces.pdf"))
print(lattice::xyplot(t1))
dev.off()

write.csv(retrieve_traces_results(t1), 
          file.path(OUTPUT, "homogeneous-traces.csv"),
          quote = F, row.names = F)

ft1 <- show_model_check_table(res1)
flextable::save_as_image(ft1, file.path(OUTPUT, "homogeneous-table.png"))

################################################################################
################################################################################
################################################################################

hmc2 <- HouseholdModelClass$new("chain-binomial-heterogeneous")
hmc2$load_param_list(c(cb_params, list(
  inf_child = 1., com_sus_child = 1.,
  com_sus_woman = 1.0, sus_child = 1.0, 
  sus_woman = 1., child_age = 12)))

res2 <- retrieve_expected_cases_by_size(
  hmc2,
  dataset = dataset,
  size = size,
  n_burnin = n_burnin,
  thin = thin, 
  n_adapt = n_adapt,
  seed = seed
)

t2 <- res2$traces
pdf(file.path(OUTPUT, "heterogeneous-traces.pdf"))
print(lattice::xyplot(t2))
dev.off()

write.csv(retrieve_traces_results(t2), 
          file.path(OUTPUT, "heterogeneous-traces.csv"),
          quote = F, row.names = F)

ft2 <- housepi::show_model_check_table(res2)
flextable::save_as_image(ft2, file.path(OUTPUT, "heterogeneous-table.png"))
