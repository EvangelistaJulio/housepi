list(
  beta = list(
    x = 0.05, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  gamma = list(
    x = 0.85, prior = "unif", args = list(min = 0, max = 4),
    rule = "norm", rule_args = list(sd = 1)
  ),
  sus_vacc = list(
    x = 0.15, prior = "unif", args = list(min = 0, max = 5),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  inf_vacc = list(
    x = 0.4, prior = "unif", args = list(min = 0, max = 5),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_adult_child = list(x = 0.75),
  cont_child_adult = list(x = 1.1),
  cont_child_child = list(x = 1.25),
  study_period = list(x = 60),
  child_age = list(x = 15)
)
