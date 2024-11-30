list(
  alpha = list(
    value = 0.15, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  beta = list(
    value = 0.22, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  inf_child = list(
    value = 1.25, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  comsus_child = list(
    value = 1.3, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  comsus_woman = list(
    value = 1.0, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  sus_child = list(
    value = 0.9, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  sus_woman = list(
    value = 1.2, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  child_age = list(value = 12)
)
