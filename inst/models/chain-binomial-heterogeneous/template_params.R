list(
  alpha = list(
    value = 0.15, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_scale = 1
  ),
  beta = list(
    value = 0.22, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_scale = 1
  ),
  inf_child = list(
    value = 1.25, prior = "lnorm", args = list(sdlog = 3),
    rule = "lnorm", rule_scale = 1
  ),
  com_sus_child = list(
    value = 1.3, prior = "lnorm", args = list(sdlog = 3),
    rule = "lnorm", rule_scale = 1
  ),
  com_sus_woman = list(
    value = 1.0, prior = "lnorm", args = list(sdlog = 3),
    rule = "lnorm", rule_scale = 1
  ),
  sus_child = list(
    value = 0.9, prior = "lnorm", args = list(sdlog = 3),
    rule = "lnorm", rule_scale = 1
  ),
  sus_woman = list(
    value = 1.2, prior = "lnorm", args = list(sdlog = 3),
    rule = "lnorm", rule_scale = 1
  ),
  child_age = list(value = 12)
)
