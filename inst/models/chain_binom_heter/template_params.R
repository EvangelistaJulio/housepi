list(
  alpha = list(
    value = 0.15, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  beta = list(
    value = 0.22, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_adult_child = list(
    value = 0.75, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_child_adult = list(
    value = 1.1, prior = "unif",  args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_child_child = list(
    value = 1.25, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  inf_child = list(
    value = 1.25, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  sus_child = list(
    value = 1.3, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  sus_women = list(
    value = 1.1, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  child_age = list(value = 12)
)
