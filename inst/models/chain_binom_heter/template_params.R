list(
  alpha = list(
    x = 0.15, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  beta = list(
    x = 0.22, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_adult_child = list(
    x = 0.75, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_child_adult = list(
    x = 1.1, prior = "unif",  args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  cont_child_child = list(
    x = 1.25, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  inf_child = list(
    x = 1.25, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  sus_child = list(
    x = 1.3, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  sus_women = list(
    x = 1.1, prior = "unif", args = list(min = 0, max = 4),
    rule = "lnorm", rule_args = list(sdlog = 1)
  ),
  child_age = list(x = 12)
)
