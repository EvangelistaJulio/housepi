list(
  alpha = list(
    value = 0.0001, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_scale = 1
  ),
  beta = list(
    value = 0.07, prior = "unif", args = list(min = 0, max = 1),
    rule = "lnorm", rule_scale = 1
  ),
  gamma = list(
    value = 0.85, prior = "unif", args = list(min = 0, max = 4),
    rule = "norm", rule_scale = 1
  ),
  size_ref = list(value = 4),
  study_period = list(value = 50)
)
