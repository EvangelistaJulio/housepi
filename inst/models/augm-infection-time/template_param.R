list(
  alpha = list(
    value = 0.0001, prior = "unif", args = c(min = 0, max = 1),
    rule = "lnorm"
  ),
  beta = list(
    value = 0.07, prior = "unif", args = c(min = 0, max = 1),
    rule = "lnorm"
  ),
  gamma = list(
    value = 0.9, prior = "unif", args = c(min = 0, max = 4),
    rule = "norm"
  ),
  size_ref = 4,
  study_period = 45
)
