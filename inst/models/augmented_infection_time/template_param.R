list(
  alpha = list(
    value = 0.0001, prior = "uniform", rule = "lognormal"
  ),
  beta = list(
    value = 0.07, prior = "uniform", rule = "lognormal"
  ),
  gamma = list(
    value = 0.9, prior = "uniform", args = c(max = 4),
    rule = "normal"
  ),
  size_ref = 4,
  study_period = 45
)
