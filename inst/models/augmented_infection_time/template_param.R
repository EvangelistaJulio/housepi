list(
  alpha = list(value = 1.e-4,
               prior = "uniform",
               proposal = "lognormal"),
  beta = list(value = 0.07,
              prior = "uniform",
              proposal = "lognormal"),
  gamma = list(value = 0.9,
               prior = "uniform",
               max = 4,
               proposal = "normal"),
  size_ref = 4,
  study_period = 45
)
