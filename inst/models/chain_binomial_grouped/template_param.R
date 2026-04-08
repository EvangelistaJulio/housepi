list(
  alpha = list(value = 0.15,
               prior = "uniform",
               proposal = "lognormal"),
  beta = list(value = 0.20,
              prior = "uniform",
              proposal = "lognormal"),
  gamma = list(value = 1,
               prior = "uniform",
               max = 4,
               proposal = "normal"),
  size_ref = 4
)
