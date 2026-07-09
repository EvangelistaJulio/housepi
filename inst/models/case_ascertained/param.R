list(
  alpha = list(value = 1.e-3,
               prior = "uniform",
               proposal = "lognormal"),
  beta = list(value = 0.15,
              prior = "uniform",
              proposal = "lognormal"),
  gamma = list(value = 1,
               prior = "uniform",
               max = 4,
               proposal = "normal"),
  hhsize_ref = 4
)
