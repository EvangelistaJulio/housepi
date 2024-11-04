list(
  alpha = list(
    value = 0.15, prior = "unif", args = list(min = 0, max = 5), 
    rule = "lnorm", rule_args = list(sdlog = 1)
  ), 
  beta = list(
    value = 0.25, prior = "unif", args = list(min = 0, max = 5), 
    rule = "lnorm", rule_args = list(sdlog = 1)
  )
)
