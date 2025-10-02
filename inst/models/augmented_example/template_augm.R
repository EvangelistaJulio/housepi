list(
  incub_period = ~onset_time@A - infec_time@A,
  infec_time = list(
    value = ~exposure_min@D,
    model = "uniform", args = c(min = ~exposure_min@D, max = ~exposure_max@D),
    rule = "normal", scale = 0.5
  ),
  onset_time = list(
    value = ~onset_max@D,
    model = "uniform", args = c(min = ~onset_min@D, max = ~onset_max@D),
    rule = "normal", scale = 0.5
  ),
  proba = list(
    func = "gamma",
    func_args = c(x = ~incub_period@A,  shape = ~shape@P, scale = ~scale@P)
  )
)
