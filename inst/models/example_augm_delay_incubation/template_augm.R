list(
  exposure = list(
    x = "exposure_min@D",
    model = "unif", args = list(min = "exposure_min@D", max = "exposure_max@D"),
    rule = "norm", rule_args = list(sd = 0.25)
  ),
  onset = list(
    x = "onset_max@D",
    model = "unif", args = list(min = "onset_min@D", max = "onset_max@D"),
    rule = "norm", rule_args = list(sd = 0.25)
  )
)