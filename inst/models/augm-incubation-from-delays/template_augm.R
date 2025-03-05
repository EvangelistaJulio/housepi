list(
  exposure = list(
    value = "exposure_min@D",
    model = "unif", args = c(min = "exposure_min@D", max = "exposure_max@D"),
    rule = "norm", rule_scale = 0.5
  ),
  onset = list(
    value = "onset_max@D",
    model = "unif", args = c(min = "onset_min@D", max = "onset_max@D"),
    rule = "norm", rule_scale = 0.5
  )
)