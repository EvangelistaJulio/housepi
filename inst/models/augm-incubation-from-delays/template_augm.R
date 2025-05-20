list(
  infec_time = list(
    value = "exposure_min@D",
    model = "unif", args = c(min = "exposure_min@D", max = "exposure_max@D"),
    rule = "norm", scale = 0.5
  ),
  onset_time = list(
    value = "onset_max@D",
    model = "unif", args = c(min = "onset_min@D", max = "onset_max@D"),
    rule = "norm", scale = 0.5
  )
  incub_period = "onset_time@A - infec_time@A"
)