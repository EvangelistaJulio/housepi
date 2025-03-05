list(
  onset_date = list(
    value = "onset_date@D",
    rule = "norm", rule_scale = 0.5,
    rule_from = "onset_date@D"
  ),
  incub_period = list(
    init_rand = TRUE, model = "gamma",
    args = c(shape = "inc_shape@P", scale = "inc_scale@P"),
    rule = "gibbs", max = 14
  )
)