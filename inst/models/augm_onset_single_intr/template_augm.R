list(
  onset_date = list(
    value = "onset_date@D",
    rule = "norm", rule_args = list(sd = 0.25), rule_from = "onset_date@D"
  ),
  incub_period = list(
    init_rand = TRUE, model = "gamma", args = list(shape = "inc_shape@P", scale = "inc_scale@P"),
    rule = "gibbs", range = list(max = 14)
  )
)