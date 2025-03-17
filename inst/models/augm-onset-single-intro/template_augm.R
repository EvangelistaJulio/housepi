list(
  onset_date = list(
    value = "onset_date@D",
    rule = "norm", scale = 0.5,
    from = "onset_date@D"
  ),
  incub_period = list(
    init_rand = TRUE, model = "reparam_gamma",
    args = c(mean = 2, sd = 1),
    max = 14, rule = "gibbs"
  )
)
