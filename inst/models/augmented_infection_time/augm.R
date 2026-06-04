list(
  infec_time = ~ onset_time@A - incub_period@A,
  onset_time = list(prior = "normal",
                    mean = ~ onset_date@D),
  incub_period = list(prior = "gamma",
                      shape = 2, scale = 1,
                      max = 14)
)
