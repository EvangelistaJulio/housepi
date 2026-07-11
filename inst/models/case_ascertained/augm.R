list(
  infec_time = list(
    prior = "uniform",
    min = ~ start_date@D - 10, max = ~ start_date@D + study_period@D
  )
)
