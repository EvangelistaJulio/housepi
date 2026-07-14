list(
  infec_time = list(
    prior = "uniform",
    min = ~ start_date@D - 10, max = ~ end_date@D
  )
)
