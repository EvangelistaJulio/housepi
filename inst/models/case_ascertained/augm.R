list(
  infec_time = list(
    model = "uniform",
    min = ~ start_date@D - 10, max = ~ end_date@D
  )
)
