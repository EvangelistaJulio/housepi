list(
  infec_time = list(prior = "uniform",
                    min = -14, max = ~ end_followup@D,
                    proposal = "normal")
)
