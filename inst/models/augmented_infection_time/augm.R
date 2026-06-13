list(
  infec_time = list(prior = "uniform",
                    min = ~ start_followup@D - 14,
                    max = ~ end_followup@D,
                    proposal = "normal")
)
