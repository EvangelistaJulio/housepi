# Print behavior works as expected

    Code
      print(hhmodel)
    Message
      i Household model chain binomial (template).
      v Parameters definitions are complete.
      ! Group rates:
      There are 3 missing group rate definitions: 'community', 'infectivity' and 'susceptibility'.
      They can be specified for model definitions.

---

    Code
      print(hhmodel2)
    Message
      i Household model chain binomial (template).
      x Parameters:
      There are 4 missing parameter definitions: 'alpha', 'beta', 'gamma' and 'size_ref'.
      They are required by model definitions.
      ! Group rates:
      There are 3 missing group rate definitions: 'community', 'infectivity' and 'susceptibility'.
      They can be specified for model definitions.

---

    Code
      print(new_housepi_model(model = "blank"))
    Message
      i Household model blank (template).

# Nested definitions works as expected

    Code
      print(hhmodel)
    Message
      i Household model blank (template).
      v Augmented definitions are complete.

