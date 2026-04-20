# Print behavior works as expected

    Code
      print(hhmodel)
    Message
      i Household model chain binomial (template).
      v Parameters definitions are complete.

---

    Code
      print(hhmodel2)
    Message
      i Household model chain binomial (template).
      x Parameters:
      There are 4 missing parameter definitions: 'alpha', 'beta', 'gamma' and 'size_ref'.
      They are required by model definitions.

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

