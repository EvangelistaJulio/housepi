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
      There are 4 missing parameter definitions: 'alpha', 'beta', 'gamma', 'size_ref'.
      They are required by model definitions.

---

    Code
      print(new_housepi_model("blank"))
    Message
      i Household model blank (template).

