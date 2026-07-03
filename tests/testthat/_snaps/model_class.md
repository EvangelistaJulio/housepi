# Print behavior works as expected

    Code
      print(hhmodel)
    Message
      i Household model "final size" (`template`).
      v Parameters

---

    Code
      print(hhmodel2)
    Message
      i Household model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "size_ref".
        They are required by model definitions.

---

    Code
      print(new_housepi_model(model = "blank"))
    Message
      i Household model "blank" (`template`).

---

    Code
      print(hhmodel_d)
    Message
      i Household model "case ascertained" (`template`).
      Model is set to discrete time mode.
      x Data: there is 1 missing data definition "start_date".
        It is required by model definitions.
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "size_ref".
        They are required by model definitions.
      x Augmented: there is 1 missing augmented definition "infec_time".
        It is required by model definitions.
      ! Group rates: using default group rates.
        Available: "community", "susceptibility" and "infectivity".
      ! Contact rates: using default contact rates.
      ! Infectivity profile: using templated profile.
        A custom choice can be provided.

---

    Code
      print(hhmodel_ip)
    Message
      i Household model "case ascertained" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "size_ref".
        They are required by model definitions.
      x Augmented: there is 1 missing augmented definition "infec_time".
        It is required by model definitions.
      ! Group rates: using default group rates.
        Available: "community", "susceptibility" and "infectivity".
      ! Contact rates: using default contact rates.
      ! Infectivity profile: using templated profile.
        A custom choice can be provided.
      v Data

# Nested definitions works as expected

    Code
      print(hhmodel)
    Message
      i Household model "blank" (`template`).
      v Augmented

