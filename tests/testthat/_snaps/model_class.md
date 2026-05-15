# Print behavior works as expected

    Code
      print(hhmodel)
    Message
      i Household model 'chain binomial' (`template`).
      v Parameters: definitions complete.

---

    Code
      print(hhmodel2)
    Message
      i Household model 'chain binomial' (`template`).
      x Parameters: there are 4 missing parameter definitions:
      'alpha', 'beta', 'gamma' and 'size_ref'.
      They are required by model definitions.

---

    Code
      print(new_housepi_model(model = "blank"))
    Message
      i Household model 'blank' (`template`).

---

    Code
      print(hhmodel_d)
    Message
      i Household model 'augmented infection time' (`template`).
      Model is set to discrete time mode.
      v Data: definitions complete.
      x Parameters: there are 4 missing parameter definitions:
      'alpha', 'beta', 'gamma' and 'size_ref'.
      They are required by model definitions.
      x Augmented: there is 1 missing augmented definition: 'infec_time'.
      It is required by model definitions.
      ! Group rates: there are 3 missing group rate definitions:
      'community', 'infectivity' and 'susceptibility'.
      They can be specified for model definitions.
      ! Contact rates: using default contact rates.
      A custom choice can be provided.
      ! Infectivity profile: using templated profile. A custom choice can be provided.

---

    Code
      print(hhmodel_ip)
    Message
      i Household model 'augmented infection time' (`template`).
      v Data: definitions complete.
      x Parameters: there are 4 missing parameter definitions:
      'alpha', 'beta', 'gamma' and 'size_ref'.
      They are required by model definitions.
      x Augmented: there is 1 missing augmented definition: 'infec_time'.
      It is required by model definitions.
      ! Group rates: there are 3 missing group rate definitions:
      'community', 'infectivity' and 'susceptibility'.
      They can be specified for model definitions.
      ! Contact rates: using default contact rates.
      A custom choice can be provided.
      ! Infectivity profile: using templated profile. A custom choice can be provided.

# Nested definitions works as expected

    Code
      print(hhmodel)
    Message
      i Household model 'blank' (`template`).
      v Augmented: definitions complete.

