# A model prints with and without data and parameters

    Code
      print(new_housepi_model(model = "test"))
    Message
      i Household transmission model "test" (`template`).

---

    Code
      print(new_housepi_model(model = "final-size"))
    Message
      i Household transmission model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma", and "hh_size_ref".
        Required by transmission model.

---

    Code
      print(new_housepi_model(flu_1918_frederick, "final-size"))
    Message
      i Household transmission model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma", and "hh_size_ref".
        Required by transmission model.

---

    Code
      print(hh_model)
    Message
      i Household transmission model "final size" (`template`).
      v Parameters

