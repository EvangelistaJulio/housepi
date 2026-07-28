# Print behavior works as expected

    Code
      print(new_housepi_model(model = "test"))
    Output
      i Household transmission model "test" (`template`).

---

    Code
      print(new_housepi_model(model = "final_size"))
    Output
      i Household transmission model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "hhsize_ref".
        Required by transmission model.

---

    Code
      print(new_housepi_model(flu_1918_frederick, "final_size"))
    Output
      i Household transmission model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "hhsize_ref".
        Required by transmission model.

---

    Code
      print(add_parameter(new_housepi_model(flu_1918_frederick, "final_size"),
      import = "template"))
    Output
      i Household transmission model "final size" (`template`).
      v Parameters

