# Print behavior works as expected

    Code
      print(new_housepi_model(model = "test"))
    Output
      i Household model "test" (`template`).

---

    Code
      print(new_housepi_model(model = "final_size"))
    Output
      i Household model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "hhsize_ref".
        They are required by model definitions.

---

    Code
      print(new_housepi_model(flu_1918_frederick, "final_size"))
    Output
      i Household model "final size" (`template`).
      x Parameters: there are 4 missing parameter definitions: "alpha", "beta",
        "gamma" and "hhsize_ref".
        They are required by model definitions.

---

    Code
      print(add_parameter(new_housepi_model(flu_1918_frederick, "final_size"),
      import = "template"))
    Output
      i Household model "final size" (`template`).
      v Parameters

