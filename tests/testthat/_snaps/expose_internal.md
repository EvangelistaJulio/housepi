# Exposed object print behavior is stable

    Code
      print(expose_function("pow", x2 = 3))
    Output
      i Exposed function "pow".
      `x2` = 3

---

    Code
      print(expose_distribution("normal", mean = 2, sd = 3))
    Output
      i Exposed distribution "normal".
      `mean` = 2
      `sd` = 3

---

    Code
      print(expose_infectivity_profile("gamma", shape = 2, scale = 1))
    Output
      i Exposed infectivity profile "gamma".
      `shape` = 2
      `scale` = 1

---

    Code
      print(expose_function("sin"))
    Output
      i Exposed function "sin".

