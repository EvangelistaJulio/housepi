# Exposed object print behavior is stable

    Code
      print(expose_function("pow", x2 = 3))
    Message
      i Exposed function "pow".
      `x2` = 3

---

    Code
      print(expose_distribution("normal", mean = 2, sd = 3))
    Message
      i Exposed distribution "normal".
      `mean` = 2
      `sd` = 3

---

    Code
      print(expose_infec_profile("gamma", shape = 2, scale = 1))
    Message
      i Exposed infectivity profile "gamma".
      `shape` = 2
      `scale` = 1

---

    Code
      print(expose_function("sin"))
    Message
      i Exposed function "sin".

