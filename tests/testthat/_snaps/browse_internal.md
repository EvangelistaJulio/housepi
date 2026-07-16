# Available names listing is stable

    Code
      browse_functions()
    Message
      i Available functions:
      * "ceil"
      * "cos"
      * "div"
      * "exp"
      * "floor"
      * "log"
      * "max"
      * "min"
      * "minus"
      * "pow"
      * "prod"
      * "round"
      * "sin"
      * "sqrt"
      * "sub"
      * "sum"
      * "tan"

---

    Code
      browse_distributions()
    Message
      i Available distributions:
      * "beta"
      * "binomial"
      * "categorical"
      * "exponential"
      * "gamma"
      * "geometric"
      * "lognormal"
      * "negative_binomial"
      * "normal"
      * "poisson"
      * "rayleigh"
      * "uniform"
      * "weibull"

---

    Code
      browse_infectivity_profiles()
    Message
      i Available infectivity profiles:
      * "beta"
      * "binomial"
      * "categorical"
      * "exponential"
      * "ferretti20"
      * "gamma"
      * "geometric"
      * "lognormal"
      * "negative_binomial"
      * "normal"
      * "poisson"
      * "rayleigh"
      * "uniform"
      * "weibull"

# Exposed object print behavior is stable

    Code
      print(browse_functions("pow", x2 = 3))
    Output
      i Exposed function "pow".
      `x2` = 3

---

    Code
      print(browse_distributions("normal", mean = 2, sd = 3))
    Output
      i Exposed distribution "normal".
      `mean` = 2
      `sd` = 3

---

    Code
      print(browse_infectivity_profiles("gamma", shape = 2, scale = 1))
    Output
      i Exposed infectivity profile "gamma".
      `shape` = 2
      `scale` = 1

---

    Code
      print(browse_functions("sin"))
    Output
      i Exposed function "sin".

