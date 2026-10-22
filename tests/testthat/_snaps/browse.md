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
    Message
      i Exposed function "pow".
        x2: 3
        Domain: continuous on (-Inf, Inf)

---

    Code
      print(browse_distributions("normal", mean = 2, sd = 3))
    Message
      i Exposed distribution "normal".
        mean: 2
        sd: 3
        Domain: continuous on (-Inf, Inf)

---

    Code
      print(browse_infectivity_profiles("gamma", shape = 2, scale = 1))
    Message
      i Exposed infectivity profile "gamma".
        shape: 2
        scale: 1
        Domain: continuous on [0, Inf)

---

    Code
      print(browse_functions("sin"))
    Message
      i Exposed function "sin".
        Domain: continuous on (-Inf, Inf)

---

    Code
      print(browse_distributions("poisson", lambda = 3, min = 0.5))
    Message
      i Exposed distribution "poisson".
        lambda: 3
        Domain: discrete on [0.5, Inf)

---

    Code
      print(browse_distributions("poisson", lambda = 3, loc = 2))
    Message
      i Exposed distribution "poisson".
        lambda: 3
        loc: 2
        Domain: discrete on [2, Inf)

---

    Code
      print(browse_distributions("gamma", shape = 2, scale = 1, loc = 3))
    Message
      i Exposed distribution "gamma".
        shape: 2
        scale: 1
        loc: 3
        Domain: continuous on [3, Inf)

# browse_templates branches inform correctly

    Code
      browse_templates()
    Message
      i Available template models:
      * `case_ascertained`
      * `final_size`
      * `test`

