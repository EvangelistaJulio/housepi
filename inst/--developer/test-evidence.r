
library(tidyverse)

dataset <- data.frame(house_id = 1, id = 1, status = 1, x = 0)

# Exponential

expected <- function(x) log(1 + 1 / x) - 1 / (x + 1)
log(expected(dataset$x))

tibble(x = seq(0.01, 10, 0.01),
       y = log(expected(x))) |>
  ggplot(aes(x = x, y = y)) +
  geom_line()

results <- new_housepi_model("augmented_example") |>
  add_parameter("theta", prior = "uniform") |>
  load_augmented_list(
    .hidden = list(model = "exponential", rate = theta@P,
                   rule = "sample"),
    proba = list(func = "exponential",
                 func_args = c(x = x@D, rate = .hidden))
  ) |>
  evidence(dataset, steps = 5000, rungs = 100, chains = 3)
results$estimate

results$detail |>
  as_tibble() |>
  ggplot(aes(x = inv_temper, y = mean_loglike)) +
  geom_line()

# Gaussian

expected <- function(x) stats::dnorm(x, mean = 0, sd = sqrt(3))
log(expected(dataset$x))

tibble(x = seq(-10, 10, 0.01),
       y = log(expected(x))) |>
  ggplot(aes(x = x, y = y)) +
  geom_line()

results <- new_housepi_model("augmented_example") |>
  add_parameter("theta", prior = "normal") |>
  load_augmented_list(
    .hidden = list(model = "normal", mean = theta@P,
                   rule = "sample", min = -10, max = 10),
    proba = list(func = "normal",
                 func_args = c(x = x@D, mean = .hidden))
  ) |>
  evidence(dataset, steps = 10000, rungs = 100, chains = 4)
results$estimate + log(20)
log(expected(dataset$x))

results$detail |>
  as_tibble() |>
  ggplot(aes(x = inv_temper, y = mean_loglike)) +
  geom_line()
