library(housepi)
library(ggplot2)

hhsimul <- new_housepi_model("augmented_infection_time", discr = TRUE) |>
  load_parameter_template() |>
  load_augmented_list(
    infec_time = list(model = "lognormal", discr = TRUE,
                      args = c(meanlog = 4, sdlog = 0.5),
                      max = 365),
    incub_period = list(model = "gamma", discr = TRUE,
                        args = c(shape = ~inc_shape@P,  scale = ~inc_scale@P),
                        max = 14),
    onset_time = ~infec_time + incub_period
  ) |>
  load_parameter_list(
    inc_shape = 2,
    inc_scale = 1.1,
    append = TRUE
  ) |>
  add_data_wide(dplyr::mutate(housepi::survey, status = NA))

simulated <- hhsimul |>
  simulate(obs_process = function(x) {
    x$onset_date <- ifelse(x$status > 0, floor(x$onset_time + 0.25), NA)
    x$incub_period <- x$onset_time <- x$infec_time <- NULL
    x
  })

hhmodel <- hhsimul |>
  # loading definitions from list drops previous by default
  load_augmented_list(
    infec_time = ~onset_time - incub_period,
    onset_time = list(model = "normal", discr = TRUE,
                      args = c(mean = ~onset_date@D, sd = 1)),
    incub_period = list(model = "gamma", discr = TRUE,
                        args = c(shape = ~inc_shape@P, scale = ~inc_scale@P))
  ) |>
  add_data_wide(simulated)
hhmodel$summarise_cases()

hhmodel$infec_profile <- list(name ="gamma_meanvar", args = c(mean = 3, var = 1))

results <- hhmodel |>
  infer(chains = 3, steps = 10000, adapt = 2000,
        burnin = 1000, thin = 10, progress = TRUE)

# select parameters used for the simulation
simul_parameters <- get_parameters(hhsimul) |>
  dplyr::filter(!const) |>
  dplyr::select(name, value) |>
  dplyr::rename(parameter = name, simulated = value)

# use bayesplot to plot the posterior draws
get_param_draws(results) |>
  bayesplot::mcmc_trace(facet_args = list(ncol = 1)) +
  geom_hline(data = simul_parameters, aes(yintercept = simulated),
             colour = "coral", linewidth = 1)
