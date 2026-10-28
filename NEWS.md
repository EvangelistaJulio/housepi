# housepi 0.10.0

First tagged version.

* Household transmission models are declared with `new_housepi_model()`, from the built-in
  `final_size` and `case_ascertained` models or from custom C++ and discrete-time models
  compiled at runtime.
* Household data are imported in wide or long format (`add_data_wide()`, `add_data_long()`).
* Parameters take priors and proposals written in a small DSL (`add_parameter()`).
* Group and contact rates (`add_group_rate()`, `add_contact_rate()`) and infectivity
  profiles (`add_infectivity_profile()`) shape transmission between members.
* Unobserved quantities are handled by data augmentation: `add_latent()`, `add_observation()`,
  `add_status_process()` and `add_reporting_process()`.
* `run_inference()` fits the model with Bayesian MCMC over several chains,
  `run_simulation()` draws outbreaks from it and `evaluate_log_lik()` scores a parameter set.
* Draws are read with the `get_*_draws()` family and summarised with `summarise_draws()`.
* Fits are validated with posterior predictive checks (`check_predictive()`), parameter
  recovery checks (`check_recovery()`) and simulation-based calibration (`check_calibration()`).
