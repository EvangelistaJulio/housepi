# housepi (development version)

* Articles, README and help pages are cleaned up, and the custom model example is a community-only final-size model.
* `measure_geom_mean_hh_size()` is removed.
* The upper tail ESS in `summarise_draws()` and in the convergence diagnostics counts draws at or above
  the 95% quantile, so a parameter with 5% or more of its draws on the top of its support no longer gets
  an `NA` tail ESS and a convergence warning.
* `summarise_draws()` with `.args` and its default summaries no longer fails with an unused-argument error.

# housepi 0.11.0

Pre-release.

* Household transmission models are declared with `new_housepi_model()`, from the built-in
  `final_size` and `case_ascertained` models or from custom C++ and discrete-time models
  compiled at runtime.
* Household data are attached as a dataset, one row per member (`add_dataset()`).
* Parameters take priors and proposals written in a small DSL (`add_parameter()`).
* Group and contact rates (`add_group_rate()`, `add_contact_rate()`) and infectivity
  profiles (`add_infectivity_profile()`) shape transmission between members.
* Unobserved quantities are handled by data augmentation: `add_latent()`, `add_observation()`
  and `add_status_process()`.
* `run_inference()` fits the model with Bayesian MCMC over several chains, run at the same time when
  `threads` allows. `run_simulation()` draws outbreaks from it and `evaluate_log_lik()` scores a parameter set.
* Draws are read with the `get_*_draws()` family and summarised with `summarise_draws()`.
* Fits are validated with predictive checks (`check_predictive()`) and simulation-based calibration (`check_calibration()`), against the prior or posterior.
  Performance at known parameter values is evaluated with `check_performance()`.
