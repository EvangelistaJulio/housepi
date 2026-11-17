# housepi <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![release](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/badges/release.svg?key_text=release&key_width=50)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/releases)
[![pipeline status](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/pipeline.svg)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![coverage report](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/coverage.svg?job=coverage)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![license](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mmmi-pasteur.pages.pasteur.fr/housepi)
[![HAL](https://img.shields.io/badge/HAL-pasteur--05869643-blue.svg)](https://pasteur.hal.science/pasteur-05869643)
<!-- badges: end -->

housepi is an R package for the statistical analysis and modelling
of infectious disease transmission within households.
Household transmission models are specified declaratively,
fitted to data with Bayesian MCMC and used to simulate outbreaks.
The package support data augmentation, infectivity profiles, and individual group and contact rates.
At the same time it support standard Bayesian workflow through predictive checks and simulation studies.

## Installation

housepi needs a working C++17 toolchain, check that it is available with:
``` r
install.packages(c("pkgbuild", "remotes"))
pkgbuild::has_build_tools(debug = TRUE)
```
The same toolchain compiles custom and discrete-time models at runtime.

Install it from [GitLab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi):
``` r
remotes::install_gitlab("mmmi-pasteur/housepi", host = "gitlab.pasteur.fr")
```

## Getting started

Start from a public dataset shipped with the package 
(e.g ["Spanish" flu in Frederick, USA, 1918](https://mmmi-pasteur.pages.pasteur.fr/housepi/reference/flu_1918_frederick.html)),
simulate an outbreak using a templated model 
(e.g. [final-size model](https://mmmi-pasteur.pages.pasteur.fr/housepi/reference/template_final_size)),
associate standard parameter definitions,
then fit the model to the simulated study:
``` r
library(housepi)

# fix the seed to ensure reproducibility on the current machine
set.seed(42)

# define household transmission model from a template or a custom script
hh_model <- new_housepi_model(flu_1918_frederick, "final-size")

# associate custom parameter definitions or use templated choices
hh_model <- hh_model |>
  add_parameter(import = "template")

# simulate the outbreak at the given parameters
simulated_data <- hh_model |>
  run_simulation(c(alpha = 0.1, beta = 0.3, gamma = 1))

# run the inference using multi-threading
hh_fit <- hh_model |>
  add_dataset(simulated_data) |>
  run_inference(steps = 2000, chains = 2, threads = 2)

# extract the results using standard Bayesian packages
hh_fit |>
  posterior::summarise_draws() 
```

Use
[Bayesian data augmentation](https://mmmi-pasteur.pages.pasteur.fr/housepi/articles/data_augmentation)
to deal with unobserved variables, as infection times in 
[case-ascertained studies](https://mmmi-pasteur.pages.pasteur.fr/housepi/reference/template_case_ascertained).

[Get started](https://mmmi-pasteur.pages.pasteur.fr/housepi/articles/housepi)
and learn more about data format, parameter definitions, inference and simulation.

## Citation

> Camaglia F, Cristancho-Fajardo L, Cortier T, Tsang TK, Bosetti P, Brault A, Cowling BJ, Cauchemez S (2026).
> Seamless household epidemiological modelling with housepi.
> HAL preprint [pasteur-05869643](https://pasteur.hal.science/pasteur-05869643).

Here for [BibTeX entry](https://mmmi-pasteur.pages.pasteur.fr/housepi/authors.html#citation).

## Contributing

Report bugs and request features on the [GitLab issue tracker](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/work_items), following the [contributing guide](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/blob/main/CONTRIBUTING.md).

## AI disclosure

Claude models (Anthropic) were used to develop the code and documentation of this package.
All code was validated against standard statistical tests and previous implementations.
The authors take full responsibility for its content.

## License

Copyright 2024-2026 Institut Pasteur (Paris), France.
Distributed under the terms of the GNU General Public License, version 3 or later.
For details, see [LICENSE](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/blob/main/LICENSE.md).
