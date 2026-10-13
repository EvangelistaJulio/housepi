# housepi <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![release](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/badges/release.svg?key_text=release&key_width=50)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/releases)
[![pipeline status](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/pipeline.svg)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![coverage report](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/coverage.svg?job=coverage)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![license](https://img.shields.io/badge/license-GPL--3.0--or--later-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mmmi-pasteur.pages.pasteur.fr/housepi)
<!-- badges: end -->

*housepi* is an R package for the statistical analysis and modelling
of infectious disease transmission within households.
Household transmission models are specified declaratively, fitted to data with Bayesian MCMC
and used to simulate outbreaks.
They support data augmentation, infectivity profiles, and individual group and contact risks.
Validation is supported through simulations studies and posterior predictive checks.

> [!IMPORTANT]
> *housepi* has not been released yet: the API may change between versions until the first release.

## Installation

*housepi* is compiled from source, so it needs a working C++17 toolchain.
The same toolchain compiles custom and discrete-time models at runtime.
Check that it is available with:
``` r
install.packages(c("pkgbuild", "remotes"))
pkgbuild::has_build_tools(debug = TRUE)
```
Then install the development version from [GitLab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi):
``` r
remotes::install_gitlab("mmmi-pasteur/housepi", host = "gitlab.pasteur.fr")
```

## Getting started

Fit the final-size model to the 1918 influenza outbreak in Frederick (USA),
shipped with the package, using the parameter definitions of the model template:
``` r
library(housepi)

fit <- flu_1918_frederick |>
  new_housepi_model("final_size") |>
  add_parameter(import = "template") |>
  run_inference()

summarise_draws(fit)
```
Learn more about data format, parameter definitions, inference and simulation,
in the guide [Get started](https://mmmi-pasteur.pages.pasteur.fr/housepi/articles/housepi).

## Citation

To cite *housepi* in publications, run `citation("housepi")`.

## License

Copyright 2024-2026 Institut Pasteur (Paris), France.
Distributed under the terms of the GNU General Public License, version 3 or later.
For details, see [LICENSE](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/blob/main/LICENSE.md).
