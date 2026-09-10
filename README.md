# housepi <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![release](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/badges/release.svg?key_text=release&key_width=50)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/releases)
[![pipeline status](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/pipeline.svg)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![coverage report](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/coverage.svg?job=coverage)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mmmi-pasteur.pages.pasteur.fr/housepi)
<!-- badges: end -->

housepi is an open source R package for the statistical analysis
and modelling of epidemiological data from household outbreaks.
Household transmission surveys are a standard route to estimating
key epidemiological parameters.
The relative ease of sampling household data makes it a readily available
source of information during the early phases of an epidemic
(e. g. longitudinal studies).
Additionally, retrospective serological surveys can provide critical insights
into the history of transmission within a population (final size studies).

> [!IMPORTANT]
> The package has not been released yet. Certain options are still under development.

## Installation

Prior to installation, check that `devtools` is installed and test that a working C++ toolchain is present by running:
``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
pkgbuild::has_build_tools(debug = TRUE)
```

Download the source code on [GitLab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi), set the main folder as your working directory and run
```r
devtools::install()
```

## License

GPLv3 © Institut Pasteur, Paris (France).
Distributed under the terms of the GNU General Public License, Version 3.
For details, see [LICENSE](LICENSE.md).
