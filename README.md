# housepi <img src="man/figures/logo.svg" align="right" alt="" width="120" />

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![release](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/badges/release.svg?key_text=release&key_width=50)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/releases)
[![pipeline status](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/pipeline.svg)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![coverage report](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/badges/main/coverage.svg?job=coverage)](https://gitlab.pasteur.fr/mmmi-pasteur/housepi/-/commits/main)
[![license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://mmmi-pasteur.pages.pasteur.fr/housepi)
<!-- badges: end -->

We introduce "housepi", an open source R package for the statistical analysis
and modelling of epidemiological data from household outbreaks.
Mathematical modelling has proven to be a valuable tool for estimating
key epidemiological parameters from household transmission surveys.
The relative ease of sampling household data makes it a readily available
source of information during the early phases of an epidemic
(e. g. longitudinal studies).
Additionally, retrospective serological surveys can provide critical insights
into the history of transmission within a population (final size studies).
The package aims to combine ease of use with flexibility, supporting both
straightforward analyses and advanced epidemiological investigations.

> [!IMPORTANT]
> The package is still under development and certain options might need further testing.

## Installation

Prior to installation, check that `devtools` is installed and test that a working C++ toolchain is present by running:
``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
pkgbuild::has_build_tools(debug = TRUE)
```

> [!NOTE]
> The package is currently under development as a private repository on GitLab.
> You need access permissions and a valid GitLab token to install it.

Download the source code on [GitLab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi), set the main folder as your working directory and run
```r
devtools::install()
```

## License

GPLv3 © Institut Pasteur, Paris (France).
Distributed under the terms of the GNU General Public License, Version 3.
For details, see [LICENSE](LICENSE.md).

> [!CAUTION]
> The package has not been released. Sharing the package or parts of its code without consent of the authors is strictly forbidden.
