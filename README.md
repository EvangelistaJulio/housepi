# housepi <img src="man/figures/logo.png" align="right" alt="" width="120" />

We introduce "housepi", an open source R package for the statistical analysis
and modeling of epidemiological data from household outbreaks.
Mathematical modeling has proven to be a valuable tool for estimating
key epidemiological parameters from household transmission surveys.
The relative ease of sampling household data makes it a readily available
source of information during the early phases of an epidemic
(e. g. longitudinal studies).
Additionally, retrospective serological surveys can provide critical insights
into the history of transmission within a population (finite size studies).
The package aims to combine ease of use with flexibility, supporting both
straightforward analyses and advanced epidemiological investigations.

::: caution
The package is still under development.
:::

## Installation

Prior to installation, check that `devtools` is installed and test that a working C++ toolchain is present by running:
``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
pkgbuild::has_build_tools(debug = TRUE)
```

::: note
The package is currently under development as a private repository on GitLab.
You need access permissions and a valid GitLab token to install it.
:::

If you have SSH access set up, just run the following
``` r
install.packages("remotes")
remotes::install_git("git@gitlab.pasteur.fr:mmmi-pasteur/housepi.git", ref = "main")
```
Otherwise download the source code on [GitLab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi), set the main folder as your working directory and run
```r
devtools::install()
```

## License

housepi is distributed under the terms of the GNU General Public License, Version 3 (GPLv3).
See the LICENSE file for details.

::: warning
The package has not been released. Sharing the package or parts of its code without consent of the authors is strictly forbidden.
:::
