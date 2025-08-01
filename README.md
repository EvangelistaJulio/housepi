# housepi 

## Overview

`housepi` is an R package designed for statistical analysis of disease transmission within households. 
The package is a toolkit to perform Bayesian inference and modeling using the implemented template models or user defined models.
By combining ease of use and customizability, housepi supports both off-the-shelf analysis and advanced epidemiological investigation. 

**IMPORTANT:** 
The package has not been released.
Sharing the package or parts of its code without consent of the authors is strictly forbidden.

## Installation

The package is currently under development as a private repository on [gitlab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi).

Prior to installation, check that `devtools` is installed and test that a working C++ toolchain is present by running:

``` r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
pkgbuild::has_build_tools(debug = TRUE)
```

To install `housepi`, open Rstudio, download and unzip the source code from the webpage to your working directory, then run the following command in your R console:

``` r
devtools::install()
```

## Vignettes

In order to build locally your vignettes you will need to install the following packages:

``` r
# packages required for building vignettes
for (pkg in c("ggplot2", "flextable", "bayesplot", "dplyr", "loo")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
```

## License
GNU GENERAL PUBLIC LICENSE, Version 3.

**IMPORTANT:** 
The package has not been released.
Sharing the package or parts of its code without consent of the authors is strictly forbidden.
