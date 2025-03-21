# housepi 

## Description

`housepi` is an R package designed for statistical analysis of disease transmission within households. 
The package is a toolkit to perform Bayesian inference and modeling using the implemented template models or user defined models.
By combining ease of use and customizability, housepi supports both off-the-shelf analysis and advanced epidemiological investigation. 

## Installation

The package is currently under development as a private repository on [gitlab](https://gitlab.pasteur.fr/mmmi-pasteur/housepi).

Prior to installation, test that a working C++ toolchain is present onyour machine:

```{r}
pkgbuild::has_build_tools(debug = TRUE)
```

To install `housepi`, download the source as a zip file from the webpage to your curernt working directory, then run the following command in your R console:

```{r}
install.packages("housepi-main.zip", repos = NULL, type = "source")
```

## License

**IMPORTANT:** 
The package has not been released.
Sharing the package or parts of its code without consent of the authors is strictly forbidden.