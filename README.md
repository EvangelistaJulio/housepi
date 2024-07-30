# housepi <img src="man/figures/logo.png" align="right" height="130"/>

## Description

`housepi` is an R package designed for the statistical analysis and modeling of epidemiological data from household outbreaks. By combining ease of use and customizability, housepi supports both straightforward analysis and advanced epidemiological investigation of disease transmission within households. The package is a toolkit to perform Bayesian inference and modeling using the implemented classical models or user defined models.

### Key Features

-   **Ease of Use**: Load data, model configuration, and run inference with straightforward commands.
-   **Customizability**: Create custom models starting from scratch or modifying default templates.
-   **Model Validation**: The support for both inference and simulations estimation and model validation.
-   **Rcpp Integration**: The package uses Rcpp for efficient computation, making it suitable for handling large and complex datasets.

## Installation

To install `housepi`, use the following command in R:

```{r}
# Missing a non private repository.
```

## Getting started

Here is a quick guide to get started. For a more comprehensive tutorial, please refer to the package vignettes and help files.

Start by loading the library:

```{r}
library(housepi)
```

### Generate a Model

To generate default epidemiological model (e. g. chain binomial model), run the following:

```{r}
hemi <- generate("chain_binom")
```

Alternatively, to customize an epidemiological model from a default one:

```{r}
# Define your model path
my.model.path <- "path/to/your_model"

# Use the template of heterogeneous chain binomial
skeleton("chain_binom_heter", where=my.model.path)
```

Edit the implementation "model.cpp" and/or "model.h" and generate:

```{r}
hemi <- generate(my.model.path)
```

### Build the Model

At any moment you can load the necessary definitions (e. g., database, parameters, setup) and build the model:

```{r}
build(hemi, 
      data = "data.csv", 
      params = "params.csv", 
      setup = "setup.csv")
```

The definitions can be loaded from both a data frame or a file.

### Infer the Model

The package uses the Metropolis-Hastings algorithm to infer the model parameters. Run the Monte Carlo Markov Chain with the following command:

```{r}
mcmc <- infer(hemi)
```

This generate a `coda` mcmc object to be used for extract the inference results (see [documentation](https://rdrr.io/cran/coda/man/mcmc.html)).

### Simulate Data from the Model

Produce a simulated dataset according to the model by running:

```{r}
df <- simulate(hemi)
```

## License

[[MISSING]]
