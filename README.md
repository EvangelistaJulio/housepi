# housepi <img src="man/figures/logo.png" align="right" height="130"/>

## Description

`housepi` is an R package designed for the statistical analysis and modeling of epidemiological data from household outbreaks. By combining ease of use and customizability, housepi supports both straightforward analysis and advanced epidemiological investigation of disease transmission within households. The package is a toolkit to perform Bayesian inference and modeling using the implemented classical models or user defined models.

### Key Features

-   **Ease of Use**: Load data, model configuration, and run inference with straightforward commands.
-   **Customizability**: Create custom models starting from scratch or modifying default templates.
-   **Model Validation**: Support for inference, simulations and model validation.
-   **Rcpp Integration**: Use of Rcpp for efficient computation, making it suitable for handling large and complex datasets.

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
hemi <- generate_model("chain_binom")
```

Alternatively, you can provided a custom model (see documentation).

### Build the Model

At any moment you can load the necessary definitions to build the model (e. g., database, parameters, setup):

```{r}
build_model(hemi,
            data = "data.csv", 
            params = "params.csv", 
            setup = "setup.csv")
```

The definitions can be loaded from both a file or a dataframe.

### Infer the Model

The package uses the Metropolis-Hastings algorithm to infer the model parameters. Run the Monte Carlo Markov Chain with the following command:

```{r}
mcmc <- run_inference(hemi)
```

This generate a `coda` mcmc object to be used for extract the inference results (see [documentation](https://rdrr.io/cran/coda/man/mcmc.html)).

### Simulate Data from the Model

Produce a simulated dataset according to the model by running:

```{r}
hemi$simulate()
df <- hemi$data()
```

## License

[[MISSING]]
