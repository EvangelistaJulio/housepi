# housepi

**housepi** is an R package designed for the statistical analysis and modeling of epidemiological data from household outbreaks.
By combining ease-of-use with customizability, housepi supports both novice users and advanced researchers in their epidemiological investigations of disease transmission within households.
The package is a toolkit to perform Bayesian inference and modeling using the implemented classical models or user defined models.

## Key Features

- **Ease of Use**: Load data, model configuration, and run inference with straightforward commands.
- **Customizability**: Create custom models starting from scratch or modifying default templates.
- **Model Validation**: The support for both inference and simulations estimation and model validation.
- **Rcpp Integration**: The package uses Rcpp for efficient computation, making it suitable for handling large and complex datasets.

## Installation

To install **housepi**, use the following command in R:

```{r}
# Missing a non private repository.
```

## Getting started

Here is a quick guide to get started. For a more comprehensive tutorial, please refer to the package vignettes and help files.

Start by loading the library:

```{r}
library(housepi)
```

### Generating a Model

To generate default epidemiological model (e. g. chain binomial model), run the following:

```{r}
hemi <- generate("chain_binom")
```

Alternatively, to customize an epidemiological model from a default one:

```{r}
# Define your model path
model_path <- "path/to/your_model"

# Use the template of heterogeneous chain binomial
skeleton("chain_binom_heter", where=model_path)
```

Edit the implementation "model.cpp" and/or "model.h" and generate:

```{r}
hemi <- generate(model_path)
```

### Building the Model

At any moment you can load the necessary definitions (e.g., database, parameters, setup) and build the model:

```{r}
build(hemi, 
      data = "data.csv", 
      params = "params.csv", 
      setup = "setup.csv")
```

The definitions can be loaded from both a data frame or a file.

### Running Inferece

The package uses the Metropolis-Hastings algorithm to infer the model parameters. Run the Monte Carlo Markov Chain with the following command:

```{r}
mcmc <- infer(hemi)
```

This generate a coda mcmc object to be used for extract the inference results (see [documentation](https://rdrr.io/cran/coda/man/mcmc.html)).

### Simulating Data

Produce a simulated dataset according to the model by runnning:

```{r}
df <- simulate(hemi)
```

## License

[[MISSING]]