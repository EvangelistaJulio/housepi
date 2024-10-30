# housepi 

## Description

`housepi` is an R package designed for statistical analysis of disease transmission within households. 
The package is a toolkit to perform Bayesian inference and modeling using the implemented template models or user defined models.
By combining ease of use and customizability, housepi supports both off-the-shelf analysis and advanced epidemiological investigation. 

## Installation

To install `housepi`, use the following command in R:

```{r}
# Missing a non private repository.
```

## Getting started

Here is a just a quick guide, for a more comprehensive tutorial refer to the package vignettes and help files.

Start by loading the library:

```{r}
library(housepi)
```

### Use a Template Model

For a list of available template models, run 

```{r}
# get template list
hemi <- household_model_template()

# get a template description, e. g. chain binomial model
hemi <- household_model_template("chain_binom")
```

To generate an instance of a template model class, run the following:

```{r}
hemi <- HouseholdModelClass("chain_binom")
```

Alternatively, you can provided a custom model (see documentation).

### Associate parameters

Associate the parameters, e. g. from the template definitions

```{r}
hemi$param.template()
```

Alternatively, you can provided custom parameter definitions (see documentation).

### Infer the Model

The package uses the Metropolis-Hastings algorithm to infer the model parameters. 
Run the Monte Carlo Markov Chain with the following command:

```{r}
hemi$infer("data.csv", num_iter = 10000)
```

Extract the inference results from the traces:

```{r}
hemi$traces.summary()
```
