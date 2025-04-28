
source("R/utils.R")
source("R/template.R")
source("R/augmented.R")
source("R/parameter.R")
source("R/infectivity.R")
source("R/group_rate.R")
source("R/expected_cases.R")
source("R/model_generator.R")
source("R/model_wrapper.R")
source("R/model_class.R")
source("R/sampler_class.R")

hhmodel <- HouseholdModel$new("chain-binomial")
hhmodel$load_param_template()
hhmodel$infer(
  dataset = housepi::survey, # import the dataset as a data.frame
  n_steps = 20000,  # number of mcmc steps (per chain)
  n_burnin = 1000,  # steps to burn
  thin = 8,         # delay between sampling steps
  n_adapt = 1000,   # finite adaptation steps before sampling
  n_chains = 4,     # number of independent chains
  seed = 6043,      # internal random number generator seed
  n_threads = 1     # number of threads
)

#FIXME: why is this?


loo::waic(simplify2array(lapply(hhmodel$traces, as.matrix)))
x <- simplify2array(lapply(hhmodel$traces, as.matrix))
typeof(x)

library(posterior)

# Imagine you have 3 chains of 1000 iterations each for 2 parameters
chain1 <- data.frame(a = rnorm(10), b = rnorm(10))
chain2 <- data.frame(a = rnorm(10), b = rnorm(10))
chain3 <- data.frame(a = rnorm(10), b = rnorm(10))

# Turn into a draws_array
draws <- as_draws_array(list(chain1, chain2, chain3))

# Now you can easily use it:
posterior::rhat(draws)

bayesplot::mcmc_trace(draws)  # traceplot
loo::loo(array(draws, dim = c(3000, 2)))  # careful: loo expects matrix
