library(testthat)

test_that("Model class inconsistent calls", {
  # blank household model object
  hhmodel <- HouseholdModel$new("blank")

  # traces
  expect_error(hhmodel$traces)

  # household traces
  expect_error(hhmodel$loo())
  expect_error(hhmodel$waic())

  # inference without definitions
  expect_error(hhmodel$infer())
})

test_that("Model class blank calls", {

  expect_warning(HouseholdModel$new())

  # blank household model object
  hhmodel <- HouseholdModel$new("blank")

  # blank default dataset
  dataset <- housepi:::template_blank_dataset_()
  expect_true(all(hhmodel$dataset == dataset))

  # correct datalist call
  datalist <- data.frame(id = dataset$id[1],
                         time = seq(3),
                         measure = seq(3)**3)
  hhmodel$datalist <- datalist
  expect_true(all(hhmodel$datalist == datalist))
  
  # factorization 
  hhmodel$dataset <- housepi::survey
  expect_true(all(hhmodel$dataset == housepi::survey))

  # templates
  expect_warning(hhmodel$load_param_template())
  expect_warning(hhmodel$load_augm_template())

  # zero loglikelihood
  expect_warning(x <- hhmodel$loglikelihood())
  expect_equal(x, 0)

  # blank simulation
  expect_error(hhmodel$simulate())
})

test_that("Chain binomial models are coherent", {
  # loglikelihood from chain binomial model
  hhmodel1 <- HouseholdModel$new("chain-binomial")

  hhmodel1$load_param_list(alpha = 0.15, beta = 0.22)
  ll1 <- hhmodel1$loglikelihood(housepi::survey)

  # loglikelihood from heterogeneous chain binomial model
  # with homogenous rates
  hhmodel2 <- HouseholdModel$new("chain-binomial-heterogeneous")

  hhmodel2$load_param_list(alpha = 0.15, beta = 0.22)
  ll2 <- hhmodel2$loglikelihood(housepi::survey)

  expect_equal(ll1, ll2)
})

test_that("Loglikelihood traces and pointwise traces are coherent", {
  
  hhmodel <- HouseholdModel$new("augm-infection-time")
  hhmodel$load_param_template()
  hhmodel$load_augm_template()
  
  hhmodel$load_group_rate_list(infectivity = "sex", 
                               susceptibility = "sex")
  hhmodel$infer(housepi::survey, sample_size = 10, thin = 5,
                n_adapt = 50, n_burnin = 50)
  total <- hhmodel$traces$loglikelihood
  poitwise <- hhmodel$house_loglike_traces |> 
    posterior::as_draws_matrix() |> rowSums() |> unname()
  expect_equal(total, poitwise)
})

test_that("Augmentation example on incubation from delays is coherent", {
  # loglikelihood from simulated delays
  shape <- 2.0
  scale <- 1.1
  dataset <- data.frame(status = 1,
                        ip = rgamma(223, shape = shape, scale = scale))
  exact <- sum(log(stats::dgamma(dataset$ip, shape = shape, scale = scale)))

  # generate an instance of the household model class
  hhmodel <- HouseholdModel$new("augm-example")
  hhmodel$load_augm_list(
    proba = list(func = "gamma",
                 func_args = c(x = "ip@D", shape = shape, scale = scale)),
    log_proba = "log(proba)"
    )

  # warning for no structure
  expect_warning(x <- hhmodel$loglikelihood(dataset))
  expect_equal(x, exact)
})

test_that("Arithmetic function with undefined x return NaN", {
  hhmodel <- HouseholdModel$new("augm-example")

  # NA values
  hhmodel$dataset <- data.frame(house_id = seq(3), id = seq(3),
                                status = 1, test = c(5.5, NA, -1.5))
  # NaN + something = NaN
  hhmodel$load_augm_list(noise = list(model = "norm"),
                         log_proba = "test@D + noise")
  # NA loglikelihood
  expect_error(hhmodel$infer())

  # log(NaN) = NaN
  hhmodel$load_augm_list(log_proba = "log(test@D)")
  expect_error(hhmodel$infer())
})

test_that("Augmented function example returns expected distribution", {
  hhmodel <- housepi::HouseholdModel$new("blank")

  # define Box-Muller transformation
  hhmodel$load_augm_list(
    norm = list(model = "norm"),
    u1 = list(model = "unif"),
    u2 = list(model = "unif"),
    log_u1 = "log(u1)",
    two_log_u1 = "2 * log_u1",
    m_two_log_u1 = "0 - two_log_u1",
    sqrt_m_two_log_u1 = "sqrt(m_two_log_u1)",
    two_pi_u2 = "6.283185 * u2",
    cos_two_pi_u2 = "cos(two_pi_u2)",
    z0 = "sqrt_m_two_log_u1 * cos_two_pi_u2"
  )

  hhmodel$infer(n_steps = 5e5, progress = 0)
  sample_norm <- hhmodel$traces$norm
  sample_boxmul <- hhmodel$traces$z0

  expect_gt(stats::ks.test(sample_norm, sample_boxmul)$p.value, 0.05)
  expect_gt(stats::t.test(sample_norm, sample_boxmul)$p.value, 0.05)
  expect_gt(stats::var.test(sample_norm, sample_boxmul)$p.value, 0.05)
})

test_that("Augmented infection times model has correct loglikelihood", {
  
  infec_profile_inst <- function(x, param) {
    dgamma(x, shape = param$shape, scale = param$scale)
  }
  
  infec_profile_cumul <- function(x, param) {
    pgamma(x, shape = param$shape, scale = param$scale)
  }
  
  transmission_rate <- function(house, param) {
    param$beta * (param$size_ref / dim(house)[1]) ^ param$gamma
  }
  
  loglike_healthy <- function(house, param) {
    i0 <- which(house$infec_time == min(house$infec_time, na.rm = TRUE))
    member_0 <- house[i0, ]
    
    alpha_cumul <- param$alpha * (param$study_period + member_0$incub_period)
    beta <- transmission_rate(house, param)
    
    output <- 0
    for (i in which(house$status == 0)) {
      roi_cumul <- 0
      for (k in which(house$status == 1)) {
        t_infec_k <- house[k, "infec_time"]
        exposure <- param$study_period + member_0$onset_time - t_infec_k
        roi_cumul <- roi_cumul + infec_profile_cumul(exposure, param)
      }
      output <- output - alpha_cumul - beta * roi_cumul
    }
    
    return(output)
  }
  
  loglike_infectious <- function(house, param) {
    output <- 0.
    i0 <- which(house$infec_time == min(house$infec_time, na.rm = TRUE))
    member_0 <- house[i0, ]
    t_intro <- member_0$infec_time
    
    alpha <- param$alpha
    beta <- transmission_rate(house, param)
    
    for (i in which(house$status == 1)) {
      if (i != i0) {
        t_infec_i <- house[i, "infec_time"]
        roi_inst_i <- 0.
        roi_cumul_i <- 0.
        for (k in which(house$status == 1)) {
          t_infec_k <- house[k, "infec_time"]
          gentime <- t_infec_i - t_infec_k
          if (gentime > 0) {
            roi_inst_i <- roi_inst_i + infec_profile_inst(gentime, param)
            roi_cumul_i <- roi_cumul_i + infec_profile_cumul(gentime, param)
          }
        }
        
        roi_inst_i <- roi_inst_i * beta + alpha
        roi_cumul_i <- roi_cumul_i * beta + alpha * (t_infec_i - t_intro)
        output <- output + log(roi_inst_i) - roi_cumul_i
      }
    }
    
    return(output)
  }
  
  house <- data.frame(id = sprintf("hh_01_%s", seq(1, 7)),
                      house_id = "hh_01",
                      status = c(1, 1, 1, 0, 0, 0, 0),
                      infec_time = c(0.0, 2.7, 3.6, NA, NA, NA, NA),
                      incub_period = c(2.5, 3.1, 2.0, NA, NA, NA, NA))
  house$onset_time <- house$infec_time + house$incub_period
  
  param <- list(alpha = 1e-04, beta = 0.07, gamma = 0.9,
                size_ref = 4, study_period = 45, shape = 2, scale = 1.1)
  
  hhmodel <- HouseholdModel$new("augm-infection-time")
  hhmodel$dataset <- house
  hhmodel$param_defs <- param
  hhmodel$augm_defs <- list(infec_time = "infec_time@D",
                            incub_period = "incub_period@D",
                            onset_time = "onset_time@D")
  hhmodel$add_infec_profile(name = "dgamma", 
                            args = c(shape = "shape@P", scale = "scale@P"))
  
  expect_equal(hhmodel$loglikelihood(), 
               loglike_healthy(house, param) + loglike_infectious(house, param))
})
