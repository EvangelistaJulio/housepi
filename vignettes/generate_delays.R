
#' @title
#' Generate Individual Delays
#' @description
#' Reported individual periods of exposure and symptom onset.
#' Simulated assuming an incubation period gamma distributed.
#' Date are assumed shifted so that exposure period starts from day 0.
#'
#' @param size `integer(1)` Number of individuals to simulate.
#' @param mean `numeric(1)` Mean of the incubation period distribution.
#' @param df `numeric(1)` Std. deviation of the incubation period distribution.
#' @param report_error `numeric(1)` Probability of reporting error on symptom onset date.
#' @param proba_duration `numeric()` Probabilities of duration of exposure.
#'
#' @return `data.frame`\cr
#' \item{exposure_min}{Day exposure starts}
#' \item{exposure_max}{Day exposure ends}
#' \item{onset_min}{Day symptom onset starts}
#' \item{onset_max}{Day symptom onset ends}
#' \item{status}{Individual status }
#'
#' @keywords internal
generate_delays_report <- function(size,
                                   mean = 2.2,
                                   var = 1.4,
                                   report_error = 0.73,
                                   proba_duration = c(0.28, 0.34, 0.26, 0.12)) {
  shape <- mean * mean / var
  scale <- var / mean

  simulate_delay_data <- function(x) {
    # exposure interval
    t1_min <- 0
    t1_max <- sample(length(proba_duration), 1, prob = proba_duration)

    # exact time of infection
    t1 <- runif(1, t1_min, t1_max)
    # onset time + incubation period
    t2 <- t1 + rgamma(1, shape = shape, scale = scale)

    # onset interval
    t2_min <- as.integer(t2)
    t2_max <- t2_min + 1.

    # reporting errors on the interval of t2
    if (runif(1) > report_error) {
      if (runif(1) < 0.5) t2_min <- t2_min - 1.
      else t2_max <- t2_max + 1.
    }

    c(t1_min = t1_min, t1_max = t1_max,
      t2_min = t2_min, t2_max = t2_max)
  }

  full_dataset <- as.data.frame(t(sapply(1:size, simulate_delay_data)))

  # final dataset
  data.frame(status = 1,
             exposure_min = full_dataset$t1_min,
             exposure_max = full_dataset$t1_max - 1,
             onset_min = full_dataset$t2_min,
             onset_max = full_dataset$t2_max - 1)
}
