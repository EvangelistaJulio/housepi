
#' @title
#' Generate Individual Delays
#' @description
#' Reported individual periods of exposure and symptom onset.
#' Simulated assuming an incubation period gamma distributed.
#' Date are assumed shifted so that exposure period starts from day 0.
#'
#' @return `data.frame`\cr
#' \item{exposure_min}{Day exposure starts}
#' \item{exposure_max}{Day exposure ends}
#' \item{onset_min}{Day symptom onset starts}
#' \item{onset_max}{Day symptom onset ends}
#' \item{status}{Individual status }
#'
#' @keywords internal
generate_delays <- function(n = 223,
                            mean = 2.2,
                            var = 2.5,
                            report_error = 0.73) {
  # time 1
  value_duration_1 <- 1:4
  proba_duration_1 <- c(0.28, 0.34, 0.26, 0.12)

  shape <- mean ^2 / var
  scale <- var / mean
  # incubation period
  delay_value <- function() return(rgamma(1, shape = shape, scale = scale))

  simulate_delay_data <- function(x) {
    # exposure interval
    t1_min <- 0
    t1_max <- value_duration_1[rmultinom(1, 1, proba_duration_1)[, 1] == 1]

    # exact time of infection
    t1 <- runif(1, t1_min, t1_max)
    # onset time
    t2 <- t1 + delay_value()

    # onset interval
    t2_min <- as.integer(t2)
    t2_max <- t2_min + 1.

    # reporting errors on the interval of t2
    if (runif(1) > report_error) {
      if (runif(1) < 0.5) {
        t2_min <- t2_min - 1.
      } else {
        t2_max <- t2_max + 1.
      }
    }

    return(c(t1 = t1,
             t2 = t2,
             t1_min = t1_min,
             t1_max = t1_max,
             t2_min = t2_min,
             t2_max = t2_max))
  }

  full_dataset <- as.data.frame(t(sapply(1:n, simulate_delay_data)))

  #-----------final dataset----------
  data.frame(
    status = 1,
    exposure_min = full_dataset$t1_min,
    exposure_max = full_dataset$t1_max - 1,
    onset_min = full_dataset$t2_min,
    onset_max = full_dataset$t2_max - 1
  )
}
