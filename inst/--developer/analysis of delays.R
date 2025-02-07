# fix random seed
set.seed(1)

# time 1
value_duration_1 <- 1:4
proba_duration_1 <- c(0.28, 0.34, 0.26, 0.12)

# time 2
proba_exact_duration2 <- 0.73

# incubation period
delay_value <- function() {
  return(rgamma(1, shape = 2.0, scale = 1.1))
}

simulate_delay_data <- function(x) {
  # exposure interval
  t1_min <- 0
  t1_max <- value_duration_1[rmultinom(1, 1, proba_duration_1)[ ,1] == 1]

  # exact time of infection
  t1 <- runif(1, t1_min, t1_max)  
  # onset time
  t2 <- t1 + delay_value()

  # onset interval
  t2_min <- as.integer(t2)
  t2_max <- t2_min + 1.
  
  # reporting errors on the interval of t2
  if(runif(1) > proba_exact_duration2){
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

full_dataset <- as.data.frame(t(sapply(1:223, simulate_delay_data)))

#-------checks--------
# mean(full_dataset$t2 - full_dataset$t1)
# sd(full_dataset$t2 - full_dataset$t1)
# 
# av_t1 <- 0.5 * (full_dataset$t1_min + full_dataset$t1_max)
# av_t2 <-  0.5 * (full_dataset$t2_min + full_dataset$t2_max)
# mean(av_t2 - av_t1)

#-----------final dataset----------
delays <- data.frame(exposure_min = full_dataset$t1_min,
                     exposure_max = full_dataset$t1_max - 1, # censoring
                     onset_min = full_dataset$t2_min,
                     onset_max = full_dataset$t2_max - 1) # censoring
