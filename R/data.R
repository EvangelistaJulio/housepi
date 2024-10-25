#' Household Survey 
#'
#' Household data from a simulated demographic structure with a simulated outbreak.
#' We assume an infectious and preventable disease. 
#'
#' @format ## `survey`
#' A data frame with 2,249 rows and 6 columns:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{house_id}{Household identifier}
#'   \item{status}{Serological test result}
#'   \item{age}{Age}
#'   \item{sex}{Sex: `1` women, `2` men}
#'   \item{vacc_status}{Vaccinated before survey}
#' }
"survey"

#' Individual Delays
#'
#' Reported individual periods of exposure and symptom onset.
#' Simulated assuming an incubation period gamma distributed.
#' Date are shifted so that exposure period starts from day 0.
#'
#' @format ## `delays`
#' A data frame with 223 rows and 4 columns:
#' \describe{
#'   \item{exposure_min}{Day exposure starts}
#'   \item{exposure_max}{Day exposure ends}
#'   \item{onset_min}{Day symptom onset starts}
#'   \item{onset_max}{Day symptom onset ends}
#' }
"delays"