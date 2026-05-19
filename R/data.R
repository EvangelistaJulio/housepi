## housepi - Epidemiological Analysis of Household Outbreaks.
##
## Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
## Distributed under the terms of the GNU General Public License (GPLv3).
## See the LICENSE file for details.
##

#' Minimal 'blank' household survey for testing
#'
#' Household data with a single household of 2 individuals.
#'
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{house_id}{Household identifier}
#'   \item{status}{Serological test result}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name blank_house
#' @usage data(blank_house)
#' @format A data frame with 2 rows and 2 columns
NULL

#' Reported "Spanish" flu cases in Frederick (USA), 1918
#'
#' Reported cases of "Spanish" flu in Frederick (USA), in 1918.
#'
#' \describe{
#'   \item{house_id}{Household identifier}
#'   \item{status}{Reported case}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name flu_1918_frederick
#' @usage data(flu_1918_frederick)
#' @format A data frame with 2,207 rows and 2 columns
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_1918_baltimora_frederick.R", package = "housepi")}
#'
#' @references
#' \insertRef{fraser_influenza_2011}{housepi}
NULL

#' Reported "Spanish" flu cases in Baltimora (USA), 1918
#'
#' Reported cases of "Spanish" flu in Baltimora (USA), in 1918.
#'
#' \describe{
#'   \item{house_id}{Household identifier}
#'   \item{status}{Reported case}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name flu_1918_baltimora
#' @usage data(flu_1918_baltimora)
#' @format A data frame with 28,977 rows and 2 columns
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_1918_baltimora_frederick.R", package = "housepi")}
#'
#' @references
#' \insertRef{fraser_influenza_2011}{housepi}
NULL

#' SARS-CoV-2 2021 household survey from the Israeli COVID-19 Family Study
#'
#' Prospective cohort conducted during COVID-19 Delta variant wave, 2021.
#' Follow-up started with first PCR positive test, symptoms onset were reported.
#' Contacts averaged 3 tests during follow-up.
#'
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{house_id}{Household identifier}
#'   \item{status}{PCR test result}
#'   \item{test_date}{Date of first positive test}
#'   \item{vacc_status}{Vaccinated before survey}
#'   \item{age}{Age}
#'   \item{isolation}{Whether isolated}
#'   \item{sex}{Sex}
#'   \item{asymptomatic}{Whether asymptomatic}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sarscov2_2021_icofs
#' @usage data(sarscov2_2021_icofs)
#' @format A data frame with 4,168 rows and 11 columns
#' @source <https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt>
#'
#' @details
#' Original data from household transmission study conducted in Israel, 2021.
#'
#' Processing steps:
#' \enumerate{
#'   \item Isolation status converted from "no quarantine" string to logical
#'   \item Positive test status derived from presence of first_pos date
#'   \item Test dates converted to days since 2021-01-01
#'   \item Selected relevant columns and renamed for consistency
#' }
#'
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_sarscov2_2021_icofs.R", package = "housepi")}
#'
#' @references
#' \insertRef{cortier_factors_2025}{housepi}
NULL

#' Flu 2014-15 household survey from Matsumoto, Japan
#'
#' Household survey in Matsumoto city for parents having children in primary school.
#' Report on influenza episodes of family members during the season 2014-15.
#' Original observational study from Uchida et al. (2017).
#'
#' \describe{
#'   \item{house_id}{Household identifier}
#'   \item{member}{Family member}
#'   \item{status}{Symptoms reported}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name flu_2014_matsumoto
#' @usage data(flu_2014_matsumoto)
#' @format A data frame with 47,213 rows and 3 columns
#' @source <https://raw.githubusercontent.com/akira-endo/HHstudy_FluMatsumoto2014-15/master/HHdata_FluMatsumoto2014-15.csv>
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_2014_matsumoto.R", package = "housepi")}
#'
#' @references
#' \insertRef{uchida_effectiveness_2017}{housepi}
#'
#' \insertRef{endo_fine_scale_2019}{housepi}
NULL
