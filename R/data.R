#
#        __
#       / /_  ____  __  __ ____ __  ____  __
#      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
#     / / / / /_/ / /_/ (__  )  __/ /_/ / /
#    /_/ /_/\____/\____/____/\___/ ____/_/
#                               /_/
#
# housepi - Epidemiological Analysis of Household Outbreaks.
#
# Author: Francesco Camaglia
# Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
#
# Distributed under the terms of the GNU General Public License (GPLv3).
# See the LICENSE file for details.
#
# data.R
#

#' Household Testing Survey
#'
#' Household data with a simulated outbreak to be used for housepi vignettes.
#' Note: this will be dropped in future versions of the package.
#'
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{house_id}{Household identifier}
#'   \item{status}{Serological test result}
#'   \item{age}{Age}
#'   \item{sex}{Sex}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name survey
#' @usage data(survey)
#' @format A data frame with 2,249 rows and 6 columns
NULL

#' Flu 1918 Frederick
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
#' @source <https://academic.oup.com/aje/article/174/5/505/117666>
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_1918_baltimora_frederick.R", package = "housepi")}
#'
#' @references
#' Fraser, et al. (2011)
#' Influenza Transmission in Households During the 1918 Pandemic.
#' American Journal of Epidemiology 174, 505-514.
#' [doi:10.1093/aje/kwr122](https://doi.org/10.1093/aje/kwr122)
NULL

#' Flu 1918 Baltimora
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
#' @source <https://academic.oup.com/aje/article/174/5/505/117666>
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_1918_baltimora_frederick.R", package = "housepi")}
#'
#' @references
#' Fraser, et al. (2011)
#' Influenza Transmission in Households During the 1918 Pandemic.
#' American Journal of Epidemiology 174, 505-514.
#' [doi:10.1093/aje/kwr122](https://doi.org/10.1093/aje/kwr122)
NULL

#' SARS-CoV-2 2021 Israeli Covid-19 Family Study
#'
#' Prospective cohort conducted during COVID-19 Delta variant wave in Israel, 2021.
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
#'   \item{isolation}{Wheter isolated}
#'   \item{sex}{Sex}
#'   \item{asymptomatic}{Wheter asymptomatic}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name sarscov2_2021_icofs
#' @usage data(sarscov2_2021_icofs)
#' @format A data frame with 4,168 rows and 9 columns
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
#' Cortier, et al. (2025)
#' Factors Associated With the Transmission of the Delta Severe Acute Respiratory Syndrome Coronavirus 2 Variant in Households:
#' The Israeli COVID-19 Family Study (ICoFS)
#' The Journal of Infectious Diseases; jiaf001.
#' [doi:10.1093/infdis/jiaf001](https://doi.org/10.1093/infdis/jiaf001)
NULL

#' Flu 2014-15 Matsumoto
#'
#' Household survey in Matsumoto city (Japan) for parents having children in primary school.
#' Report on influenza episodes of family members during the season 2014-15.
#' Original observational study from Uchida, et al. (2017).
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
#' @format A data frame with 42,045 rows and 3 columns
#' @source <https://raw.githubusercontent.com/akira-endo/HHstudy_FluMatsumoto2014-15/master/HHdata_FluMatsumoto2014-15.csv>
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_2014_matsumoto.R", package = "housepi")}
#'
#' @references
#' Endo, et al. (2019)
#' Fine-scale family structure shapes influenza transmission risk in households:
#' Insights from primary schools in Matsumoto city, 2014/15.
#' PLOS Computational Biology;15: e1007589.
#' [doi:10.1371/journal.pcbi.1007589](https://doi.org/10.1371/journal.pcbi.1007589)
NULL
