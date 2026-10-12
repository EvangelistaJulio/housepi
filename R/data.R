## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Minimal household survey for testing
#'
#' Household data with a single household of 2 individuals.
#'
#' \describe{
#'   \item{id}{Individual identifier;}
#'   \item{hh_id}{Household identifier;}
#'   \item{status}{Serological test result.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @examples
#' dummy_household
#' new_housepi_model(dummy_household, "test")
#' @name dummy_household
#' @usage data(dummy_household)
#' @format A data frame with 2 rows and 3 columns
NULL

#' Reported "Spanish" flu cases in Frederick (USA), 1918
#'
#' Reported cases of "Spanish" flu in Frederick (USA), in 1918.
#'
#' \describe{
#'   \item{hh_id}{Household identifier;}
#'   \item{status}{Reported case.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @examples
#' measure_ar_by_hh_size(flu_1918_frederick)
#' @name flu_1918_frederick
#' @usage data(flu_1918_frederick)
#' @format A data frame with 2,207 rows and 2 columns
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_1918_baltimore_frederick.R", package = "housepi")}
#'
#' @references
#' \insertRef{fraser_influenza_2011}{housepi}
NULL

#' Reported "Spanish" flu cases in Baltimore (USA), 1918
#'
#' Reported cases of "Spanish" flu in Baltimore (USA), in 1918.
#'
#' \describe{
#'   \item{hh_id}{Household identifier;}
#'   \item{status}{Reported case.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @examples
#' measure_ar_by_hh_size(flu_1918_baltimore)
#' @name flu_1918_baltimore
#' @usage data(flu_1918_baltimore)
#' @format A data frame with 28,977 rows and 2 columns
#'
#' @details
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_1918_baltimore_frederick.R", package = "housepi")}
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
#'   \item{id}{Individual identifier;}
#'   \item{hh_id}{Household identifier;}
#'   \item{status}{Follow-up outcome: 0 escaped, 1 symptomatic, 2 asymptomatic;}
#'   \item{age}{Age class, adult or child;}
#'   \item{sex}{Sex;}
#'   \item{start_date}{Follow-up start: household's first positive PCR;}
#'   \item{end_date}{Follow-up end: 21 days after \code{start_date};}
#'   \item{test_date}{Date of the first positive test within the follow-up;}
#'   \item{vacc_status}{Vaccination status at the follow-up start;}
#'   \item{recovered}{Whether previously recovered from infection;}
#'   \item{isolated}{Whether isolated;}
#'   \item{asymptomatic}{Whether asymptomatic.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @examples
#' head(sarscov2_2021_icofs)
#' hh_model <- new_housepi_model(sarscov2_2021_icofs, "case_ascertained")
#' summarise_house_cases(hh_model)
#' @name sarscov2_2021_icofs
#' @usage data(sarscov2_2021_icofs)
#' @format A data frame with 4,168 rows and 12 columns
#' @source <https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt>
#'
#' @details
#' Each household is followed for 21 days from its first positive SARS-CoV-2 PCR test
#' as in the original analysis. Study starts from July 2021.
#' A member first testing positive after \code{end_date} counts as escaped.
#' In 135 households several members test positive on the first day.
#' Vaccination status (\code{vacc_status}) reflects the doses effective
#' at the start of the follow-up, assuming a 7-day delay for effectiveness.
#'
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_sarscov2_2021_icofs.R", package = "housepi")}
#'
#' @references
#' \insertRef{cortier_factors_2025}{housepi}
NULL

#' Flu 2014-15 household survey from Matsumoto, Japan
#'
#' Household survey in Matsumoto city for families with children in primary school.
#' Report on influenza episodes of family members during the season 2014-15.
#' Original observational study from Uchida et al. (2017).
#'
#' \describe{
#'   \item{hh_id}{Household identifier;}
#'   \item{member}{Family member;}
#'   \item{status}{Symptoms reported.}
#' }
#'
#' @docType data
#' @keywords datasets
#' @examples
#' table(flu_2014_matsumoto$member, flu_2014_matsumoto$status)
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
