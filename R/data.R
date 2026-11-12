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

#' SARS-CoV-2 2021 case-ascertained household study from Israel (ICoFS)
#'
#' Households with a SARS-CoV-2 case during the Delta variant wave,
#' recruited around Sheba Medical Center, Israel, in July and August 2021
#' for the Israeli COVID-19 Family Study.
#' Members were tested by PCR from the detection of the first case, and uninfected contacts
#' averaged 3 tests during follow-up.
#' Symptoms were self-reported in a phone interview 20 to 60 days after the first case.
#'
#' \describe{
#'   \item{id}{Individual identifier;}
#'   \item{hh_id}{Household identifier;}
#'   \item{status}{Follow-up outcome: 0 escaped, 1 symptomatic, 2 asymptomatic
#'   (no symptoms reported during follow-up);}
#'   \item{age}{Age class: `adult` if older than 12 years (adolescents and adults),
#'   `child` otherwise;}
#'   \item{sex}{Sex, `F` or `M`;}
#'   \item{end_date}{Follow-up end: 21 days after the household's first positive PCR;}
#'   \item{test_date}{Date of the first positive test within the follow-up;}
#'   \item{vacc_status}{Vaccination status: `unvaccinated`, `1dose`, `2doses_leq90`,
#'   `2doses_gtr90` (second dose within or more than 90 days earlier) or `3doses`,
#'   counting the doses received more than 7 days before the household's first positive PCR;}
#'   \item{recovered}{Whether previously infected: reported in the interview,
#'   or tested positive before July 2021;}
#'   \item{isolated}{Whether the household took isolation measures,
#'   full or partial quarantine of its infected members.}
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
#' @format A data frame with 4,168 rows and 10 columns
#' @source <https://gitlab.pasteur.fr/mmmi-pasteur/hh_israel_sarscov2_delta2/-/raw/main/codes/db_filtered.txt>
#'
#' @details
#' Each household is followed for 21 days from its first positive SARS-CoV-2 PCR test
#' as in the original analysis.
#' A positive test before July 2021 counts as a prior infection, not as a case.
#' A member first testing positive after \code{end_date} counts as escaped.
#' In 135 households several members test positive on the first day.
#' In 2 households isolation differs between members.
#'
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_sarscov2_2021_icofs.R", package = "housepi")}
#'
#' Adults gave informed consent to the study \insertCite{cortier_factors_2025}{housepi}.
#' The source repository distributes the data under the GNU Lesser General Public License, version 3.
#'
#' @references
#' \insertRef{cortier_factors_2025}{housepi}
NULL

#' Flu 2014-15 household survey from Matsumoto, Japan
#'
#' Questionnaire survey of influenza in the households of primary school students
#' in Matsumoto city, Japan, during the 2014-15 season.
#' At the end of the season, parents of the students of all 29 public primary schools
#' reported the family members living with the student and which of them had influenza.
#' Original observational study from Uchida et al. (2017), analysed by Endo et al. (2019).
#'
#' \describe{
#'   \item{hh_id}{Survey identifier, one per student: a household with several
#'   students appears once per student;}
#'   \item{member}{Family member: `student` (the surveyed student), `sibling`,
#'   `father`, `mother` or `other` (mostly grandparents);}
#'   \item{status}{`1` if reported to have had influenza during the season, `0` otherwise.
#'   Students' cases were diagnosed at a medical institution, the other members' were
#'   reported by parents without a required diagnosis.}
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
#' Where more members were reported infected than the household held, the number of
#' members is raised to the number of infected, as in Endo et al. (2019).
#'
#' To reproduce the data processing, see the script at:
#' \code{system.file("data-raw", "make_flu_2014_matsumoto.R", package = "housepi")}
#'
#' @references
#' \insertRef{uchida_effectiveness_2017}{housepi}
#'
#' \insertRef{endo_fine_scale_2019}{housepi}
NULL
