## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Test
#'
#' @description
#' A minimal placeholder model used to test the housepi infrastructure without committing to a specific transmission mechanism.
#' The observation likelihood is constant for every household, and simulation only fills latent variables.
#' This template is intended for sanity checks, demos and end-to-end tests.
#'
#' @section Parameters:
#' No parameters are required by the model.
#'
#' @section Transmission:
#' \deqn{p^\mathrm{transm} = 1}
#'
#' @seealso [browse_templates()]
#' @name template_test
#' @keywords internal
NULL

#' Final Size
#'
#' @description
#' The final-size model is a Reed-Frost chain binomial: transmission proceeds in discrete generations,
#' and each generation's new infections are independent binomial draws among the remaining susceptibles.
#' The index cases are not known a priori: the first generation is infected from the community.
#' Each later generation is infected only by the members infected in the previous one.
#' The likelihood of a household sums over every transmission chain that ends in the observed final size.
#'
#' @section Parameters:
#' - \code{alpha} \eqn{\alpha}{alpha} : baseline community transmission risk
#' - \code{beta} \eqn{\beta}{beta} : baseline person-to-person transmission risk
#' - \code{gamma} \eqn{\gamma}{gamma} : household-size scaling of person-to-person risk
#' - \code{hh_size_ref} \eqn{n_\mathrm{ref}}{n_ref} : size of reference household
#'
#' @section Group Risks:
#' - \code{community} \eqn{\mu^\mathrm{com}}{mu^com} : predisposition to infection from community relative to reference
#' - \code{susceptibility} \eqn{\mu^\mathrm{sus}}{mu^sus} : predisposition to infection relative to reference
#' - \code{infectivity} \eqn{\mu^\mathrm{inf}}{mu^inf} : predisposition to transmission relative to reference
#'
#' @section Transmission:
#' \deqn{p^\mathrm{transm} = p^\mathrm{esc} \, p^\mathrm{inf}}{p^transm = p^esc p^inf}
#'
#' for a household of size \eqn{n} where member \eqn{i} has status \eqn{s_i = 1}{s_i = 1} if infected, \eqn{s_i = 0}{s_i = 0} if escaped,
#' and \eqn{n^{+}}{n^+} members are infected.
#'
#' - Escaped contribution :
#' \eqn{p^\mathrm{esc} = \prod_{i:\,s_i = 0} q^\mathrm{com}_{i} \prod_{j:\,s_j = 1} q^\mathrm{pp}_{ji}}{p^esc = prod_{i: s_i = 0} q^com_i prod_{j: s_j = 1} q^pp_ji}
#' - Infected contribution \insertCite{longini_household_1982,addy_generalized_1991}{housepi} :
#' \eqn{p^\mathrm{inf} = 1 - \sum_{m=0}^{n^{+}-1} P_m}{p^inf = 1 - sum_{m = 0}^{n^+ - 1} P_m}
#' - Chains of smaller final size :
#' \eqn{P_m = \sum_{k=1}^{\binom{n^{+}}{m}} p^\mathrm{transm}\left(s^{(m,k)}\right)}{P_m = sum_{k = 1}^{choose(n^+, m)} p^transm(s^(m,k))},
#' where \eqn{s^{(m,k)}}{s^(m,k)} keeps the \eqn{k}-th choice of \eqn{m} infected among the \eqn{n^{+}}{n^+} and sets the others to escaped
#'
#' with escape probabilities
#'
#' - from the community :
#' \eqn{q^\mathrm{com}_{i} = \exp\left(-\alpha \, \mu^\mathrm{com}_{i}\right)}{q^com_i = exp(-alpha mu^com_i)}
#' - from member \eqn{j} to member \eqn{i} :
#' \eqn{q^\mathrm{pp}_{ji} = \exp\left(-\beta \left(\tfrac{n_\mathrm{ref}}{n}\right)^{\gamma} \mu^\mathrm{inf}_{j} \, \mu^\mathrm{ctc}_{ji} \, \mu^\mathrm{sus}_{i}\right)}{q^pp_ji = exp(-beta (n_ref / n)^gamma mu^inf_j mu^ctc_ji mu^sus_i)}
#' - \eqn{\mu^\mathrm{ctc}_{ji}}{mu^ctc_ji} : contact risk from member \eqn{j} to member \eqn{i}, equal to 1 unless set by [add_contact_risk()]
#'
#' @section Model Quantities:
#' [inspect_house_state()] reports in `model_quantities`:
#' - \code{community_escape} : \eqn{q^\mathrm{com}_{i}}{q^com_i}, per member
#' - \code{p2p_escape} : \eqn{q^\mathrm{pp}_{ji}}{q^pp_ji}, infector \eqn{j} by row and infectee \eqn{i} by column, 1 on the diagonal
#'
#' @section Notes:
#' In a household of size \eqn{n} the baseline person-to-person risk is \eqn{\beta (n_\mathrm{ref}/n)^{\gamma}}{beta (n_ref / n)^gamma} \insertCite{cauchemez_bayesian_2004}{housepi},
#' so \eqn{\beta}{beta} is the baseline risk in a household of size \eqn{n_\mathrm{ref}}{n_ref}.
#' The choice of \eqn{n_\mathrm{ref}}{n_ref} leaves the likelihood unchanged, but \eqn{\beta}{beta} and \eqn{\gamma}{gamma} can become strongly correlated and chains mix slowly.
#' Choose \eqn{n_\mathrm{ref}}{n_ref} from [measure_geom_mean_hh_size()] to keep them nearly uncorrelated.
#'
#' @references
#' \insertAllCited{}
#' @seealso [browse_templates()]
#' @name template_final_size
#' @keywords internal
NULL

#' Case Ascertained
#'
#' @description
#' The case-ascertained model assumes that each household member \eqn{i} is subject to a force of infection \eqn{{\lambda}_{i}(t)}{lambda_i(t)},
#' where time \eqn{t} can take continuous or discrete values.
#' A susceptible member can get infected at time/day \eqn{t_{i}}{t_i} according to an inhomogeneous Poisson process of rate equal to the force of infection.
#' The household is followed from its first infection, the primary case, which enters the likelihood only as the start of exposure.
#'
#' @section Data:
#' - \code{start_date} : first day of follow-up of the member
#' - \code{end_date} \eqn{t_\mathrm{end}}{t_end} : last day of follow-up of the member, included
#' - \code{index_case} : 1 for the index case of the household, 0 for every other member
#'
#' All three need a value for every member, escaped ones included.
#' Each household has exactly one index case, so households with several are dropped beforehand.
#' In simulation the index case is the primary case.
#' An observation standing for them is defined at household level, or from data columns only.
#'
#' @section Parameters:
#' - \code{alpha} \eqn{\alpha}{alpha} : baseline community instantaneous transmission risk
#' - \code{beta} \eqn{\beta}{beta} : baseline person-to-person instantaneous transmission rate
#' - \code{gamma} \eqn{\gamma}{gamma} : household-size scaling of person-to-person risk
#' - \code{hh_size_ref} \eqn{n_\mathrm{ref}}{n_ref} : size of reference household
#'
#' @section Latent:
#' - \code{infec_time} \eqn{t_{i}}{t_i} : time of infection of infectious member
#'
#' Infection times are generated by transmission: they take no `model`, only `min` and `max`.
#' In simulation, transmission sets every infection time but the primary case's, which then requires `value`.
#'
#' @section Group Risks:
#' - \code{community} \eqn{\mu^\mathrm{com}}{mu^com} : predisposition to infection from community relative to reference
#' - \code{susceptibility} \eqn{\mu^\mathrm{sus}}{mu^sus} : predisposition to infection relative to reference
#' - \code{infectivity} \eqn{\mu^\mathrm{inf}}{mu^inf} : predisposition to transmission relative to reference
#'
#' @section Transmission:
#' \deqn{p^\mathrm{transm} = p^\mathrm{esc} \, p^\mathrm{inf}}{p^transm = p^esc p^inf}
#'
#' for a household of size \eqn{n} where member \eqn{i} has status \eqn{s_i = 1}{s_i = 1} if infected, \eqn{s_i = 0}{s_i = 0} if escaped.
#'
#' - Force of infection :
#' \deqn{\lambda_{i}(t) = \alpha \, \mu^\mathrm{com}_{i} + \beta \left(\tfrac{n_\mathrm{ref}}{n}\right)^{\gamma} \mu^\mathrm{sus}_{i} \sum_{j:\,t_{j} < t} \mu^\mathrm{inf}_{j} \, \mu^\mathrm{ctc}_{ji} \, f_{j}(t - t_{j})}{lambda_i(t) = alpha mu^com_i + beta (n_ref / n)^gamma mu^sus_i sum_{j: t_j < t} mu^inf_j mu^ctc_ji f_j(t - t_j)}
#' - \eqn{f_{j}}{f_j} : infectivity profile of member \eqn{j} as a function of time since its infection, set by [add_infectivity_profile()]
#' - \eqn{\mu^\mathrm{ctc}_{ji}}{mu^ctc_ji} : contact risk from member \eqn{j} to member \eqn{i}, equal to 1 unless set by [add_contact_risk()]
#' - \eqn{t_{0}}{t_0} : introduction time, earliest infection time in the household
#' - \eqn{i_{0}}{i_0} : primary case, infected at \eqn{t_{0}}{t_0} (on ties, the first in household order)
#'
#' \subsection{Real time}{
#' \itemize{
#' \item Infected contribution :
#' \eqn{p^\mathrm{inf} = \prod_{i:\,s_{i} = 1,\, i \neq i_{0}} \lambda_{i}(t_{i}) \exp\left(-\int_{t_{0}}^{t_{i}} \lambda_{i}(t) \, \mathrm{d}t\right)}{p^inf = prod_{i: s_i = 1, i != i_0} lambda_i(t_i) exp(-integral_{t_0}^{t_i} lambda_i(t) dt)}
#' \item Escaped contribution :
#' \eqn{p^\mathrm{esc} = \prod_{i:\,s_{i} = 0} \exp\left(-\int_{t_{0}}^{t_\mathrm{end} + 1} \lambda_{i}(t) \, \mathrm{d}t\right)}{p^esc = prod_{i: s_i = 0} exp(-integral_{t_0}^{t_end + 1} lambda_i(t) dt)}
#' }
#' }
#' \subsection{Discrete time}{
#' \itemize{
#' \item Infected contribution :
#' \eqn{p^\mathrm{inf} = \prod_{i:\,s_{i} = 1,\, i \neq i_{0}} \left(1 - \exp\left(-\lambda_{i}(t_{i})\right)\right) \exp\left(-\sum_{t = t_{0} + 1}^{t_{i} - 1} \lambda_{i}(t)\right)}{p^inf = prod_{i: s_i = 1, i != i_0} (1 - exp(-lambda_i(t_i))) exp(-sum_{t = t_0 + 1}^{t_i - 1} lambda_i(t))}
#' \item Escaped contribution :
#' \eqn{p^\mathrm{esc} = \prod_{i:\,s_{i} = 0} \exp\left(-\sum_{t = t_{0} + 1}^{t_\mathrm{end}} \lambda_{i}(t)\right)}{p^esc = prod_{i: s_i = 0} exp(-sum_{t = t_0 + 1}^{t_end} lambda_i(t))}
#' }
#' }
#'
#' @section Model Quantities:
#' [inspect_house_state()] reports in `model_quantities` the terms of \eqn{\lambda_{i}(t)}{lambda_i(t)}:
#' - \code{community_rate} : \eqn{\alpha \, \mu^\mathrm{com}_{i}}{alpha mu^com_i}, per member
#' - \code{p2p_rate} : \eqn{\beta \left(\tfrac{n_\mathrm{ref}}{n}\right)^{\gamma} \mu^\mathrm{inf}_{j} \, \mu^\mathrm{ctc}_{ji} \, \mu^\mathrm{sus}_{i}}{beta (n_ref / n)^gamma mu^inf_j mu^ctc_ji mu^sus_i}, multiplying \eqn{f_{j}(t - t_{j})}{f_j(t - t_j)}, infector \eqn{j} by row and infectee \eqn{i} by column, 0 on the diagonal
#'
#' @section Notes:
#' In a household of size \eqn{n} the baseline person-to-person risk is \eqn{\beta (n_\mathrm{ref}/n)^{\gamma}}{beta (n_ref / n)^gamma} \insertCite{cauchemez_bayesian_2004}{housepi},
#' so \eqn{\beta}{beta} is the baseline risk in a household of size \eqn{n_\mathrm{ref}}{n_ref}.
#' The choice of \eqn{n_\mathrm{ref}}{n_ref} leaves the likelihood unchanged, but \eqn{\beta}{beta} and \eqn{\gamma}{gamma} can become strongly correlated and chains mix slowly.
#' Choose \eqn{n_\mathrm{ref}}{n_ref} from [measure_geom_mean_hh_size()] to keep them nearly uncorrelated.
#'
#' @references
#' \insertAllCited{}
#' @seealso [browse_templates()]
#' @name template_case_ascertained
#' @keywords internal
NULL
