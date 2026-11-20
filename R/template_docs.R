## housepi - Epidemiological Analysis of Household Outbreaks.
## Copyright 2024-2026 Institut Pasteur (Paris), France.
## SPDX-License-Identifier: GPL-3.0-or-later
##

#' Test
#'
#' @description
#' A minimal placeholder model with no transmission mechanism, for sanity checks, demos and end-to-end tests of housepi.
#' The observation likelihood is constant for every household, and simulation only fills latent variables.
#'
#' @section Parameters:
#' No parameters are required by the model.
#'
#' @section Transmission:
#' \deqn{p^\mathrm{transm} = 1}
#'
#' @seealso [browse_templates()]
#' @name template_test
NULL

#' Final Size
#'
#' @description
#' The final-size model follows \insertCite{longini_household_1982;textual}{housepi}.
#' Each member escapes infection from the community, and from each infected member of the household,
#' with its own escape probability.
#' The final size depends only on these escape probabilities, not on when infections occurred \insertCite{ball_unified_1986}{housepi}.
#' The likelihood of a household sums over every transmission chain that ends in the observed final size.
#' Infected members share a rate partition when their group and contact rates match.
#' When all members share one partition, the model is a Reed-Frost chain binomial \insertCite{frost_conceptions_1976}{housepi}:
#' the community infects the first generation, and each later generation is infected by the previous one.
#' A household is refused when its infected members fall into too many partitions to sum over.
#' Up to 15 partitions are accepted when each holds a single member, fewer when partitions hold several.
#'
#' @section Parameters:
#' * \code{alpha} \eqn{\alpha}{alpha}: baseline community transmission rate, cumulated over the outbreak
#' * \code{beta} \eqn{\beta}{beta}: baseline person-to-person transmission rate, cumulated over one generation
#' * \code{gamma} \eqn{\gamma}{gamma}: household-size scaling of person-to-person rate
#' * \code{hh_size_ref} \eqn{n_\mathrm{ref}}{n_ref}: size of reference household
#'
#' @section Group Rates:
#' * \code{community} \eqn{\mu^\mathrm{com}}{mu^com}: member factor on the community rate
#' * \code{susceptibility} \eqn{\mu^\mathrm{sus}}{mu^sus}: member factor on the person-to-person rate it receives
#' * \code{infectivity} \eqn{\mu^\mathrm{inf}}{mu^inf}: member factor on the person-to-person rate it exerts
#'
#' @section Transmission:
#' \deqn{p^\mathrm{transm} = p^\mathrm{esc} \, p^\mathrm{inf}}{p^transm = p^esc p^inf}
#'
#' for a household of size \eqn{n} where member \eqn{i} has status \eqn{s_i = 1}{s_i = 1} if infected, \eqn{s_i = 0}{s_i = 0} if escaped,
#' and \eqn{n^{+}}{n^+} members are infected.
#'
#' * Escaped contribution:
#' \eqn{p^\mathrm{esc} = \prod_{i:\,s_i = 0} q^\mathrm{com}_{i} \prod_{j:\,s_j = 1} q^\mathrm{pp}_{ji}}{p^esc = prod_{i: s_i = 0} q^com_i prod_{j: s_j = 1} q^pp_ji}
#' * Infected contribution \insertCite{longini_household_1982,addy_generalized_1991}{housepi}:
#' \eqn{p^\mathrm{inf} = 1 - \sum_{m=0}^{n^{+}-1} p_m}{p^inf = 1 - sum_{m = 0}^{n^+ - 1} p_m}
#' * Chains of smaller final size:
#' \eqn{p_m = \sum_{k=1}^{\binom{n^{+}}{m}} p^\mathrm{transm}\left(s^{(m,k)}\right)}{p_m = sum_{k = 1}^{choose(n^+, m)} p^transm(s^(m,k))},
#' the probability that an outbreak confined to the \eqn{n^{+}}{n^+} infected members infects exactly \eqn{m} of them,
#' where \eqn{s^{(m,k)}}{s^(m,k)} holds the statuses of the infected members only,
#' the \eqn{k}-th choice of \eqn{m} of them infected and the others escaped
#'
#' with escape probabilities
#'
#' * from the community:
#' \eqn{q^\mathrm{com}_{i} = \exp\left(-\alpha_{i}\right)}{q^com_i = exp(-alpha_i)},
#' with \eqn{\alpha_{i} = \alpha \, \mu^\mathrm{com}_{i}}{alpha_i = alpha mu^com_i} the community rate of member \eqn{i}
#' * from member \eqn{j} to member \eqn{i}:
#' \eqn{q^\mathrm{pp}_{ji} = \exp\left(-\beta \left(\tfrac{n_\mathrm{ref}}{n}\right)^{\gamma} \mu^\mathrm{inf}_{j} \, \mu^\mathrm{ctc}_{ji} \, \mu^\mathrm{sus}_{i}\right)}{q^pp_ji = exp(-beta (n_ref / n)^gamma mu^inf_j mu^ctc_ji mu^sus_i)}
#' * \eqn{\mu^\mathrm{ctc}_{ji}}{mu^ctc_ji}: contact rate from member \eqn{j} to member \eqn{i}, equal to 1 unless set by [add_contact_rate()]
#'
#' @section Model Quantities:
#' [inspect_state()] reports in `model_quantities`:
#' * \code{community_escape}: \eqn{q^\mathrm{com}_{i}}{q^com_i}, per member
#' * \code{p2p_escape}: \eqn{q^\mathrm{pp}_{ji}}{q^pp_ji}, infector \eqn{j} by row and infectee \eqn{i} by column, 1 on the diagonal
#'
#' @section Notes:
#' In a household of size \eqn{n} the baseline person-to-person rate is \eqn{\beta (n_\mathrm{ref}/n)^{\gamma}}{beta (n_ref / n)^gamma} \insertCite{cauchemez_bayesian_2004}{housepi},
#' so \eqn{\beta}{beta} is the baseline rate in a household of size \eqn{n_\mathrm{ref}}{n_ref}.
#'
#' @references
#' \insertAllCited{}
#' @seealso [browse_templates()]
#' @name template_final_size
NULL

#' Case Ascertained
#'
#' @description
#' The case-ascertained model assumes that each household member \eqn{i} is subject to an instantaneous risk of infection \eqn{\lambda_{i}(t)}{lambda_i(t)},
#' where time \eqn{t} is continuous, or counts days \eqn{d} in discrete time.
#' A susceptible member can get infected at time \eqn{t_{i}}{t_i}, or on day \eqn{d_{i}}{d_i},
#' according to an inhomogeneous Poisson process of rate equal to the risk of infection.
#' The household is followed from its first infection, the primary case, which enters the likelihood only as the start of exposure.
#'
#' @section Data:
#' * \code{end_date} \eqn{d_\mathrm{end}}{d_end}: last day of follow-up of the member, included
#'
#' It needs a value for every member, escaped ones included.
#' An observation standing for it is defined at household level, or from data columns only.
#' In simulation, the primary case is drawn uniformly among the members \insertCite{cauchemez_bayesian_2004}{housepi}, and infected on day 0.
#' \code{end_date} must count days from day 0. Calendar dates are not supported.
#' In discrete time, alongside the primary case, other members can be infected from the community on day 0.
#' In continuous time the primary case infection time is uniform within day 0.
#'
#' @section Parameters:
#' * \code{alpha} \eqn{\alpha}{alpha}: baseline community instantaneous transmission rate
#' * \code{beta} \eqn{\beta}{beta}: baseline person-to-person instantaneous transmission rate
#' * \code{gamma} \eqn{\gamma}{gamma}: household-size scaling of person-to-person rate
#' * \code{hh_size_ref} \eqn{n_\mathrm{ref}}{n_ref}: size of reference household
#'
#' @section Latent:
#' * \code{infec_time} \eqn{t_{i}}{t_i}: time of infection of an infected member, day \eqn{d_{i}}{d_i} in discrete time
#'
#' Infection times are generated by transmission: they take no `model`, only `min` and `max`.
#' Simulation needs no `infec_time` definition, and ignores it.
#'
#' The likelihood is unchanged when every time shifts by the same amount.
#' Inference samples infection times between `min` and `max`,
#' while observations inform about the timing.
#'
#' @section Group Rates:
#' * \code{community} \eqn{\mu^\mathrm{com}}{mu^com}: member factor on the community rate
#' * \code{susceptibility} \eqn{\mu^\mathrm{sus}}{mu^sus}: member factor on the person-to-person rate it receives
#' * \code{infectivity} \eqn{\mu^\mathrm{inf}}{mu^inf}: member factor on the person-to-person rate it exerts
#'
#' @section Transmission:
#' \deqn{p^\mathrm{transm} = p^\mathrm{esc} \, p^\mathrm{inf}}{p^transm = p^esc p^inf}
#'
#' for a household of size \eqn{n} where member \eqn{i} has status \eqn{s_i = 1}{s_i = 1} if infected, \eqn{s_i = 0}{s_i = 0} if escaped.
#'
#' * Instantaneous risk of infection:
#' \deqn{\lambda_{i}(t) = \alpha_{i} + \beta \left(\tfrac{n_\mathrm{ref}}{n}\right)^{\gamma} \mu^\mathrm{sus}_{i} \sum_{j:\,t_{j} < t} \mu^\mathrm{inf}_{j} \, \mu^\mathrm{ctc}_{ji} \, f_{j}(t - t_{j})}{lambda_i(t) = alpha_i + beta (n_ref / n)^gamma mu^sus_i sum_{j: t_j < t} mu^inf_j mu^ctc_ji f_j(t - t_j)}
#' * \eqn{\alpha_{i} = \alpha \, \mu^\mathrm{com}_{i}}{alpha_i = alpha mu^com_i}: community rate of member \eqn{i}
#' * \eqn{f_{j}}{f_j}: infectivity profile of member \eqn{j} as a function of time since its infection, set by [add_infectivity_profile()]
#' * \eqn{\mu^\mathrm{ctc}_{ji}}{mu^ctc_ji}: contact rate from member \eqn{j} to member \eqn{i}, equal to 1 unless set by [add_contact_rate()]
#' * \eqn{t_{0}}{t_0}: introduction time, earliest infection time in the household, that of the primary case
#'
#' \subsection{Real time}{
#' \itemize{
#' \item Infected contribution:
#' \eqn{p^\mathrm{inf} = \prod_{i:\,s_{i} = 1,\, t_{i} > t_{0}} \lambda_{i}(t_{i}) \exp\left(-\int_{t_{0}}^{t_{i}} \lambda_{i}(t) \, \mathrm{d}t\right)}{p^inf = prod_{i: s_i = 1, t_i > t_0} lambda_i(t_i) exp(-integral_{t_0}^{t_i} lambda_i(t) dt)}
#' \item Escaped contribution:
#' \eqn{p^\mathrm{esc} = \prod_{i:\,s_{i} = 0} \exp\left(-\int_{t_{0}}^{t_\mathrm{end}} \lambda_{i}(t) \, \mathrm{d}t\right)}{p^esc = prod_{i: s_i = 0} exp(-integral_{t_0}^{t_end} lambda_i(t) dt)},
#' where the follow-up ends at \eqn{t_\mathrm{end} = d_\mathrm{end} + 1}{t_end = d_end + 1}, the end of its last day
#' }
#' }
#' \subsection{Discrete time}{
#' The primary case is infected on day \eqn{d_{0}}{d_0}, the earliest infection day in the household.
#' Members infected on day \eqn{d_{0}}{d_0} form the set \eqn{\mathcal{I}_{0}}{I_0}, which holds the primary case.
#' On that day only the community infects, so \eqn{\lambda_{i}(d_{0}) = \alpha_{i}}{lambda_i(d_0) = alpha_i}.
#' Any member of \eqn{\mathcal{I}_{0}}{I_0} may be the primary case, so the first-day contribution sums over them
#' \insertCite{cauchemez_household_2009,zhang_frequent_2026}{housepi}, and the transmission probability becomes \eqn{p^\mathrm{transm} = p^\mathrm{esc} \, p^\mathrm{inf} \, p^\mathrm{prim}}{p^transm = p^esc p^inf p^prim}.
#' \itemize{
#' \item First-day contribution:
#' \eqn{p^\mathrm{prim} = \sum_{i \in \mathcal{I}_{0}} \prod_{j \in \mathcal{I}_{0},\, j \neq i} \left(1 - e^{-\lambda_{j}(d_{0})}\right)}{p^prim = sum_{i in I_0} prod_{j in I_0, j != i} (1 - exp(-lambda_j(d_0)))}
#' \item Infected contribution:
#' \eqn{p^\mathrm{inf} = \prod_{i:\,s_{i} = 1,\, d_{i} > d_{0}} \left(1 - e^{-\lambda_{i}(d_{i})}\right) \exp\left(-\sum_{d = d_{0}}^{d_{i} - 1} \lambda_{i}(d)\right)}{p^inf = prod_{i: s_i = 1, d_i > d_0} (1 - exp(-lambda_i(d_i))) exp(-sum_{d = d_0}^{d_i - 1} lambda_i(d))}
#' \item Escaped contribution:
#' \eqn{p^\mathrm{esc} = \prod_{i:\,s_{i} = 0} \exp\left(-\sum_{d = d_{0}}^{d_\mathrm{end}} \lambda_{i}(d)\right)}{p^esc = prod_{i: s_i = 0} exp(-sum_{d = d_0}^{d_end} lambda_i(d))}
#' }
#' }
#'
#' @section Model Quantities:
#' [inspect_state()] reports in `model_quantities` the terms of \eqn{\lambda_{i}(t)}{lambda_i(t)}:
#' * \code{community_rate}: \eqn{\alpha_{i}}{alpha_i}, per member
#' * \code{p2p_rate}: \eqn{\beta \left(\tfrac{n_\mathrm{ref}}{n}\right)^{\gamma} \mu^\mathrm{inf}_{j} \, \mu^\mathrm{ctc}_{ji} \, \mu^\mathrm{sus}_{i}}{beta (n_ref / n)^gamma mu^inf_j mu^ctc_ji mu^sus_i}, multiplying \eqn{f_{j}(t - t_{j})}{f_j(t - t_j)}, infector \eqn{j} by row and infectee \eqn{i} by column, 0 on the diagonal
#'
#' @section Notes:
#' In a household of size \eqn{n} the baseline person-to-person rate is \eqn{\beta (n_\mathrm{ref}/n)^{\gamma}}{beta (n_ref / n)^gamma} \insertCite{cauchemez_bayesian_2004}{housepi},
#' so \eqn{\beta}{beta} is the baseline rate in a household of size \eqn{n_\mathrm{ref}}{n_ref}.
#'
#' @references
#' \insertAllCited{}
#' @seealso [browse_templates()]
#' @name template_case_ascertained
NULL
