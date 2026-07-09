\docType{class}
\title{Case Ascertained}
\name{housepi template : Case Ascertained}
\description{
This model assumed that each household member \eqn{i} is suject to a force of infection \eqn{{\lambda}_{i}(t)}{lambda_i(t)},
where time \eqn{t} can take continuous or discrete values.
A susceptible member can get infected at time/day \eqn{t_{i}} according to a inhomogeneous Poisson process of rate equal to the force infection.
\if{html}{\out{<hr>}}
}
\section{Data}{
\itemize{
\item \code{end_followup} \eqn{t_\mathrm{end}} : end date of followup study
}
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} \eqn{\alpha} : community transmission rate
\item \code{beta} \eqn{\beta} : person to person transmission rate
\item \code{gamma} \eqn{\gamma} : parameter to model dilution of contacts
\item \code{hhsize_ref} \eqn{n_{\ast}} : size of reference household (required constant)
}
\if{html}{\out{<hr>}}
}
\section{Augmented}{
\itemize{
\item \code{infec_time} \eqn{t_{i}} : time of infection of infectious member
}
\if{html}{\out{<hr>}}
}
\section{Group Rates}{
\itemize{
\item \code{community} \eqn{\mu^\mathrm{com}} : predisposition to infection from community relative to reference
\item \code{susceptibility} \eqn{\mu^\mathrm{sus}}: predisposition to infection relative to reference
\item \code{infectivity} \eqn{\mu^\mathrm{inf}} : predisposition to transmission relative to reference
}
\if{html}{\out{<hr>}}
}
\section{Transmission}{

\deqn{p^\mathrm{transm} = p_\mathrm{I} \cdot p_\mathrm{S}}

If \eqn{s_i} is the status of household member \eqn{i}, and

\itemize{
\item Force of infection for a household of size \eqn{n} : \deqn{
{\lambda}_{i}(t) = \mu^\mathrm{com}_{i} {\alpha} + \mu^\mathrm{sus}_{i} {\beta} \left(\tfrac{n_{\ast}}{n}\right)^{\gamma} \sum_{j\,:\,t_{j}<t_{i}} \mu^\mathrm{inf}_{j} f_{j}(t)
}
\item \code{infectivity profile} \eqn{f_{i}} : can be specified
\item \eqn{t_{0}} : introduction time (earliest infection time in the household)
}

\subsection{Real time}{
\itemize{
\item Infectious contribution :
\eqn{
p_\mathrm{I} =\prod_{{i}:s_{i}>0} \lambda_{i}(t_{i})\exp{-\int_{t_{0}}^{t_{i}}\lambda_{i}(t)\mathrm{d}t}
}
\item Susceptible contribution :
\eqn{
p_\mathrm{S} = \prod_{{i}:s_{i}<0} \exp{-\int_{t_{0}}^{t_\mathrm{end}}\lambda_{i}(t)\mathrm{d}t}
}
}
}
\subsection{Discrete time}{
\itemize{
\item Infectious contribution :
\eqn{
p_\mathrm{I} = \prod_{{i}:s_{i}>0}\left(1-\exp{-\lambda_{i}(t_{i})}\right)\exp{-\sum_{t={t}_{0}}^{t_{i}-1} \lambda_{i}(t)}
}
\item Susceptible contribution :
\eqn{
p_\mathrm{S} =\prod_{{i}:s_{i}<0}\exp{-\sum_{t={t}_{0}}^{t_\mathrm{end}}\lambda_{i}(t)}
}
}
}
\if{html}{\out{<hr>}}
}
\section{References}{
\insertRef{cauchemez_bayesian_2004}{housepi}
}
