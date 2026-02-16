\docType{class}
\title{Augmented Infection Time}
\name{housepi template : Augmented Infection Time}
\description{
This model assumed that each household member \eqn{i} is suject to a force of infection \eqn{{\lambda}_{i}(t)}{lambda_i(t)},
where time \eqn{t} can take continuous or discrete values.
A susceptible member can get infected at time/day \eqn{t_{i}} according to a inhomogeneous Poisson process of rate equal to the force infection.
\if{html}{\out{<hr>}}
}
\section{Data}{
\itemize{
\item \code{end_followup} \eqn{t_{\rm end}} : end date of followup study
}
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} \eqn{\alpha} : community transmission rate
\item \code{beta} \eqn{\beta} : person to person transmission rate
\item \code{gamma} \eqn{\gamma} : parameter to model the dilution of contacts
\item \code{size_ref} \eqn{n_{\ast}} : the size of the reference household (required constant)
}
\if{html}{\out{<hr>}}
}
\section{Group Rates}{
\itemize{
\item \code{community} \eqn{{\mu}^{(\rm com)}} : predisposition to infection from community relative to reference
\item \code{susceptibility} \eqn{{\mu}^{(\rm sus)}}: predisposition to infection relative to reference
\item \code{infectivity} \eqn{{\mu}^{(\rm inf)}} : predisposition to transmission relative to reference
}
\if{html}{\out{<hr>}}
}
\section{Augmented}{
\itemize{
\item \code{infec_time} \eqn{t_{i}} : time of infection of infectious member
\item \code{onset_date} : time of symptom onset of infectious member
\item \code{incub_period} : time between infection and symptom onset
}
\if{html}{\out{<hr>}}
}
\section{Likelihood}{

\deqn{L = P_{I} \cdot P_{S}}

If \eqn{s_i} is the status of household member \eqn{i}, and

\itemize{
\item Force of infection for a household of size \eqn{n} : \deqn{
{\lambda}_{i}(t) = {\mu}^{(\rm com)}_{i} {\alpha} + {\mu}^{(\rm sus)}_{i} {\beta} \left(\frac{n_{\ast}}{n}\right)^{\gamma} \sum_{j\,:\,t_{j}<t_{i}} {\mu}^{(\rm inf)}_{j} f_{j}(t)
}
\item \code{infectivity profile} \eqn{f_{i}} : can be specified
\item \eqn{t_{0}} : introduction time (earliest infection time in the household)
}

\subsection{Real time}{
\itemize{
\item Infectious contribution :
\eqn{
P_{I} =\prod_{{i}:s_{i}>0} \lambda_{i}(t_{i})\exp{-\int_{t_{0}}^{t_{i}}\lambda_{i}(t){\rm d}t}
}
\item Susceptible contribution :
\eqn{
P_{S} = \prod_{{i}:s_{i}<0} \exp{-\int_{t_{0}}^{t_{\rm end}}\lambda_{i}(t){\rm d}t}
}
}
}
\subsection{Discrete time}{
\itemize{
\item Infectious contribution :
\eqn{
P_{ I} = \prod_{{i}:s_{i}>0}\left(1-\exp{-\lambda_{i}(t_{i})}\right)\exp{-\sum_{t={t}_{0}}^{t_{i}-1} \lambda_{i}(t)}
}
\item Susceptible contribution :
\eqn{
{P}_{S} =\prod_{{i}:s_{i}<0}\exp{-\sum_{t={t}_{0}}^{t_{\rm end}}\lambda_{i}(t)}
}
}
}
\if{html}{\out{<hr>}}
}
\section{References}{
Cauchemez S, Carrat F, Viboud C, Valleron AJ, Boelle PY.
A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data.
Statist Med. 2004;23: 3469-3487.
\url{https://doi.org/10.1002/sim.1912}
}
