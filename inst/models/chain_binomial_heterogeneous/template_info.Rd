\docType{class}
\title{Chain Binomial Heterogeneous}
\name{housepi template : Chain Binomial Heterogeneous}
\description{
The chain binomial model assumes transmission occurs at discrete generation time units.
The heterogeneous model assumes risk of transmission to depend on individuals features.
The index cases are not known a priori and represent the initial generation.
See [Chain Binomial](hhmodel_chain_binomial.html) for more details.
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\itemize{
\item Community transmission probability \eqn{i}{i}:
\eqn{1 - \exp(-{\mu}^{(\rm com)}_{i} {\alpha})}{1 - exp(-mu_com_i * alpha)}
\item Person-to-person transmission probability \eqn{k}{k} to \eqn{i}{i}:
\eqn{1 - \exp(-{\mu}^{(\rm sus)}_{i} {\beta} {\mu}^{(\rm inf)}_{k})}{1 - exp(-mu_sus_i * beta * mu_inf_k)}
}
\if{html}{\out{<hr>}}
}
\section{Group Rates}{
\itemize{
\item \code{infectivity} \eqn{{\mu}^{(\rm inf)}}{mu_inf}: individual predisposition of transmission relative to reference;
\item \code{susceptibility} \eqn{{\mu}^{(\rm sus)}}{mu_sus}: individual predisposition of infection relative to reference;
\item \code{community} \eqn{{\mu}^{(\rm com)}}{mu_com}: individual predisposition of infection from community relative to reference;
}
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} \eqn{\alpha}{alpha} : community transmission rate (baseline);
\item \code{beta} \eqn{\beta}{beta} : person to person transmission rate (baseline).
}
\if{html}{\out{<hr>}}
}
\section{References}{
Haber M, Longini IM, Cotsonis GA. 
Models for the Statistical Analysis of Infectious Disease Data. 
Biometrics. 1988;44: 163.
\url{https://doi.org/10.2307/2531904}
}
