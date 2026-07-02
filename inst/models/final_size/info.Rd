\docType{class}
\title{Final Size}
\name{housepi template : Final Size}
\description{
The final size model assumes transmission occurs at discrete generation time units.
The homogeneous model assumes risk of transmission to be equal for all individuals.
The index cases are not known a priori and represent the initial generation.
Infections in this generation are modeled with according community transmission probability.
Transmission does not account for community exposure after the initial generation.
The person-to-person transmission probability  in the successive generation depends only on the infectious individuals of the previous generation.
Given a final size of positive cases within household, all possible chains of transmission are then considered.
Finally, the likelihood of an observed household outbreak is the sum of the probability of each chain.
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} \eqn{\alpha}{alpha} : community transmission rate
\item \code{beta} \eqn{\beta}{beta} : person to person transmission rate
\item \code{gamma} \eqn{\gamma} : parameter to model dilution of contacts
\item \code{size_ref} \eqn{n_{\ast}} : size of reference household (required constant)
}
\if{html}{\out{<hr>}}
}
\section{Transmission}{
\deqn{p^\mathrm{transm} = p_\mathrm{I}(n_\mathrm{I}) \, p_\mathrm{S}(n, n_\mathrm{I})}
for a household of size \eqn{n} with
\itemize{
\item Number of sero-positive household members : \eqn{n_\mathrm{I}}
\item Sero-positive contribution :
\eqn{p_\mathrm{I}(n_\mathrm{I})
=
1 - \sum_{k=0}^{n_\mathrm{I} - 1} \binom{n_\mathrm{I}}{k}
    p_\mathrm{I}(k) \, p_\mathrm{S}(n_\mathrm{I}, k)}
\item Sero-negative contribution :
\eqn{p_\mathrm{S}(n, n_\mathrm{I})
= \left(q^\mathrm{com} \, (q^\mathrm{p2p})^{n_\mathrm{I}}\right)^{n - n_\mathrm{I}}}
}
where
\itemize{
\item Community escape probability :
\eqn{q^\mathrm{com} = \exp(-\alpha)}
\item Person-to-person escape probability :
\eqn{q^\mathrm{p2p} = \exp(-\beta)}
}
\if{html}{\out{<hr>}}
}
\section{References}{
\insertRef{longini_household_1982}{housepi}
}
