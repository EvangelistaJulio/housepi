\docType{class}
\title{Chain Binomial}
\name{housepi template : Chain Binomial}
\description{
The chain binomial model assumes individual contacts to be homogeneous, independent and discrete.
Generation time units are discrete.
The probability of being infected in the successive generation depends only ont the infectious individuals of the previous generation.
Index cases are unknown and due to transmission by community.
Given a final size of positive cases within household, all possible chains of transmission are considered.
The likelihood is the sum of the probability of each chain.
Missing data and asymptomatic status are not modeled.
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\itemize{
\item Community transmission probability: \eqn{1 - \exp(-{\alpha})}{1 - exp(-alpha)}
\item Person-to-person transmission probability: \eqn{1 - \exp(-{\beta})}{1 - exp(-beta)}
}
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} : community transmission rate;
\item \code{beta} : person to person transmission rate.
}
\if{html}{\out{<hr>}}
}
\section{References}{
1. Longini IM, Koopman JS.
Household and Community Transmission Parameters from Final Distributions of Infections in Households.
Biometrics. 1982;38: 115-126.
\url{https://doi.org/10.2307/2530294}
}