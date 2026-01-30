\docType{class}
\title{Chain Binomial}
\name{housepi template : Chain Binomial}
\description{
The chain binomial model assumes transmission occurs at discrete generation time units.
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
}
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\deqn{L = P^{(+)}(n_{+}) P^{(-)}(n, n_{+})}
for a household of size \eqn{n} with
\itemize{
\item Number of sero-positive household members : \eqn{n_{+}}
\item Sero-positive contribution :
\eqn{P^{(+)}(n_{+}) 
= 
1 - \sum_{k=0}^{n_{+} - 1} \binom{n_{+}}{k} P^{(+)}(k) P^{(-)}(n_{+}, k)}
\item Sero-negative contribution :
\eqn{P^{(-)}(n, n_{+}) = (B Q^{n_{+}})^{n - n_{+}}}
}
where
\itemize{
\item Community escape probability :
\eqn{B = \exp(-\alpha)}
\item Person-to-person escape probability :
\eqn{Q = \exp(-\beta)}
}
\if{html}{\out{<hr>}}
}
\section{References}{
Longini IM, Koopman JS.
Household and Community Transmission Parameters from Final Distributions of Infections in Households.
Biometrics. 1982;38: 115-126.
\url{https://doi.org/10.2307/2530294}
}
