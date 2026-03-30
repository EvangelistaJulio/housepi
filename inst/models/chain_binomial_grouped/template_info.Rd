\docType{class}
\title{Chain Binomial Grouped}
\name{housepi template : Chain Binomial Grouped}
\description{
The grouped chain binomial model assumes transmission occurs at discrete generation time units.
Like the heterogeneous model, risk of transmission depends on individual features (group rates).
However, the algorithm exploits type degeneracy: when multiple individuals share the same
(community, susceptibility, infectivity) profile, the algorithm enumerates distinct
type-composition sequences per generation instead of all individual permutations.
This reduces computational cost from O(N!) to a function of the number of distinct types.
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} \eqn{\alpha} : community transmission rate
\item \code{beta} \eqn{\beta} : person to person transmission rate
}
\if{html}{\out{<hr>}}
}
\section{Group Rates}{
\itemize{
\item \code{community} \eqn{{\mu}^{(\rm com)}} : predisposition to infection from community relative to reference
\item \code{susceptibility} \eqn{{\mu}^{(\rm sus)}} : predisposition to infection relative to reference
\item \code{infectivity} \eqn{{\mu}^{(\rm inf)}} : predisposition to transmission relative to reference
}
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\deqn{L = P_{I}(\mathbf{x}, \mathbf{s}) P_{S}(\mathbf{x}, \mathbf{s})}
The sero-positive contribution is computed via generation-based enumeration
with type-composition weights (multinomial coefficients).
The sero-negative contribution is identical to the heterogeneous model.
\if{html}{\out{<hr>}}
}
\section{References}{
Endo A, et al.
Introduction to particle Markov-chain Monte Carlo for disease dynamics modellers.
Epidemics. 2019;29: 100363.
\url{https://doi.org/10.1371/journal.pcbi.1007589}
}
