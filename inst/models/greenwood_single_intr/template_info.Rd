\title{Greenwood with Single Introduction}
\name{housepi template : Greenwood with Single Introduction}

\description{
The Greenwood model assumes individual contacts to be homogeneous, independent and discrete.
Generation time units are discrete.
The probability of being infected in the successive generation depends only on the presence of at least one infectious individual in the previous generation.
The index case is assumed unique and unknown.
Given a final size of positive cases within household, all possible chains of transmission are considered.
The likelihood is the sum of the probability of each chain.
Missing data and asymptomatic status are not modeled.
\if{html}{\out{<hr>}}
}

\section{Likelihood}{
    \itemize{
    \item Person-to-person transmission probability: \eqn{1 - \exp(-{\beta})}
    }
\if{html}{\out{<hr>}}
}

\section{Parameters}{
    \itemize{
        \item \code{beta} : person to person transmission rate.
        }
\if{html}{\out{<hr>}}
}

\section{References}{
    1. Halloran ME, Longini IM, Struchiner CJ.
    Design and Analysis of Vaccine Studies.
    New York, NY: Springer New York; 2010.
    \url{https://doi.org/10.1007/978-0-387-68636-3}
}