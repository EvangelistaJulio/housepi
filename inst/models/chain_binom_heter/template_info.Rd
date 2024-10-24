\title{Chain Binomial Heterogeneous}
\name{housepi template : Chain Binomial Heterogeneous}
\description{
MISSING DESCRIPTION.
Missing data and asymptomatic status are not modeled.
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
    \itemize{
    \item Hazard of infection: \deqn{
        {\lambda}_{ki}(t) = {\mu}^{(\rm sus)}_{i} {\beta} {\mu}^{(\rm inf)}_{k} {\mu}^{(\rm con)}_{k{\to}i}
        }
    \item `Susceptibility` \eqn{{\mu}^{(\rm sus)}}: \tabular{rrr}{
                       |\tab \strong{male}        |\tab \strong{female}      \cr
        \strong{adult} |\tab \emph{reference}     |\tab \code{sus_women} \cr
        \strong{child} |\tab \code{sus_child} |\tab \code{sus_child} 
        }
    \item `Infectivity` \eqn{{\mu}^{(\rm inf)}}: \tabular{rrr}{
                       |\tab \strong{FALSE}   |\tab \strong{TRUE}       \cr
        \strong{child} |\tab \emph{reference} |\tab \code{inf_child} 
        }
    \item `Contact` \eqn{{\mu}^{(\rm con)}}: \tabular{rrr}{
        \strong{from \\ to} |\tab \strong{adult}               |\tab \strong{child}               \cr
        \strong{adult}      |\tab \emph{reference}             |\tab \code{cont_adult_child} \cr
        \strong{child}      |\tab \code{cont_child_adult} |\tab \code{cont_child_child}
        }
    }
\if{html}{\out{<hr>}}
}

\section{Parameters}{
    \itemize{
        \item \code{alpha} : community transmission rate;
        \item \code{beta} : person to person transmission rate;
        \item \code{sus_women} : susceptibility of women;
        \item \code{sus_child} : susceptibility of children;
        \item \code{inf_child} : infectivity of children;
        \item \code{cont_adult_child} (c) : contact rate from adult to child;
        \item \code{cont_child_adult} (c) : contact rate from child to adult;
        \item \code{cont_child_child} (c) : contact rate from adult to child;
        \item \code{child_age} (c) : children/adult age threshold.
        }
\if{html}{\out{<hr>}}
}

\section{References}{
    1. MISSING REFERENCE
}