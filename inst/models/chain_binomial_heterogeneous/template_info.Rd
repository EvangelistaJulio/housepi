\docType{class}
\title{Chain Binomial Heterogeneous}
\name{housepi template : Chain Binomial Heterogeneous}
\description{
MISSING DESCRIPTION.
Missing data and asymptomatic status are not modeled.
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} : community transmission rate;
\item \code{beta} : person to person transmission rate;
\item \code{sus_woman} : susceptibility of women;
\item \code{sus_child} : susceptibility of children;
\item \code{inf_child} : infectivity of children;
\item \code{child_age} (c) : children/adult age threshold.
}
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\itemize{
\item Hazard of infection:
\deqn{
{\lambda}_{ki}(t) = {\mu}^{(\rm sus)}_{i} {\beta} {\mu}^{(\rm inf)}_{k}
}{
lambda_ki(t) = mu_sus_i * beta * mu_inf_k
}
\item `Infectivity` \eqn{{\mu}^{(\rm inf)}}{mu_inf}: 
\tabular{rrr}{
                       |\tab \strong{FALSE}   |\tab \strong{TRUE}       \cr
        \strong{child} |\tab \emph{reference} |\tab \code{inf_child} 
}
\item `Susceptibility` \eqn{{\mu}^{(\rm sus)}}{mu_sus}: 
\tabular{rrr}{
                       |\tab \strong{male}        |\tab \strong{female}  \cr
        \strong{adult} |\tab \emph{reference}     |\tab \code{sus_woman} \cr
        \strong{child} |\tab \code{sus_child}     |\tab \code{sus_child} 
}
}
\if{html}{\out{<hr>}}
}
\section{References}{
1. MISSING REFERENCE
}
