\docType{class}
\title{Augmented Infection Time}
\name{housepi template : Augmented Infection Time}
\description{
This model assumed that each household member \eqn{i}{i} is suject to a force of infection \eqn{{\lambda}_{i}(t)}{lambda_i(t)},
where time \eqn{t}{t} can take continuous or discrete values.
A susceptible member can get infected at time/day \eqn{t_{i}}{t_i} according to a inhomogeneous Poisson process of rate equal to the force infection.
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\itemize{
\item Force of infection: \deqn{
{\lambda}_{i}(t) = {\mu}^{(\rm sus)}_{i} {\beta} \left(\frac{n_{\ast}}{n}\right)^{-{\gamma}} \sum_{k\,:\,t_{k}<t_{i}} {\mu}^{(\rm inf)}_{k} f_{k}(t)
}{
lambda_i(t) = mu_sus_i * beta * (n_star/n)^(-gamma) sum_{k: t_k < t_i} mu_inf_k f_k(t)
};
\item Infectivity profile: \eqn{f_{k}}{f_k}.
}
\if{html}{\out{<hr>}}
}
\section{Group Rates}{
\itemize{
\item \code{infectivity} \eqn{{\mu}^{(\rm inf)}}{mu_inf}: individual predisposition of transmission relative to reference;
\item \code{susceptibility} \eqn{{\mu}^{(\rm sus)}}{mu_sus}: individual predisposition of infection relative to reference.
}
\if{html}{\out{<hr>}}
}
\section{Augmented}{
\itemize{
\item \code{infec_time} : time of symptom onset of infectious individual;
\item \code{onset_date} : time of symptom onset of infectious individual;
\item \code{incub_period} : time between infection and symptom onset.
}
Note: by including an external end of the followup, \code{onset_date} and \code{incub_period} won't be required anymore.
This change will also drop \code{study_period} (see below).
\if{html}{\out{<hr>}}
}
\section{Parameters}{
\itemize{
\item \code{alpha} \eqn{\alpha}{alpha} : community transmission rate;
\item \code{beta} \eqn{\beta}{beta} : person to person transmission rate (baseline);
\item \code{size_ref} (c) \eqn{n_{\ast}}{n_star} : the size of the reference household (it must be constant);
\item \code{gamma} \eqn{\gamma}{gamma} : parameter to model the dilution of contacts;
\item \code{study_period} (c) : number of days of the followup period after inclusion (it must be constant).
}
\if{html}{\out{<hr>}}
}
\section{References}{
Cauchemez S, Carrat F, Viboud C, Valleron AJ, Boelle PY.
A Bayesian MCMC approach to study transmission of influenza: application to household longitudinal data.
Statist Med. 2004;23: 3469-3487.
\url{https://doi.org/10.1002/sim.1912}
}
