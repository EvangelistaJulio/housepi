\docType{class}
\title{Augmented Example}
\name{housepi template : Augmented Example}
\description{
Dummy model used to illustrate the usage of augmented variables in housepi.
Each positive status individual has an associated augmented variable \code{proba} to be specified.
This example does not model transmission dynamics.
\if{html}{\out{<hr>}}
}
\section{Parameters}{
No parameters are required by the model.
\if{html}{\out{<hr>}}
}
\section{Augmented}{
\itemize{
\item \code{proba} \eqn{p} : probability associated to each member.
The definition should imply \eqn{p \in [0,\,1]}.
}
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
\deqn{L = \sum_{s_i > 0} p_i}
where \eqn{s_i} is the status of household member \eqn{i}.
}


