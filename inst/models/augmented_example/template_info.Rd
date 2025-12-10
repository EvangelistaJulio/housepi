\docType{class}
\title{Augmented Example}
\name{housepi template : Augmented Example}
\description{
Dummy model used to illustrate the usage of augmented variables in housepi.
Each positive status individual has an associated augmented variable \code{proba} to be specified.
This example does not model transmission dynamics.
\if{html}{\out{<hr>}}
}
\section{Likelihood}{
Missing.
\if{html}{\out{<hr>}}
}
\section{Augmented}{
\itemize{
\item \code{proba} \eqn{p}{p}: an augmented variable (or function) to specify a probability associated to the individual.
This variable should be defined so that \eqn{p \in [0,\,1]}{p in [0, 1]}.
}
\if{html}{\out{<hr>}}
}
\section{Parameters}{
No parameter is required by the model.
Parameters could still be required by user definitions of augmented variables.
}
