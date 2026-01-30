\docType{class}
\title{Chain Binomial Heterogeneous}
\name{housepi template : Chain Binomial Heterogeneous}
\description{
The chain binomial model assumes transmission occurs at discrete generation time units.
The heterogeneous model assumes risk of transmission to depend on individuals features.
The index cases are not known a priori and represent the initial generation.
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
\deqn{L = P^{(+)}(\mathbf{x}, \mathbf{s}) P^{(-)}(\mathbf{x}, \mathbf{s})}
with household members \eqn{\mathbf{x} = \lbrace x_i \rbrace_{i=1}^{n}}, 
statuses \eqn{\mathbf{s} = \lbrace s_i \rbrace_{i=1}^{n}}
and sero-positive members \eqn{\mathbf{x'} = \lbrace x_i \rbrace_{s_i > 0}}
\itemize{
\item Sero-positive contribution :
\deqn{
P^{(+)}(\mathbf{x'}) =
1 - \sum_{k=0}^{\vert\mathbf{x'}\vert - 1} \sum_{m = 1}^\binom{\vert\mathbf{x'}\vert}{k} 
P^{(+)}(\mathbf{x'}, \mathbf{s'}^{(k,\,m)}) P^{(-)}(\mathbf{x'}, \mathbf{s'}^{(k,\,m)})}
where \eqn{\mathbf{s'}^{(k,\,m)}} is the \eqn{m}-th combination of \eqn{k} positive statuses over \eqn{\vert\mathbf{x'}\vert} members.
\item Sero-negative contribution :
\eqn{P^{(-)}(\mathbf{x}, \mathbf{s}) = 
\prod_{i : s_i<0}B(x_i) \prod_{j : s_j > 0} Q(x_j, x_i)
}
}
where
\itemize{
\item Community escape probability \eqn{i} :
\eqn{B(x_i) = \exp(-{\mu}^{(\rm com)}_{x_i} {\alpha})}
\item Person-to-person escape probability \eqn{i} from \eqn{j} :
\eqn{Q(x_j, x_i) = \exp(-{\mu}^{(\rm sus)}_{x_i} {\beta} {\mu}^{(\rm inf)}_{x_j})}
}
\if{html}{\out{<hr>}}
}
\section{References}{
Haber M, Longini IM, Cotsonis GA. 
Models for the Statistical Analysis of Infectious Disease Data. 
Biometrics. 1988;44: 163.
\url{https://doi.org/10.2307/2531904}
}
