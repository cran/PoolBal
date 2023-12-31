% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alternatives.R
\name{rBetaH4}
\alias{rBetaH4}
\alias{rBetaH3}
\title{Generate realizations of beta alternative distributions}
\usage{
rBetaH4(a, b = 1/w + a * (1 - 1/w), w = (1 - a)/(b - a), M = 2, N = 10)

rBetaH3(
  a,
  b = 1/w + a * (1 - 1/w),
  w = (1 - a)/(b - a),
  eta = 0.5,
  M = 2,
  N = 10
)
}
\arguments{
\item{a}{first beta parameter, numeric between 0 and infinity}

\item{b}{second beta parameter, numeric  between 0 and infinity}

\item{w}{UMP parameter between 0 and 1}

\item{M}{number of p-values per realization}

\item{N}{number of realizations}

\item{eta}{numeric between 0 and 1, proportion of non-null
tests per sample}
}
\value{
An N by M matrix of simulated p-values.
}
\description{
These functions can be used to generate samples of
p-values all following a beta distribution (H4) or following
either uniform or beta distributions according to proportion
eta (H3).
}
\details{
Alternatives

These functions are provided as a convenience, and
support a/b (shape1/shape2) or a/w specification of beta
parameters.
}
\section{Functions}{
\itemize{
\item \code{rBetaH4()}: iid Beta(a,w) p-values

\item \code{rBetaH3()}: M*eta iid Beta(a,w) p-values, others uniform

}}
\examples{
rBetaH4(a = 0.5, b = 1.5, M = 10, N = 100)
rBetaH3(a = 0.5, b = 1.5, eta = 0.5, M = 10, N = 100)
}
\author{
Chris Salahub
}
