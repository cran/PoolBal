% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centrality.R
\name{estimatePrb}
\alias{estimatePrb}
\title{Compute the marginal rejection level}
\usage{
estimatePrb(
  poolFun,
  alpha = 0.05,
  b = 1,
  M = 2,
  interval = c(0, b),
  poolArgs = list(),
  ...
)
}
\arguments{
\item{poolFun}{function accepting a vector of p-values}

\item{alpha}{numeric between 0 and 1}

\item{b}{numeric, the value of the M - 1 repeated p-values}

\item{M}{integer, how many p-values are there?}

\item{interval}{two numerics giving the bounds of root-searching}

\item{poolArgs}{(optional) additional named arguments for poolFun}

\item{...}{additional arguments to uniroot}
}
\value{
The uniroot output.
}
\description{
Estimates the marginal rejection level for an
arbitrary pooled p-value function.
}
\details{
The marginal rejection level is the maximum p-value
in a single test less than b which still results in rejection of
the null when all other tests have a p-value of b.

This function is essentially a wrapper for uniroot, and accepts
a pooling function which takes a numeric vector as its first
argument and potentially other arguments given in poolArgs and
returns a single value. Using this pooling function, a specified
dimension M and a rejection level alpha, uniroot searches for the
root to poolFun - alpha along one margin when all other p-values
are equal to b.
}
\examples{
stopool <- function(p) pnorm(sum(qnorm(p, lower.tail = FALSE))/ sqrt(length(p)), lower.tail = FALSE)
estimatePrb(stopool, 0.05, M = 10, interval = c(.Machine$double.eps, 1))
estimatePrb(stopool, 0.05, M = 10, b = 0.5, interval = c(.Machine$double.eps, 1))
}
\author{
Chris Salahub
}
