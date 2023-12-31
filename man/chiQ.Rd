% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chi.R
\name{chiQ}
\alias{chiQ}
\title{Chi-squared centrality quotient}
\usage{
chiQ(kappa, M, alpha = 0.05)
}
\arguments{
\item{kappa}{numeric between 0 and infinity}

\item{M}{integer sample size greater than 0}

\item{alpha}{numeric between 0 and 1}
}
\value{
A numeric between 0 and 1.
}
\description{
Computes the centrality quotient of the chi-square
pooled p-value.
}
\details{
The centrality quotient of a pooled p-value measures the
relative preference it gives to p-values all sharing the same
level of evidence over a single test with strong evidence relative
to others. For the chi-square pooled p-value, this is a
conditional probability which this function computes.
}
\examples{
chiQ(2, 10, 0.05)
chiQ(2, 20, 0.05)
chiQ(0.5, 20, 0.05)
}
\author{
Chris Salahub
}
