% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/asafGmodeling.R
\name{deconvTweed}
\alias{deconvTweed}
\title{A Function for Computing the Bayes Rule for a \code{deconv} model fit}
\usage{
deconvTweed(x, object)
}
\arguments{
\item{x}{data for which to compute the Empirical Bayes Rule}

\item{object}{a \code{\link[deconvolveR]{deconv}} model fit}
}
\description{
Computes the Bayes Rule, as estimated using the  function,
with the option \code{family = "Normal"}.
}
