% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parametricExampleOptim.R
\name{estimateTruncParamModel}
\alias{estimateTruncParamModel}
\title{A function for estimating parmetric model 1 (truncation)}
\usage{
estimateTruncParamModel(x, threshold, sigmaGrid = NULL, maxSigma = NULL,
  nIntPoints = 10^4, seed = NULL, knownSigma = NULL)
}
\arguments{
\item{x}{the observed z-scores}

\item{threshold}{the selection threshold for the selection rule
\code{x < threshold[1]} or \code{x > threshold[2]}. If a scalar is provided,
then threshold is set to \code{c(-abs(threshold), abs(threshold))}.}

\item{sigmaGrid}{an optional grid of stadard deviation values over which
to evaluate the likelihood}

\item{maxSigma}{an optional maximal value for the standard deviation of mu}

\item{nIntPoints}{number of samples to take for numerical integration}

\item{seed}{an optional seed}

\item{knownSigma}{the standard deviation of the distribution of the
normal means is known if it is known (for computing the true bayes rule)}
}
\description{
A function for estimating the model which assumes that
mu ~ N(0, b) and y ~ TN(mu, 1, threshold). The selection rule is
\code{x < threshold[1]} or \code{x > threshold[2]}.
}
