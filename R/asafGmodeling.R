#' Computes Empirical Bayes Estimates via Deconvolution for Truncated Observations
#'
#' This function has a similar functionality as \code{\link[deconvolveR]{deconv}}, excpet
#' that it is meant to handle the case in which only the selected are observed.
#'
#' @param x the observed (truncated) sample of z-scores
#'
#' @param threshold the threshold used to screen the observations
#'
#' @param meanValues a grid of discrete mean parameter values
#'
#' @param twoSided whether the selection was one-sided or two-sided.
#' if TRUE, then a selection rule of abs(x) > abs(threshold) is assumed. Otherwise,
#' a selection rule of x > threshold is assumed.
#'
#' @param splineDegree the number of degrees for the spline fit in the deconvolution
#'
#' @param binWidth the function discretizes the sample into bins, this is the bin width.
#' Take precedence over nBins.
#'
#' @param nBins number of bins to use in discretization
#'
#' @seealso \code{\link{predict.truncDeconv}}
#'
#' @importFrom deconvolveR deconv
#' @importFrom splines ns
#' @export
truncDeconv <- function(x, threshold, meanValues = NULL ,twoSided = TRUE,
                        splineDegree = 10,
                        binWidth = 0.1,
                        nBins = 100) {
  # Some preliminaries
  signs <- sign(x)
  if(twoSided) {
    x <- abs(x)
    threshold <- abs(threshold)
  }

  # Defining Bins
  if(is.null(binWidth)) {
    if(is.null(nBins)) {
      stop("binWidth or nBins must be specified!")
    } else {
      if(nBins < 5) {
        stop("nBins must be larger than 5!")
      }
      binWidth <- (max(x) - threshold) / (nBins - 2)
    }
  }
  maxVal <- max(x)
  breaks <- seq(from = threshold - binWidth/2, to = maxVal + binWidth/2, by = binWidth)
  breaks[length(breaks)] <- max(breaks[length(breaks)], maxVal + binWidth / 2)

  # Setting up spline basis
  if(is.null(meanValues)) {
    meanValues <- seq(from = -maxVal/2, to = maxVal, by = 0.1)
  }
  Q <- ns(meanValues, df = splineDegree)

  # Computing histogram and P matrix
  marginalHistogram <- hist(x, breaks = breaks, plot = FALSE)
  y <- marginalHistogram$counts
  P <- matrix(nrow = length(y), ncol = length(meanValues))
  for(j in 1:ncol(P)) {
    for(i in 1:nrow(P)) {
      P[i, j] <- pnorm(breaks[i + 1], meanValues[j], 1) - pnorm(breaks[i], meanValues[j], 1)
    }
    selectionProb <- pnorm(threshold, meanValues[j], 1, lower.tail = FALSE)
    if(twoSided) {
      selectionProb <- selectionProb + pnorm(-threshold, meanValues[j], 1)
    }
    P[, j] <- P[, j] / selectionProb
  }

  # Estimating prior using the deconvolveR package
  truncDeconv <- deconv(tau = meanValues, P = P, Q = Q, y = y, pDegree = deg)
  class(truncDeconv) <- "truncDeconv"
  truncDeconv$x <- x * signs
  truncDeconv$binCounts <- y
  truncDeconv$binCenters <- marginalHistogram$mids
  return(truncDeconv)
}


deconvCompCondExp <- function(index, deconvFit) {
  sum(deconvFit$stats[, 'theta'] * deconvFit$stats[, "g"] * deconvFit$P[index, ]) / sum(deconvFit$stats[, "g"] * deconvFit$P[index, ])
}

#' Computes the Bayes Rule based on a \code{truncDeconv} model fit
#'
#' Computes the Bayes rule for \code{\link{truncDeconv}} model fit.
#'
#' @param object an object of class \code{\link{truncDeconv}}
#'
#' @param binned whether to return the binned Bayes rule, or the linear interpolation
#' for the data/newX
#'
#' @param newX an optional new dataset to compute the Empirical Bayes Rule for
#'
#' @export
predict.truncDeconv <- function(object, binned = FALSE, newX = NULL, ...) {
  if(is.null(newX)) {
    newX <- object$x
  }

  binnedBayes <- sapply(1:nrow(object$P), FUN = deconvCompCondExp, deconvFit = object)
  if(is.null(newX) | binned) {
    result <- data.frame(x = object$binCenters, bayes = binnedBayes)
  } else {
    signs <- sign(newX)
    newX <- abs(newX)
    centers <- object$binCenters
    if(max(newX) > max(centers)) {
      centers <- c(centers, max(newX))
      binnedBayes <- c(binnedBayes, max(newX))
    }
    bayes <- approx(x = centers, y = binnedBayes, xout = newX)$y
    result <- data.frame(x = newX * signs, bayes = bayes * signs)
  }

  return(result)
}

#' A Function for Computing the Bayes Rule for a \code{deconv} model fit
#'
#' Computes the Bayes Rule, as estimated using the  function,
#' with the option \code{family = "Normal"}.
#'
#' @param x data for which to compute the Empirical Bayes Rule
#'
#' @param object a \code{\link[deconvolveR]{deconv}} model fit
#'
#' @export
deconvTweed <- function(x, object) {
  result <- numeric(length(x))
  theta <- object$stats[, "theta"]
  gdens <- object$stats[, "g"]
  for(i in 1:length(x)) {
      numerator <-  dnorm(x[i], mean = theta) * gdens
      denominator <- sum(numerator)
      numerator <- sum(numerator * theta)
      result[i] <- numerator / denominator
  }
  return(result)
}

