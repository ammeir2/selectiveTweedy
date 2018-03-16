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
  P <- matrix(nrow = length(y), ncol = length(tau))
  for(j in 1:ncol(P)) {
    for(i in 1:nrow(P)) {
      P[i, j] <- pnorm(breaks[i + 1], tau[j], 1) - pnorm(breaks[i], tau[j], 1)
    }
    selectionProb <- pnorm(threshold, tau[j], 1, lower.tail = FALSE)
    if(twoSided) {
      selectionProb <- selectionProb + pnorm(-threshold, tau[j], 1)
    }
    P[, j] <- P[, j] / selectionProb
  }

  # Estimating prior using the deconvolveR package
  truncDeconv <- deconv(tau = tau, P = P, Q = Q, y = y, pDegree = deg)
  class(truncDeconv) <- "truncDeconv"
  truncDeconv$x <- x * signs
  truncDeconv$binCounts <- y
  truncDeconv$binCenters <- marginalHistogram$mids
  return(truncDeconv)
}


deconvCompCondExp <- function(index, deconvFit) sum(deconvFit$stats[, 'theta'] * priorEst * deconvFit$P[index, ]) / sum(priorEst * deconvFit$P[index, ])

#' @export
predict.truncDeconv <- function(object, binned = FALSE, newX = NULL, ...) {
  if(is.null(newX)) {
    newX <- object$x
  }

  binnedBayes <- sapply(1:nrow(object$P), FUN = deconvCompCondExp, object)
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

