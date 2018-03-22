truncModelIntegratedLogDens <- function(sigma, threshold, x, intGrid) {
  muDens <- dnorm(intGrid, sd = sigma) / (pnorm(-threshold, mean = intGrid) + pnorm(threshold, mean = intGrid, lower.tail = FALSE))
  xDens <- outer(x, intGrid, FUN = function(x, m) dnorm(x, m))
  logDens <- (xDens %*% muDens) %>% log() %>% sum()
  return(logDens)
}


#' A function for estimating parmetric model 1 (truncation)
#'
#' A function for estimating the model which assumes that
#' mu ~ N(0, b) and y ~ TN(mu, 1, threshold).
#'
#' @param x the observed z-scores
#' @param threshold the selection threshold for the selection rule abs(x) > threshold
#' @param maxSigma an optional maximal value for the standard deviation of mu
#' @param nIntPoints number of samples to take for numerical integration
#' @param seed an optional seed
#'
#' @export
estimateTruncParamModel <- function(x, threshold, maxSigma = NULL,
                                    nIntPoints = 10^4, seed = NULL) {
  if(!is.null(seed)) {
    set.seed(seed)
  }

  if(is.null(maxSigma)) {
    maxSigma <- (max(abs(x)) - 1)
  }

  if(maxSigma <= 0) {
    stop("maxSigma must be a number larger than zero!")
  }

  # estimating the standard deviation of the normal prior
  integrationGrid <- rnorm(nIntPoints, sd = maxSigma)
  fit <- optimize(f = truncModelIntegratedLogDens,
                  interval = c(0, maxSigma), maximum = TRUE,
                  threshold = threshold, x = x,
                  intGrid = integrationGrid)
  estSD <- fit$maximum

  # computing empirical bayes estiamates for mu
  bayesIntGrid <- rnorm(nIntPoints, sd = estSD)
  normDens <- outer(x, bayesIntGrid, function(x, m) dnorm(x, m))
  selectProb <- 1 / (pnorm(-threshold, bayesIntGrid) + 1 - pnorm(threshold, bayesIntGrid))
  denom <- (normDens %*% selectProb) %>% as.numeric()
  numIntegrator <- bayesIntGrid / selectProb
  numerator <- (normDens %*% numIntegrator) %>% as.numeric()
  bayesRule <- numerator / denom

  result <- list(x = x, bayesRule = bayesRule, estSD = fit$maximum)
  return(result)
}

censoredTruncNormDens <- function(sigma, x, threshold) {
  normDens <- dnorm(x, sd = sigma, log = TRUE)
  sProb <- pnorm(-threshold, sd = sigma) + 1 - pnorm(threshold, sd = sigma)
  return(sum(normDens - log(sProb)))
}

#' A function for estimating parametric model 2 (censoring)
#'
#' Parametric model 2 assumes that mu ~ N(0, b) and x ~ N(mu, 1), where x is discarded
#' if abs(x) < threshold.
#'
#' @param x the observed data z-scores
#' @param threshold the selection threshold for the selection rule abs(x) > threshold
#'
#' @export
estimateCensoredParamModel <- function(x, threshold) {
  maxSigma <- max(abs(x))
  estSD <- optimize(censoredTruncNormDens, interval = c(1, maxSigma),
                    maximum = TRUE,
                    x = x,
                    threshold = threshold)$maximum
  bayes <- x * (1 - 1 / estSD^2)
  estSD <- sqrt(estSD^2 - 1)
  result <- list(x = x, bayes = bayes, estSD = estSD)
  return(result)
}

truncDnormForMLE <- function(mu, x, sd, threshold) {
  normDens <- dnorm(x, mu, sd, TRUE)
  sProb <- pnorm(-threshold, mu, sd) + 1 - pnorm(threshold, mu, sd)
  return(normDens - log(sProb))
}

#' Compute the univariate conditional MLE for truncated normal data
#'
#' @param x the observed truncated normal samples
#' @param threshold the threshold for the selection rule abs(x) > threshold
#' @param xsd the standard deviation of the normal observations
#'
#' @export
univTruncNormMLE <- function(x, threshold, xsd = 1) {
  if(!(length(xsd) %in% c(1, length(x)))) {
    stop("xsd must be of either, length 1 or length(x)!")
  }

  if(length(xsd) == 1) {
    xsd <- rep(xsd, length(x))
  }

  if(any(abs(x) < threshold)) {
    stop("Some observations are below the threshold!")
  }

  if(threshold <= 0) {
    stop("threshold must be >= 0!")
  }

  mle <- numeric(length(x))
  for(i in 1:length(mle)) {
    interval <- sort(c(0, x))
    mle[i] <- optimize(truncDnormForMLE, interval = interval, maximum = TRUE,
                       x = x[i], sd = xsd[i], threshold = threshold)$maximum
  }
  return(mle)
}
