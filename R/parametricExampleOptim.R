#' A function for estimating parmetric model 1 (truncation)
#'
#' A function for estimating the model which assumes that
#' mu ~ N(0, b) and y ~ TN(mu, 1, threshold). The selection rule is
#' \code{x < threshold[1]} or \code{x > threshold[2]}.
#'
#' @param x the observed z-scores
#' @param threshold the selection threshold for the selection rule
#' \code{x < threshold[1]} or \code{x > threshold[2]}. If a scalar is provided,
#' then threshold is set to \code{c(-abs(threshold), abs(threshold))}.
#'
#' @param sigmaGrid an optional grid of stadard deviation values over which
#' to evaluate the likelihood
#' @param maxSigma an optional maximal value for the standard deviation of mu
#' @param nIntPoints number of samples to take for numerical integration
#' @param seed an optional seed
#' @param knownSigma the standard deviation of the distribution of the
#' normal means is known if it is known (for computing the true bayes rule)
#'
#' @export
estimateTruncParamModel <- function(x, threshold, sigmaGrid = NULL,
                                    maxSigma = NULL,
                                    nIntPoints = 10^4, seed = NULL,
                                    knownSigma = NULL) {
  threshold <- checkThreshold(threshold)

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
  integrationGrid <- sort(rnorm(nIntPoints, sd = 1))
  if(is.null(knownSigma)) {
    if(is.null(sigmaGrid)) {
      fit <- optimize(f = truncModelIntegratedLogDens,
                      interval = c(0, maxSigma), maximum = TRUE,
                      threshold = threshold, x = x,
                      intGrid = integrationGrid)
      estSD <- fit$maximum
    } else {
      likelihoods <- sapply(sigmaGrid, truncModelIntegratedLogDens,
                            threshold = threshold,
                            x = x,
                            intGrid = integrationGrid)
      estSD <- sigmaGrid[which.max(likelihoods)]
    }
  } else if(knownSigma <= 0) {
    stop("knownSigma must be either larger than zero, or NULL (if it is to be estimated)")
  } else {
    estSD <- knownSigma
  }

  # computing empirical bayes estiamates for mu
  bayesRule <- truncModelIntegratedLogDens(estSD, threshold, x,
                                           integrationGrid,
                                           expectation = TRUE)

  result <- list(x = x, bayesRule = bayesRule, estSD = estSD)
  return(result)
}

truncModelIntegratedLogDens <- function(sigma, threshold, x, intGrid,
                                        expectation = FALSE) {
  importanceCoef <- 1.5
  intGrid <- intGrid * sigma * importanceCoef
  ppos <- pnorm(threshold[2], mean = intGrid, lower.tail = FALSE)
  pneg <- pnorm(threshold[1], mean = intGrid)
  # if(expectation) {
  #   ppos <- pmax(ppos, 10^-2)
  # }
  muDens <- numeric(length(intGrid))
  for(i in 1:length(muDens)) {
    if(max(ppos[i], pneg[i]) < 10^-5) {
      dens <- dnorm(intGrid[i], sd = sigma, log = TRUE) - dnorm(intGrid[i], sd = sigma * importanceCoef, log = TRUE)
      pdens <- dens - pnorm(threshold[1], mean = intGrid[i], log.p = TRUE)
      ndens <- dens - pnorm(threshold[2], mean = intGrid[i], log.p = TRUE, lower.tail = FALSE)
      muDens[i] <- min(pdens, ndens) %>% exp()
    } else {
      dens <- exp(dnorm(intGrid[i], sd = sigma, log = TRUE) - dnorm(intGrid[i], sd = sigma * importanceCoef, log = TRUE))
      muDens[i] <- dens / #(dnorm(intGrid[i], sd = sigma)) /
        (pnorm(threshold[1], mean = intGrid[i]) +
           pnorm(threshold[2], mean = intGrid[i], lower.tail = FALSE))
    }
  }

  xDens <- outer(x, intGrid, FUN = function(x, m) dnorm(x, m))
  if(expectation) {
    numerator <- outer(x, intGrid, FUN = function(x, m) m * dnorm(x, m))
    numerator <- (numerator %*% muDens / length(intGrid)) %>% as.numeric()
    denominator <- (xDens %*% muDens / length(intGrid)) %>% as.numeric()
    bayes <- numerator / denominator
    return(bayes)
  } else {
    logDens <- (xDens %*% muDens) %>% log() %>% sum()
    return(logDens)
  }
}

checkThreshold <- function(threshold) {
  if(length(threshold) == 1) {
    threshold <- c(-abs(threshold), abs(threshold))
  } else if(length(threshold) == 2) {
    threshold <- sort(threshold)
  } else {
    stop("threshold must be either a scalar or a vector of length 2.")
  }
  return(threshold)
}


#' A function for estimating parametric model 2 (censoring)
#'
#' Parametric model 2 assumes that mu ~ N(0, b) and x ~ N(mu, 1),
#' where x is discarded
#' if \code{x > threshold[1]} and \code{x < threshold[2]}
#'
#' @param x the observed data z-scores
#' @param threshold the selection threshold for the selection rule
#' \code{x < threshold[1]} or \code{x > threshold[2]}. If a scalar is provided,
#' then threshold is set to \code{c(-abs(threshold), abs(threshold))}.
#'
#' @export
estimateCensoredParamModel <- function(x, threshold) {
  threshold <- checkThreshold(threshold)

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

censoredTruncNormDens <- function(sigma, x, threshold) {
  normDens <- dnorm(x, sd = sigma, log = TRUE)
  pprob <- 1 - pnorm(threshold[2], sd = sigma)
  nprob <- pnorm(threshold[1], sd = sigma)
  if(max(pprob, nprob) > 10^-4) {
    return(sum(normDens - log(nprob + pprob)))
  } else {
    pprob <- pnorm(threshold[2], sd = sigma, log.p = TRUE, lower.tail = FALSE)
    nprob <- pnorm(threshold[1], sd = sigma, log.p = TRUE)
    return(sum(normDens - max(pprob, nprob)))
  }
}


#' Compute the univariate conditional MLE for truncated normal data
#'
#' @param x the observed truncated normal samples
#' @param threshold the selection threshold for the selection rule
#' \code{x < threshold[1]} or \code{x > threshold[2]}. If a scalar is provided,
#' then threshold is set to \code{c(-abs(threshold), abs(threshold))}.
#' @param xsd the standard deviation of the normal observations
#'
#' @export
univTruncNormMLE <- function(x, threshold, xsd = 1) {
  threshold <- checkThreshold(threshold)

  if(!(length(xsd) %in% c(1, length(x)))) {
    stop("xsd must be of either, length 1 or length(x)!")
  }

  if(length(xsd) == 1) {
    xsd <- rep(xsd, length(x))
  }

  if(any(x > threshold[1] & x < threshold[2])) {
    stop("Some observations are below the threshold!")
  }

  mle <- numeric(length(x))
  for(i in 1:length(mle)) {
    interval <- sort(c(0, x))
    mle[i] <- nlm(f = truncDnormForMLE, p = x[i],
                  # interval = interval, maximum = TRUE,
                  x = x[i], sd = xsd[i], threshold = threshold)$estimate
  }
  return(mle)
}

truncDnormForMLE <- function(mu, x, sd, threshold) {
  normDens <- dnorm(x, mu, sd, TRUE)
  ppos <- pnorm(threshold[2], mean = mu, sd = sd, log.p = TRUE, lower.tail = FALSE)
  pneg <- pnorm(threshold[1], mean = mu, sd = sd, log.p = TRUE)
  if(max(exp(c(ppos, pneg))) < 10^-4) {
    logdens <- -normDens + max(ppos, pneg)
  } else {
    logdens <- -normDens + log(exp(ppos) + exp(pneg))
  }
  # print(c(x, mu, logdens))
  return(logdens)
}

