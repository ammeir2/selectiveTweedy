#' @importFrom truncnorm dtruncnorm
tnDens <- function(par, x, threshold) {
  mu <- par[1]
  sd <- par[2]
  pneg <- pnorm(-abs(threshold), mu, sd, log.p = TRUE)
  ppos <- pnorm(abs(threshold), mu, sd, lower.tail = FALSE, log.p = TRUE)
  ppos <- 1 / (1 + exp(pneg - ppos))
  pneg <- 1 - ppos
  dens <- numeric(length(x))
  dens[x > 0] <- ppos * dtruncnorm(x[x > 0], abs(threshold), Inf, mu, sd)
  dens[x < 0] <- pneg * dtruncnorm(x[x < 0], -Inf, -abs(threshold), mu, sd)
  return(dens)
}

tnLogLik <- function(params, x, threshold, weights = NULL, sdLbound) {
  params[2] <- max(params[2], sdLbound)
  dens <- tnDens(params, x, threshold) %>% log() %>% weighted.mean(weights)
  return(-dens)
}

truncNormML <- function(weights, x, threshold, sdLbound = 1) {
  params <- c(mu = mean(x), sd = sd(x))
  fit <- optim(par = params, fn = tnLogLik,
               x = x, threshold = threshold, weights = weights,
               sdLbound = sdLbound)
  fit$par[2] <- max(fit$par[2], sdLbound)
  return(fit$par)
}

mixDens <- function(x, threshold, normFit, paretoFit, probs) {
  normDens <- apply(normFit, 2, tnDens, x, threshold) %*% probs[2:length(probs)] %>% as.numeric()
  paretoDens <- paretoDens(x = x, location = threshold, scale = paretoFit$par[1], shape = paretoFit$par[2], log = FALSE) * probs[1]
  return(normDens + paretoDens)
}

#' A function for fitting a Genearlized Pareto/Truncated Normal mixture to a data
#'
#' The aim of this function is to estimate the marginal distribution of a normal, convolved
#' with an unknown prior and then truncated. The model fit can be used to compute the
#' Empirical Bayes rule using Tweedy's Formula.
#'
#' @param x the observed (truncated) z-scores
#'
#' @param threshold the threshold used for screening, the selection rule is abs(x) > abs(threshold)
#'
#' @param normComps number of truncated normal components to fit
#'
#' @param iterations number of EM iterations
#'
#' @param paretoComp whether to include a Generalzed Pareto mixture component
#'
#' @param verbose whether to print a progress bar
#'
#' @seealso \code{\link{predict.ptnMix}}
#'
#' @importFrom mixtools normalmixEM
#' @import magrittr
#' @export
paretoTruncNormMix <- function(x, threshold, normComps = 1, iterations = 100, paretoComp = FALSE,
                               verbose = TRUE) {
  # Fitting mixture of truncated normal and pareto
  paretoFit <- paretoML(x, location = threshold, barrier = 0.005)

  if(normComps > 1) {
    initClusters <- suppressMessages(normalmixEM(x, k = normComps))
  } else {
    initClusters <- list()
    initClusters$posteriors <- matrix(1, ncol = 1, nrow = length(x))
  }
  normFit <- apply(initClusters$posterior, 2, truncNormML, x, threshold, sdLbound = 1)

  probs <- rep(1 / (1 + normComps), 1 + normComps)
  if(!paretoComp) {
    probs[1] <- 0
    probs <- probs / sum(probs)
    paretoPost <- rep(1, length(x))
  }

  if(verbose) pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  for(i in 1:iterations) {
    # E-Step --------
    if(paretoComp) {
      paretoPost <- paretoDens(x, threshold, paretoFit$par[1], paretoFit$par[2], log = FALSE) * probs[1]
    }
    tnPost <- apply(normFit, 2, tnDens, x = x, threshold = threshold)
    tnPost <- (t(tnPost) * probs[2:(normComps + 1)]) %>% t()
    totDens <- rowSums(tnPost) + paretoPost
    tnPost <- tnPost / totDens
    paretoPost <- paretoPost / totDens

    # M-Step ------
    if(paretoComp) {
      paretoFit <- paretoML(x, location = threshold, weights = paretoPost, barrier = 0.005)
    }
    normFit <- apply(tnPost, 2, truncNormML, x = x, threshold = threshold, sdLbound = 1)
    probs[1] <- mean(paretoPost)
    probs[2:length(probs)] <- colMeans(tnPost)

    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) close(pb)

  posteriors <- cbind(paretoPost, tnPost)
  result <- list(x = x, threshold = threshold,
                 posteriors = posteriors,
                 probs = probs,
                 paretoFit = paretoFit,
                 normFit = normFit)
  class(result) <- "ptnMix"
  return(result)
}

#' Applies Tweedy Correction based on a Pareto/Truncated Normal Mixture
#'
#' Computes the Tweedy correction for a truncated sample based on a
#' \code{\link{paretoTruncNormMix}} model fit.
#'
#' @param object an object of class \code{ptnMix}, obtained from fitting a
#' \code{\link{paretoTruncNormMix}} model
#'
#' @export
predict.ptnMix <- function(object, ...) {
  # browser()
  x <- object$x
  threshold <- object$threshold
  paretoFit <- object$paretoFit
  normFit <- object$normFit
  mu <- normFit[1, ]
  sd <- normFit[2, ]
  posteriors <- object$posteriors

  paretoTweed <- paretoTweedy(x, location = threshold,
                              scale = paretoFit$par[1],
                              shape = paretoFit$par[2])
  paretoTweed <- x - paretoTweed
  paretoTweed <- paretoTweed * posteriors[, 1]

  normTweed <- numeric(length(x))
  nComps <- ncol(object$posteriors)
  for(i in 1:length(x)) {
    normTweed[i] <- sum(posteriors[i, -1] * (x[i] - mu) / sd^2)
  }

  tweed <- x - normTweed - paretoTweed
  return(tweed)
}





