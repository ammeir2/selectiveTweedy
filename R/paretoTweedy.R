#' The density of the generalized pareto distribution
paretoDens <- function(x, location, scale, shape, log = TRUE) {
  if(any(x < location)) {
    stop("Some observations are below the lower bound (location parameter)!")
  }
  logdens <- - log(scale) - (1 / shape + 1) * log(1 + shape * (x - location) / scale)
  if(!log) {
    return(exp(logdens))
  } else {
    return(logdens)
  }
}

#' Computes the log-likelihood of the generalized pareto distribution
paretoLogLik <- function(params, x, location = NULL, weights = NULL, barrier = 0) {
  if(is.null(weights)) {
    weights <- rep(1, length(x))
  }

  if(is.null(location)) {
    location <- min(x)
  }

  loglik <- paretoDens(x, location, params[1], params[2], TRUE)
  loglik <- sum(weights * loglik)
  penalty <- - log(params[2]) * barrier * sum(weights)
  # print(c(params, -loglik))
  return(-loglik + penalty)
}

#' Computes the gradient of the minus log-likelihood of the generalized pareto distribution
paretoGrad <- function(params, x, location = NULL, weights = NULL, barrier = 0.005) {
  if(is.null(weights)) {
    weights <- rep(1, length(x))
  }

  if(is.null(location)) {
    location <- min(x)
  }

  scale <- params[1]
  shape <- params[2]

  z <- (x - location) / scale
  dscale <- -1/scale + z / (scale^2 + scale^2 * shape * z) * (1/shape + 1)
  dshape <- 1/shape^2 * log(1 + shape * z) - (1/shape + 1) * z / (1 + shape * z)

  dscale <- sum(weights * dscale)
  dshape <- sum(weights * dshape)
  penalty <- - sum(weights) / shape * barrier
  return(-c(dscale + penalty, dshape))
}

#' Computes the MLE of a generalized pareto model
paretoML <- function(x, location = NULL, weights = NULL, barrier = 0.005) {
  if(is.null(weights)) {
    weights <- rep(1, length(x))
  }

  if(is.null(location)) {
    location <- min(x)
  }

  init <- c(sd(x), 1)
  fit <- optim(init, fn = paretoLogLik, gr = paretoGrad,
               x = x, location = location, weights = weights,
               barrier = barrier,
               lower = 10^-12,
               method = "L-BFGS-B")
  return(fit)
}

#' Computes the tweedy correction based on an estimated generalized pareto distribution
paretoTweedy <- function(x, location, scale, shape) {
  adjust <- x - shape/scale * (1/shape + 1) / (1 + shape * (x - location) / scale)
  return(adjust)
}


#' Tweedy correction for tails of distributions based on the generalized pareto distribution
#'
#' @export
paretoTweedyCorrection <- function(x, threshold, weights = NULL, barrier = 0.005) {
  xsign <- sign(x)
  x <- x * xsign

  mlfit <- paretoML(x, location = threshold, weights = weights, barrier = barrier)
  correction <- paretoTweedy(x, threshold, mlfit$par[1], mlfit$par[2])

  result <- list(x = x,
                 estimate = correction,
                 mlfit = mlfit)

  return(result)
}


truncLogNormDens <- function(mu, x, sd, threshold) {
  dens <- dnorm(x, mu, sd, log = TRUE)
  prob <- pnorm(-threshold, mu, sd) + pnorm(threshold, mu, sd, lower.tail = FALSE)
  return(dens - log(prob))
}

#' Computes the conditioal MLE for a truncated normal observation
#'
#' @export
truncNormMLE <- function(x, threshold, sd = 1) {
  opt <- optimize(f = truncLogNormDens, lower = 0, upper = abs(x),
                  x = abs(x), threshold = abs(threshold), sd = sd,
                  maximum = TRUE)
  return(opt$maximum * sign(x))
}

