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

paretoML <- function(x, location = NULL, weights = NULL, barrier = 0) {
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

n <- 10^4
threshold <- 2
mprob <- c(0.8, 0.2, 0)
muNorm <- rnorm(n, mean = 0, sd = 0.25)
muExp <- rnorm(n, mean = 0, sd = 2) #rexp(n, rate = 0.5) * (1 - 2 * rbinom(n, 1, 0.5))
muGamma <- rgamma(n, 1, 3)
memberships <- t(rmultinom(n, 1, mprob))
mu <- rowSums(memberships * cbind(muNorm, muExp, muGamma))
z <- mu + rnorm(n)

keep <- z > threshold
x <- z[keep]

fit <- paretoML(x, location = threshold, barrier = 0.005)
scale <- fit$par[1]
shape <- fit$par[2]
par(mfrow = c(1, 1))
plot(density(x), ylim = c(0, 1))
lines(sort(x), paretoDens(sort(x), threshold, scale, shape, FALSE), type = "l", col = "blue")

library(gPdtest)
gpdfit <- gpd.fit(x, "amle")
lines(sort(x), paretoDens(sort(x), threshold, gpdfit[2], gpdfit[1], FALSE), type = "l", col = "red")

mixfit <- normParetoMixEM(x, threshold, emIterations = 200, barrier = 0)
mixdens <- normParetoMixDens(x, mixfit)
lines(sort(x), normParetoMixDens(sort(x), mixfit), col = "green")


# Tweedy -------
npTweed <- normParetoTweedy(mixfit)
location <- threshold
scale <- fit$par[1]
shape <- fit$par[2]
tweedy <- paretoTweedy(x, location, scale, shape)
mean(tweedy - mu[keep])
mean(x - mu[keep])
mean(npTweed - mu[keep])










