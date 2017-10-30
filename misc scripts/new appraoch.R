library(ggplot2)
library(truncnorm)
library(ggplot2)
library(truncnorm)
n <- 200
samp <- numeric(n)
thetas <- numeric(n)
current <- 1
threshold <- c(-2, 2)
tries <- 0
while(current <= n) {
  tries <- tries + 1
  theta <- rexp(1, 0.5) * (1 - 2*rbinom(1, 1, 0.5))
  x <- rnorm(1, theta)
  if(x > threshold[2]) {
    samp[current] <- x
    thetas[current] <- theta
    current <- current + 1
  }
}
plot(thetas, samp)
abline(a = 0, b = 1)
data <- data.frame(samp = samp, theta = thetas)
ggplot(data) + geom_density(aes(x = thetas), col = "black")  +
  geom_density(aes(x = samp), col = "red", linetype = 2) +
  theme_bw() + xlab("values") +
  geom_vline(xintercept = threshold, col = "black", linetype = 3)

# Helper functions --------------
dtrunc <- function(z, thetas, thershold, log = FALSE) {
  probs <- log(pnorm(threshold[1], thetas) + pnorm(threshold[2], thetas, lower.tail = FALSE))
  densities <- dnorm(z, thetas, log = TRUE)
  if(log) {
    return(densities - probs)
  } else{
    return(exp(densities - probs))
  }
}
compExp <- function(coefs, support, thetas) {
  df <- length(coefs)
  supp <- sapply(1:df, function(i) support^i)
  dens <- supp %*% coefs
  dens <- exp(dens - max(dens))
  suffExp <- as.numeric(t(supp) %*% dens / sum(dens))
  return(suffExp)
}
optimF <- function(coefs, support, thetas) {
  df <- length(coefs)
  dens <- sapply(1:df, function(i) support^i) %*% coefs
  maxdens <- max(dens)
  dens <- exp(dens - maxdens)
  Z <- sum(dens)
  thetadens <- (sapply(1:df, function(i) thetas^i) %*% coefs - maxdens) - log(Z)
  return(-mean(thetadens))
}
optimGrad <- function(coefs, support, thetas) {
  suffStat <- colMeans(sapply(1:df, function(i) thetas^i))
  suffExp <- compExp(coefs, support, thetas)
  # print(suffExp - suffStat)
  return(suffExp - suffStat)
}
compCDF <- function(coefs, support) {
  dens <- as.numeric(sapply(1:length(coefs), function(i) support^i) %*% coefs)
  dens <- exp(dens - max(dens))
  Z <- sum(dens)
  cdf <- cumsum(dens) / Z
  return(cdf)
}
computeThetaDens <- function(coefs, support) {
  df <- length(coefs)
  supp <- sapply(1:df, function(i) support^i)
  dens <- supp %*% coefs
  dens <- exp(dens - max(dens))
  dens <- dens / sum(dens)
  return(dens)
}
sampCondTheta <- function(z, coefs, support, thetadens, threshold = 2) {
  zdens <- dtrunc(z, support, threshold)
  zdens[is.infinite(zdens)] <- 0
  dens <- zdens * thetadens
  cdf <- cumsum(dens) / sum(dens)
  return(support[min(which(runif(1) < cdf))])
}

# Loading science data ------------------
scidat <- readRDS("dataAnalysis/sciDat.rds")
samp <- scidat$orig
thetas <- scidat$rep

# Analysis ----------------------
df <- 5
support <- seq(from = -max(abs(samp)), to = max(abs(samp)), by = 0.1)
init <- samp
coefs <- rep(0, df)
fit <- optim(coefs, fn = optimF, gr = optimGrad, method = "BFGS",
             support = support, thetas = samp)
coefs <- fit$par
compExp(coefs, support, thetas)

# iterating
delay <- 1
iterations <- 2000
estimates <- rep(0, length(samp))
errors <- matrix(nrow = iterations, ncol = 2)
for(i in 1:iterations) {
  thetadens <- as.numeric(computeThetaDens(coefs, support))
  current <- sapply(samp, sampCondTheta, coefs, support, thetadens, threshold)
  if(TRUE) {
    thetasamp <- current
  } else {
    thetasamp <- c(thetasamp, current)
  }

  if(i >= delay) {
    w <- 1 / max(i - delay, 1)
    estimates <- (1 - w) * estimates + current * w
  }
  fit <- optim(coefs, fn = optimF, gr = optimGrad, method = "BFGS",
               support = support, thetas = thetasamp)
  coefs <- fit$par
  # plot(rep(samp, max(i - delay + 1, 1)), thetasamp)
  # abline(a = 0, b = 1)
  errors[i, 1] <- mean(estimates - thetas, na.rm = TRUE)
  errors[i, 2] <- sqrt(mean((estimates - thetas)^2, na.rm = TRUE))
  print(round(c(i, coefs), 3))
}

# Evaluating estimation ------------
plot(support, thetadens, type = 'l')
lines(support, 0.5 * dexp(abs(support)))
abline(a = 0, b = 1)

plot(sort(samp), sort(estimates), ylim = c(-3, 9), col = "blue", type = "l",
     lty = 2, lwd = 3)
points(samp, thetas, col = "red")
abline(a = 0, b = 1)
cond <- sapply(samp, computeConditionalEst, threshold)
mean(samp - thetas, na.rm = TRUE)
mean(estimates - thetas, na.rm = TRUE)
mean(cond - thetas, na.rm = TRUE)

sqrt(mean((samp - thetas)^2, na.rm = TRUE))
sqrt(mean((estimates - thetas)^2, na.rm = TRUE))
sqrt(mean((cond - thetas)^2, na.rm = TRUE))


forplot <- rbind(data.frame(est = samp, type = "observed"),
                 data.frame(est = na.omit(thetas), type = "rep"),
                 data.frame(est = estimates, type = "tweedy"),
                 data.frame(est = cond, type = "conditional"))
if(estFull) {
  forplot <- rbind(forplot, data.frame(est = c(samp, censored), type = "sim"))
}
ggplot(forplot) + geom_density(aes(x = est, col = type, linetype = type)) + theme_bw()

# plot(support, thetadens, type = "l")

plot(1:nrow(errors), errors[, 1], lty = 1, col = "red", type = "l", ylim = c(0, 3))
lines(1:nrow(errors), errors[, 2], lty = 1, col = "blue", type = "l")

