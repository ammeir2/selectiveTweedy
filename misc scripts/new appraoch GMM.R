# Loading Science Data --------------------

# Artificial Data ------------------------
library(ggplot2)
library(truncnorm)
library(mixtools)
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

# Loading science data ------------------
scidat <- readRDS("dataAnalysis/sciDat.rds")
samp <- scidat$orig
thetas <- scidat$rep

# Defining helper functions ---------------
computeGmmDens <- function(thetas, df, ...) {
  fit <- normalmixEM(x = thetas, k = df, ...)
  return(fit)
}
dtrunc <- function(z, thetas, thershold, log = FALSE) {
  probs <- log(pnorm(threshold[1], thetas) + pnorm(threshold[2], thetas, lower.tail = FALSE))
  densities <- dnorm(z, thetas, log = TRUE)
  if(log) {
    return(densities - probs)
  } else{
    return(exp(densities - probs))
  }
}
sampGMM <- function(nsamp, fit) {
  probs <- fit$lambda
  sigma <- fit$sigma
  mu <- fit$mu
  k <- sample.int(length(probs), nsamp, replace = TRUE, prob = probs)
  samp <- rnorm(nsamp, mean = mu[k], sd = sigma[k])
  return(samp)
}
sampFromDens <- function(z, fit, nsamp, threshold) {
  thetas <- sort(sampGMM(nsamp, fit))
  dens <- dtrunc(z, thetas, threshold)
  dens <- dens / sum(dens)
  cdf <- cumsum(dens)
  samp <- thetas[min(which(runif(1) < cdf))]
  return(samp)
}
computeZloglik <- function(z, support, thetadens, threshold) {
  zdens <- dtrunc(z, support, threshold)
  dens <- thetadens * zdens
  return(log(sum(dens)))
}
sampTnorm <- function(z, theta, threshold) {
  ppos <- pnorm(threshold[2], mean = theta, lower.tail = FALSE, log.p = TRUE)
  pneg <- pnorm(threshold[1], mean = theta, lower.tail = TRUE, log.p = TRUE)
  ppos <- 1 / (1 + exp(pneg - ppos))
  if(runif(1) < ppos) {
    samp <- rtruncnorm(1, mean = theta, a = threshold[2], b = Inf)
  } else {
    samp <- rtruncnorm(1, mean = theta, a = -Inf, b = threshold[1])
  }
  return(samp)
}
sampToSelections <- function(naccept, fit, threshold) {
  pselect <- pnorm(threshold[1]) + pnorm(threshold[2], lower.tail = FALSE)
  forsamp <- rep(NA, 10 / pselect)
  accepted <- 0
  sampInd <- 1
  while(accepted < naccept) {
    theta <- sampGMM(1, fit)
    z <- rnorm(1, mean = theta)
    if(z < threshold[1] | z > threshold[2]) {
      accepted <- accepted + 1
    } else {
      forsamp[sampInd] <- z
      sampInd <- sampInd + 1
    }
  }
  return(as.numeric(na.omit(forsamp)))
}
computeConditionalEst <- function(z, threshold) {
  dtr <- function(mu, z, threshold) {
    dtrunc(z, mu, threshold, log = TRUE)
  }
  cond <- optimize(dtr, interval = c(-z, z), z = z, threshold = threshold,
                   maximum = TRUE)$maximum
  return(cond)
}

# Optimizing ----------------------------
estFull <- TRUE
df <- 3
support <- seq(from = -2 * max(abs(samp)), to = 2 * max(abs(samp)), by = 0.2)
thetasamp <- samp
delay <- 150
estimates <- rep(0, length(samp))
misProp <- 0
probs <- 1
for(i in 1:300) {
  try(invisible(capture.output(thetadens <- computeGmmDens(thetasamp, df, arbvar = FALSE))))
  # print(thetadens$mu)
  if(!estFull) {
    current <- sapply(samp, sampFromDens, thetadens, 10^3, threshold)
    if(i <= Inf) {
      thetasamp <- current
    } else {
      thetasamp <- c(thetasamp, current)
    }
  } else {
    censored <- sampToSelections(length(samp), thetadens, threshold)
    current <- sapply(samp, sampFromDens, thetadens, 10^3, threshold = c(0,0))
    if(i <= Inf) {
      thetasamp <- c(current, censored)
    } else {
      thetasamp <- c(thetasamp, current, censored)
    }
  }

  if(i >= delay) {
    w <- 1 / max(i - delay, 1)
    estimates <- (1 - w) * estimates + current * w
    misProp <- (1 - w) * misProp + length(censored) / length(samp) * w
    # probs <- (1 - w) * probs + thetadens * w
  }
  bias <- mean((current - thetas), na.rm = TRUE)
  rmse <- sqrt(mean((current - thetas)^2, na.rm = TRUE))
  print(c(i, bias, rmse, length(censored) / (length(censored) + length(samp))))
}

# Evaluating model fit -----------------
loglik <- sapply(samp, computeZloglik, support, probs, threshold)
bic <- sum(-2 * loglik) + df * log(length(samp))
print(bic)
print(misProp)

# Evaluating estimation ------------
plot(support, thetadens, type = 'l')
lines(support, 0.5 * dexp(abs(support)))
abline(a = 0, b = 1)

plot(sort(samp), sort(estimates), ylim = c(-2, 7), type = "l",
     col = "blue", lwd = 2, lty = 2)
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
misProp * length(samp) / ((misProp * length(samp)) + length(samp))



