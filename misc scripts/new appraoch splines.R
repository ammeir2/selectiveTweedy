# Loading Science Data --------------------

# Artificial Data ------------------------
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

# Loading science data ------------------
# scidat <- readRDS("dataAnalysis/sciDat.rds")
# samp <- scidat$orig
# thetas <- scidat$rep

# Defining helper functions ---------------
computeSplineDens <- function(thetas, support, df, symmetric = FALSE) {
  if(symmetric) {
    tab <- table(abs(thetas))
  } else {
    tab <- table((thetas))
  }
  x <- as.numeric(names(tab))
  y <- log(tab / sum(tab))
  if(is.null(df)) {
    sfit <- smooth.spline(x, y)
  } else {
    sfit <- smooth.spline(x, y, df = df)
  }
  if(symmetric) {
    dens <- predict(sfit, abs(support))$y
  } else {
    dens <- predict(sfit, (support))$y
  }
  dens <- exp(dens - max(dens))
  dens <- dens / sum(dens)
  return(dens)
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
sampFromDens <- function(z, support, thetadens, threshold) {
  zdens <- dtrunc(z, support, threshold)
  dens <- thetadens * zdens
  dens <- dens / sum(dens)
  cdf <- cumsum(dens)
  samp <- support[min(which(runif(1) < cdf))]
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
sampToSelections <- function(naccept, support, thetadens, threshold) {
  pselect <- pnorm(threshold[1]) + pnorm(threshold[2], lower.tail = FALSE)
  forsamp <- rep(NA, 10 / pselect)
  accepted <- 0
  thetacdf <- cumsum(thetadens)
  sampInd <- 1
  while(accepted < naccept) {
    theta <- support[min(which(runif(1) < thetacdf))]
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
estFull <- FALSE
df <- 4
support <- seq(from = -2 * max(abs(samp)), to = 2 * max(abs(samp)), by = 0.2)
thetasamp <- samp
delay <- 1000
estimates <- rep(0, length(samp))
misProp <- 0
probs <- 1
for(i in 1:2000) {
  thetadens <- computeSplineDens(thetasamp, support, df, symmetric = TRUE)
  if(!estFull) {
    current <- sapply(samp, sampFromDens, support, thetadens, threshold)
    if(i <= Inf) {
      thetasamp <- current
    } else {
      thetasamp <- c(thetasamp, current)
    }
  } else {
    censored <- sampToSelections(length(samp), support, thetadens, threshold)
    current <- sapply(samp, sampFromDens, support, thetadens, threshold = c(0, 0))
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
    probs <- (1 - w) * probs + thetadens * w
  }
  bias <- mean((current - thetas), na.rm = TRUE)
  rmse <- sqrt(mean((current - thetas)^2, na.rm = TRUE))
  print(c(i, bias, rmse, length(censored)))
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

plot(samp, estimates)
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




