sampToNumber <- function(naccept, mixfit, lower, upper) {
  samp <- rep(NA, naccept * 10^3)
  accepted <- 0
  for(i in 1:length(samp)) {
    cluster <- sample.int(length(mixfit$lambda), 1, prob = mixfit$lambda)
    candidate <- rnorm(1, mixfit$mu[cluster], mixfit$sigma[cluster])
    if(candidate < lower | candidate > upper) {
      accepted <- accepted + 1
      if(accepted == naccept) {
        break
      }
    } else {
      samp[i] <- candidate
    }
  }

  return(samp[!is.na(samp)])
}

# Normal Mixture Model -------------------
# Generating data
par(mfrow = c(2, 2), mar = rep(3, 4))
niter <- 2 ^ 9
n <- 4000
k <- 20
nselect <- n * 0.05
musd <- runif(k, min = 0.1, max = 2)
muMean <- rnorm(k)
muprobs <- runif(k)
muprobs <- muprobs / sum(muprobs)

# k <- 3
musd <- c(0.1, 4)
muMean <- c(0, 0)
muprobs <- c(0.95, 0.05)
k <- 10

mu <- sampNormalMixture(n, muMean, musd, muprobs)
# mu <- rexp(n) * (1 - 2 * rbinom(n, 1, 0.5))
# mu <- runif(n, min = -4, max = 4)
plot(density(mu))
x <- mu + rnorm(n)
mu <- mu[order(x)]
x <- sort(x)
lines(density(x), col = "red")

# Censoring
# y <- x[c(1:nselect,  n:(n - nselect + 1))]
# true <- mu[c(1:nselect,  n:(n - nselect + 1))]
# lines(density(y), col = "blue")
# lower <- max(y[y < 0])
# upper <- min(y[y > 0])
lower <- -2
upper <- 2
y <- x[x < lower | x > upper]
true <- mu[x < lower | x > upper]


# Fitting

plot(density(x, adjust = 1), lwd = 2)
nPost <- 1
cloney <- rep(y, nPost)
# fit <- mixtools::normalmixEM(y, k = k)

curve <- 1
legend <- c()
sprop <- numeric(niter)
for(i in 1:niter) {
  if(i == 1) {
    fit <- list()
    fit$lambda <- 1
    fit$mu <- mean(y)
    fit$sigma <- sd(y)
  }

  augsamp <- sampToNumber(length(cloney), fit, lower, upper)
  augment <- c(cloney, augsamp)
  sprop[i] <- length(cloney) / length(augment)
  print(mean(sprop[max(1, i - 20):i]))

  # fit <- fitNormalMixture(augment, k = k, iter = 20)
  invisible(capture.output(mixfit <- mixtools::normalmixEM(augment, k = k, maxit = 10)))
  if(i <= Inf) {
    fit <- mixfit
  } else {
    w <- 1 / (i - (niter / 2))
    fit$mu <- mixfit$mu * w + fit$mu * (1 - w)
    fit$sigma <- mixfit$sigma * w + fit$sigma * (1 - w)
    fit$lambda <- mixfit$lambda * w + fit$lambda * (1 - w)
  }

  if(i / 2^(curve - 1) == 1 |
     i == niter) {
    if(i == niter) {
      size <- 3
    } else {
      size <- 1
    }
    legend <- c(legend, i)
    print(i)
    curve <- curve + 1
    lines(density(sampNormalMixture(n, fit$mu, fit$sigma, fit$lambda)), col = curve, lty = curve, lwd = size)
  }
}

# legend("topright", col = 0:length(legend) + 1, lty = 1:curve,
#        legend = c("target", legend))
abline(v = c(lower, upper), col = "grey")


cbind(muMean[order(muprobs)], fit$mu[order(fit$lambda)])
cbind((musd^2)[order(muprobs)], (fit$sigma^2 - 1)[order(fit$lambda)])
cbind(sort(muprobs), sort(fit$lambda))

# Adjustment ------
adjusted <- numeric(length(y))
for(i in 1:length(adjusted)) {
  adjustment <- sum(sapply(1:k, function(j) -fit$posterior[i, j] * (y[i] - fit$mu[j]) / fit$sigma[j]^2))
  adjusted[i] <- y[i] + adjustment
}

plot(density(true))
lines(density(adjusted), col = "blue")
lines(density(y), col = "red")

sqrt(mean((adjusted - true)^2)) / sqrt(mean((y - true)^2))
mean((adjusted - true) * sign(y))
mean((y - true) * sign(y))

plot(density((y - true) * sign(true)), col = "red", xlim = c(-5, 5))
lines(density((adjusted - true) * sign(true), adjust = 1), col = "blue")
abline(v = 0)

# plot(true, y)
# points(true, adjusted, col = "blue")
# abline(a = 0, b = 1)
plot(sprop, type = "l", ylim = c(0, max(sprop)))
abline(h = length(x) / length(y))
abline(h = mean(sprop[niter:(ceiling(niter/2))]), col = "red")


# Vs. Real Correction Plot
compAdjustment <- function(z, fit) {
  w <- dnorm(z, fit$mu, fit$sigma) * fit$lambda
  w <- w / sum(w)
  adjustment <- -sum((z - fit$mu) / fit$sigma^2 * w)
  return(z + adjustment)
}

trueMix <- fit
trueMix$mu <- muMean
trueMix$sigma <- sqrt(musd^2 + 1)
trueMix$lambda <- muprobs
range <- seq(from = min(x) - 1, to = max(x) + 1, by = 0.01)
range <- range[range < lower | range > upper]

par(mfrow = c(1, 1))
lrange <- range[range <= lower]
urange <- range[range >= upper]
plot(lrange, sapply(lrange, compAdjustment, trueMix), col = "blue", type = "l",
     xlim = c(min(range), max(range)), ylim = c(min(range), max(range)))
lines(urange, sapply(urange, compAdjustment, trueMix), col = "blue")
lines(lrange, sapply(lrange, compAdjustment, fit), col = "red")
lines(urange, sapply(urange, compAdjustment, fit), col = "red")
abline(v = 0, h = 0)
abline(v = c(lower, upper), col = "grey")
abline(a = 0, b = 1, lty = 2)
legend("topleft", col = c("red", "blue"), legend = c("eBayes", "Bayes"),
       lty = 1)

