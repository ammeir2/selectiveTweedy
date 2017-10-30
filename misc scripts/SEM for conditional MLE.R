# Normal Mixture Model -------------------
# Generating data
par(mfrow = c(2, 2), mar = rep(3, 4))
niter <- 256
n <- 2000
k <- 5
nselect <- n * 0.05
musd <- runif(k, min = 0.1, max = 2)
muMean <- rnorm(k)
muprobs <- runif(k)
muprobs <- muprobs / sum(muprobs)

# k <- 3
musd <- c(0.1, 4)
muMean <- c(0, 0)
muprobs <- c(0.99, 0.01)
k <- 10

mu <- sampNormalMixture(n, muMean, musd, muprobs)
# mu <- rexp(n) * (1 - 2 * rbinom(n, 1, 0.5))
plot(density(mu))
x <- mu + rnorm(n)
mu <- mu[order(x)]
x <- sort(x)
signx <- sign(x)
lines(density(x), col = "red")

# Censoring
y <- x[c(1:nselect,  n:(n - nselect + 1))]
true <- mu[c(1:nselect,  n:(n - nselect + 1))]
lines(density(y), col = "blue")
lower <- max(y[y < 0])
upper <- min(y[y > 0])

# Fitting
plot(density(x, adjust = 1), lwd = 2)
nPost <- 1
cloney <- rep(y, nPost)
# fit <- mixtools::normalmixEM(y, k = k)

curve <- 1
legend <- c()
for(i in 1:niter) {
  if(i == 1) {
    augsamp <- rnorm(n * nPost * 20, mean = mean(y), sd = sd(y))
  } else {
    augsamp <- sampNormalMixture(n * nPost * 20, fit$mu, fit$sigma, fit$lambda)
  }
  fullaug <- augsamp
  augsamp <- augsamp[augsamp > lower & augsamp < upper]
  augsamp <- augsamp[1:round(nPost * (n - nselect * 2))]
  augment <- c(cloney, augsamp)

  # fit <- fitNormalMixture(augment, k = k, iter = 20)
  fit <- mixtools::normalmixEM(augment, k = k, maxit = 10)

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
    lines(density(fullaug), col = curve, lty = curve, lwd = size)
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


