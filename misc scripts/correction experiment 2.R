n <- 10^4
select <- n * 0.49
ind <- rbinom(n, 1, .00001)
mu <- ind * rnorm(n, sd = 0.05) + (1 - ind) * rnorm(n, sd = 1)
#mu <- rexp(n, 0.5) * (1 - 2 * ind)
# mu <- runif(n, min = -3, max = 3)
#mu <- ind * rnorm(n, -1, 1) + (1 - ind) * rnorm(n, 2, 0.5)
osamp <- rnorm(n, mu)
J <- 4

upperSet <- order(osamp, decreasing = TRUE)[1:select]
umu <- mu[upperSet]
upper <- osamp[upperSet]
uthreshold <- min(upper) * select / (select + 1)
lowerSet <- order(osamp, decreasing = FALSE)[1:select]
lmu <- mu[lowerSet]
lower <- osamp[lowerSet]
lthreshold <- max(lower) * select / (select + 1)
zsamp <- c(upper, lower)
mu <- c(umu, lmu)

fit <- optimExpFamily(zsamp, threshold, degree = J,
                      stepCoef = 0.15, stepRate = 0.01,
                      delay = 100, maxiter = 5000,
                      nsamp = 400)

beta <- fit$coef
beta[abs(beta) < 10^-2] <- 0
J <- length(beta)
adjusted <- sign(zsamp) * pmin(abs(zsamp), abs(sapply(zsamp, function(x) x + sum(x^((1:J) - 1) * beta * (1:J)))))
plot(zsamp, adjusted)
abline(a = 0, b = 1)

par(mfrow = c(2, 1))
uadjusted <- adjusted[1:select]
ladjusted <- adjusted[(select + 1):length(adjusted)]
plot(density(ladjusted - lmu), ylim = c(0, 0.65))
lines(density(lower - lmu), col = "red")
abline(v = 0)
plot(density((uadjusted - umu)), ylim = c(0, 0.65))
lines(density(upper - umu), col = "red")
abline(v = 0)

par(mfrow = c(1, 1))
plot(density(mu), col = "blue")
lines(density(adjusted), col = "black")
lines(density(zsamp), col = "red")

#beta <- fit$coef
range <- seq(from = -10, to = 10, length.out = 1000)
c <- integrate(10^5, beta, 0.1)
dens <- exp(as.numeric(sapply(1:length(beta), function(x) range ^ x) %*% beta))
plot(range, dens * 1 / c, type = "l", ylim = c(0, 1))
lines(density(osamp), type = "l", col = "red")
lines(density(zsamp), type = "l", col = "blue")
densamp <- rejectSamp(10^4, beta, 0, 0, zsamp)
lines(density(densamp), col = "orange")

c(mean(ladjusted - lmu), mean(lower - lmu))
c(mean(uadjusted - umu), mean(upper - umu))

