n <- 2000
ind <- rbinom(n, 1, .00001)
mu <- ind * rnorm(n, sd = 0.05) + (1 - ind) * rnorm(n, sd = 1)
# mu <- rexp(n, 0.5) * (1 - 2 * ind)
# mu <- runif(n, min = -3, max = 3)
#mu <- ind * rnorm(n, -1, 1) + (1 - ind) * rnorm(n, 2, 0.5)
osamp <- rnorm(n, mu)
threshold <- 2.5
zsamp <- osamp[abs(osamp) > threshold]
mu <- mu[abs(osamp) > threshold]
J <- 4

fit <- optimExpFamily(zsamp, threshold, degree = J,
                      stepCoef = 0.15, stepRate = 0.7,
                      delay = 100, maxiter = 2500,
                      nsamp = 50)

beta <- fit$coef
J <- length(beta)
adjusted <- sign(zsamp) * pmin(abs(zsamp), abs(sapply(zsamp, function(x) x + sum(x^((1:J) - 1) * beta * 1:J))))
plot(zsamp, adjusted)
abline(a = 0, b = 1)

plot(density((adjusted - mu) * sign(zsamp)), ylim = c(0, 1))
lines(density((zsamp - mu) * sign(zsamp)), col = "red")
abline(v = 0)

plot(density(mu), col = "blue")
lines(density(adjusted), col = "black")
lines(density(zsamp), col = "red")

beta <- fit$coef
range <- seq(from = -10, to = 10, length.out = 1000)
c <- integrate(10^5, beta, 0.1)
dens <- exp(as.numeric(sapply(1:length(beta), function(x) range ^ x) %*% beta))
plot(range, dens * 1 / c, type = "l", ylim = c(0, 1))
lines(density(osamp), type = "l", col = "red")
lines(density(zsamp), type = "l", col = "blue")

mean((adjusted - mu) * sign(zsamp))
mean((zsamp - mu) * sign(zsamp))

