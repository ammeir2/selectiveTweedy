n <- 10^4 * 1
threshold <- 2
mprob <- c(0.1, 0.9, 0)
muNorm <- rnorm(n, mean = 0, sd = 3)
muExp <- rexp(n, 2) #rexp(n, rate = 0.5) * (1 - 2 * rbinom(n, 1, 0.5))
muGamma <- rnorm(n, mean = 1, sd = 3)
memberships <- t(rmultinom(n, 1, mprob))
mu <- rowSums(memberships * cbind(muNorm, muExp, muGamma))
z <- mu + rnorm(n)

keep <- z > threshold
x <- z[keep]

normComps <- 1
ptnFit <- paretoTruncNormMix(x, threshold, normComps, iterations = 100)
plot(x, predict(ptnFit))
abline(a = 0, b = 1)
mean((x - mu[keep]))
mean((predict(ptnFit) - mu[keep]))




plot(sort(x), mixDens(sort(x), threshold, normFit, paretoFit, probs), col = "blue", type = "l")
lines(density(x))
abline(v = threshold)







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










