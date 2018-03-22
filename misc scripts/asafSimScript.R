library('selectiveTweedy')

set.seed(500)
set.seed(7329)
par(mfrow = c(1, 1))

N <- 10000
pi0 <- .9
sd <- 1
muSD <- 3
expRate <- 2
threshold <- 2

isNull <- runif(N) <= pi0
mu <- numeric(N)
nNull <- sum(isNull)
nAlt <- N - nNull
mu[isNull] <- rexp(nNull, rate = expRate) * (1 - 2 * rbinom(nNull, 1, 0.5))
mu[!isNull] <- rnorm(nAlt, 0, muSD)
z <- rnorm(N, mu, 1)
data <- data.frame(isNull = isNull,
                   mu = mu,
                   z = z,
                   selected = z >= threshold)

plot(data$z, data$mu, col = "grey", pch = ".")
points(data$z[data$selected], data$mu[data$selected], pch = ".")
abline(a = 0, b = 1, col = "red")
abline(v = threshold, col = "black")


x <- data$z[data$selected]
splineDegree <- 10
truncDeconvFit <- truncDeconv(x, threshold = threshold,
                              meanValues = NULL,
                              twoSided = FALSE,
                              splineDegree = splineDegree,
                              binWidth = 0.1)
empBayesEst <- predict(truncDeconvFit)$bayes
points(x, empBayesEst, col='blue',cex=.25)


# eb from entire data (by deconvolution; could compare to Tweedie as well)
tau <- truncDeconvFit$stats[, "theta"]
allBayes <- deconv(tau, X = data$z, family = "Normal", pDegree = splineDegree)
allBayesEst <- deconvTweed(x, allBayes)
points(data$z[data$selected], allBayesEst, col='forestgreen', cex=.25)
legend('topleft',legend = c('estimate from truncated data', 'estimate from entire data', 'F-Model'),pch=c(1,1,1), pt.cex = .25,col=c('blue','forestgreen','burlywood3'),bty="n")

#selective MLE
abline(h=0,col='black')
abline(v=2,col='black')
# abline(v=dnorm(2)/(1-pnorm(2)),col='black')

# F-modeling
ptnFit <- paretoTruncNormMix(data$z[data$selected], threshold = 2, normComps = 3, paretoComp = FALSE)
points(data$z[data$selected], predict(ptnFit), cex = .25, col = "orange")

# compare average squared error
mean((empBayesEst-data$mu[data$selected])^2,na.rm = TRUE) # for est based on truncated data
mean((predict(ptnFit) - data$mu[data$selected])^2) # naive
mean((allBayesEst-data$mu[data$selected])^2,na.rm = TRUE) # for est based on entire data
mean((data$z[data$selected] - data$mu[data$selected])^2) # naive

# bias
mean((data$z[data$selected] - data$mu[data$selected]))
mean((empBayesEst - data$mu[data$selected]), na.rm = TRUE)
mean((predict(ptnFit) - data$mu[data$selected]))
mean((allBayesEst - data$mu[data$selected]), na.rm = TRUE)


