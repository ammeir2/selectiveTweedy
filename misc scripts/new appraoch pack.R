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
threshold <- c(-1.96, 1.96)

# Estimating -----------------
# conditional prior
gmmfit <- truncDeconvolution(samp, threshold, method = "GMM", df = 1,
                             iterations = 100, delay = 30, save_densities = FALSE,
                             progressBar = TRUE, arbvar = FALSE)
plot(sort(samp), sort(gmmfit$estimates), type = "l", ylim = c(-3, 9), lwd = 2, lty = 2)
abline(a = 0, b = 1)
expfit <- truncDeconvolution(samp, threshold, method = "exp-family", df = 4,
                             iterations = 100, delay = 30, save_densities = FALSE,
                             progress_bar = TRUE)
lines(sort(samp), sort(expfit$estimates), col = "red", lwd = 2, lty = 2)
splinefit <- truncDeconvolution(samp, threshold, method = "exp-splines", df = 8,
                                iterations = 100, delay = 30, save_densities = FALSE,
                                progress_bar = TRUE)
lines(sort(samp), sort(splinefit$estimates), col = "blue", lwd = 2, lty = 2)
legend("topleft", col = c("black", "red", "blue"),
       legend = c("GMM", "Polynomial", "Splines"), lty = 1)
points(samp, thetas)
mean(gmmfit$estimates - thetas, na.rm = TRUE)
mean(expfit$estimates - thetas, na.rm = TRUE)
mean(splinefit$estimates - thetas, na.rm = TRUE)


# Full prior
gmmfit <- truncDeconvolution(samp, threshold, method = "GMM", df = 1,
                             full_prior = TRUE,
                             iterations = 100, delay = 30, save_densities = FALSE,
                             progressBar = TRUE, arbvar = TRUE)
plot(sort(samp), sort(gmmfit$estimates), type = "l", ylim = c(-3, 9), lwd = 2, lty = 2)
abline(a = 0, b = 1)
expfit <- truncDeconvolution(samp, threshold, method = "exp-family", df = 4,
                             full_prior = TRUE,
                             iterations = 100, delay = 30, save_densities = FALSE,
                             progress_bar = TRUE)
lines(sort(samp), sort(expfit$estimates), col = "red", lwd = 2, lty = 2)
splinefit <- truncDeconvolution(samp, threshold, method = "exp-splines", df = 8,
                                full_prior = TRUE,
                                iterations = 100, delay = 30, save_densities = FALSE,
                                progress_bar = TRUE)
lines(sort(samp), sort(splinefit$estimates), col = "blue", lwd = 2, lty = 2)
legend("topleft", col = c("black", "red", "blue"),
       legend = c("GMM", "Polynomial", "Splines"), lty = 1)
points(samp, thetas)
mean(gmmfit$estimates - thetas, na.rm = TRUE)
mean(expfit$estimates - thetas, na.rm = TRUE)
mean(splinefit$estimates - thetas, na.rm = TRUE)

# Evaluating estimation ------------
estimates <- gmmfit$estimates
# estimates <- expfit$estimates
# estimates <- splinefit$estimates

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
ggplot(forplot) + geom_density(aes(x = est, col = type, linetype = type)) + theme_bw()

# plot(support, thetadens, type = "l")
expfit$estMisProp
splinefit$estMisProp
gmmfit$estMisProp


