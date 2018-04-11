library(ggplot2)
library(selectiveTweedy)
library(truncnorm)

set.seed(1)
n <- 400
muSD <- 2
threshold <- 2

# Model 1 --------------
mu <- rnorm(n, sd = muSD)
z <- numeric(n)
for(i in 1:n) {
  z[i] <- rtruncnorm(1, a = threshold, b = Inf, mean = mu[i], sd = 1)
}
modelOne <- data.frame(mu = mu, z = z, model = "Truncation",
                       threshold = threshold)

# Model 2 --------------
pselection <- pnorm(threshold, sd = sqrt(muSD^2 + 1), lower.tail = FALSE)
sample <- matrix(ncol = 2, nrow = n * 2 / pselection)
nselected <- 0
row <- 1
while(nselected < n) {
  mu <- rnorm(1, sd = muSD)
  z <- rnorm(1, mu, 1)
  sample[row, ] <- c(mu, z)
  row <- row + 1
  if(z > threshold) {
    nselected <- nselected + 1
  }
}
sample <- sample[1:(row -1), ]
modelTwo <- data.frame(mu = sample[, 1], z = sample[, 2],
                       model = "Censoring",
                       threshold = threshold)
dat <- rbind(modelOne, modelTwo)
dat$selected <- dat$z > dat$threshold
ggplot(dat) +
  geom_point(aes(x = z, y = mu, col = selected)) +
  geom_vline(xintercept = threshold, linetype = 2) +
  theme_bw() +
  facet_grid(. ~ model) +
  geom_abline(intercept = 0, slope = 1)
  # geom_vline(xintercept = 0) +
  # geom_hline(yintercept = 0)

# Censoring Estimate --------
# Model 1
x <- dat$z[dat$model == "Truncation"]
mu <- dat$mu[dat$model == "Truncation"]
censorFit <- estimateCensoredParamModel(x, threshold = c(-Inf, threshold))
censorFit$estSD
censorBayes <- censorFit$bayes
plot(x, mu)
points(x, censorBayes, col = "blue")
abline(a = 0, b = 1)
dat$censorBayes[dat$model == "Truncation"] <- censorBayes

# Model 2
x <- dat$z[dat$model == "Censoring" & dat$selected]
mu <- dat$mu[dat$model == "Censoring" & dat$selected]
censorFit <- estimateCensoredParamModel(x, threshold = c(-Inf, threshold))
censorFit$estSD
censorBayes <- censorFit$bayes
plot(x, mu)
points(x, censorBayes, col = "blue")
abline(a = 0, b = 1)
dat$censorBayes[dat$model == "Censoring" & dat$selected] <- censorBayes

# Truncation estimates ----------------------
# Model 1
x <- dat$z[dat$model == "Truncation"]
mu <- dat$mu[dat$model == "Truncation"]
truncFit <- estimateTruncParamModel(x, threshold = c(-Inf, threshold),
                                    nIntPoints = 10^3)
                                    # sigmaGrid = seq(from = 0.05, to = 2, by = 0.05))
truncFit$estSD
truncBayes <- truncFit$bayesRule
plot(x, mu)
points(x, truncBayes, col = "blue")
dat$truncBayes[dat$model == "Truncation"] <- truncBayes

# Model 2
x <- dat$z[dat$model == "Censoring" & dat$selected]
mu <- dat$mu[dat$model == "Censoring" & dat$selected]
truncFit <- estimateTruncParamModel(x, threshold = c(-Inf, threshold),
                                    nIntPoints = 10^3)
truncFit$estSD
truncBayes <- truncFit$bayesRule
plot(x, mu)
points(x, truncBayes, col = "blue")
dat$truncBayes[dat$model == "Censoring" & dat$selected] <- truncBayes

# MLE --------------------------------
# Model 1
x <- dat$z[dat$model == "Truncation"]
mu <- dat$mu[dat$model == "Truncation"]
mle <- univTruncNormMLE(x, threshold = c(-Inf, threshold), xsd = 1)
plot(x, mu, ylim = c(-5, max(x)))
points(x, mle, col = "blue")
abline(a = 0, b = 1)
dat$mle[dat$model == "Truncation"] <- mle

# Model 2
x <- dat$z[dat$model == "Censoring" & dat$selected]
mu <- dat$mu[dat$model == "Censoring" & dat$selected]
mle <- univTruncNormMLE(x, threshold = c(-Inf, threshold), xsd = 1)
plot(x, mu, ylim = c(-5, max(x)))
points(x, mle, col = "blue")
dat$mle[dat$model == "Censoring" & dat$selected] <- mle

# Oracle --------------
# Model 1
x <- dat$z[dat$model == "Truncation"]
mu <- dat$mu[dat$model == "Truncation"]
truncFit <- estimateTruncParamModel(x, threshold = c(-Inf, threshold),
                                    nIntPoints = 10^4,
                                    knownSigma = muSD)
plot(x, mu, ylim = c(-5, max(x)))
points(x, truncFit$bayesRule, col = "blue")
abline(a = 0, b = 1)
dat$oracle[dat$model == "Truncation"] <- truncFit$bayesRule

# Model 2
x <- dat$z[dat$model == "Censoring"]
mu <- dat$mu[dat$model == "Censoring"]
oracle <- x - x / (1 + muSD^2)
plot(x, mu)
points(x, oracle, col = "blue")
abline(a = 0, b = 1)
dat$oracle[dat$model == "Censoring"] <- oracle

# Plotting results -------------------------
library(reshape2)
forplot <- melt(dat, id = c("mu", "z", "model", "threshold", "selected"))
names(forplot)[6:7] <- c("Method", "Estimate")
forplot <- forplot[order(forplot$z), ]
ggplot(subset(forplot, selected)) +
  geom_point(data = subset(forplot, selected),
             aes(x = z, y = mu), alpha = 0.25, size = 1) +
  # geom_point(data = subset(forplot, !selected),
  #            aes(x = z, y = mu), alpha = 0.5, size = 0.5, color = "grey") +
  theme_bw() +
  facet_wrap(~ model) +
  geom_vline(xintercept = threshold) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_line(aes(x = z, y = Estimate, col = Method, linetype = Method),
             alpha = 1, size = 1.25) +
  ylim(-3, max(mu) + 0.5)
ggsave("figures/parametricExamplePlot.pdf", height = 4, width = 9)








