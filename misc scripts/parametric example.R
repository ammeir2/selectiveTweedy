library(ggplot2)
library(selectiveTweedy)

# Model 2 --------------
set.seed(1)
# Data according to censoring model
muSD <- 1
threshold <- 2
n <- 1000
mu <- rnorm(n, sd = muSD)
z <- rnorm(n, mu)
keep <- abs(z) > threshold
x <- z[keep]
cat("selected:", length(x), "\n")
paramDat <- data.frame(z = z, selected = keep, mu = mu)

# Computing bayes estimates
censorFit <- estimateCensoredParamModel(x, threshold)
censorFit$estSD
censorBayes <- censorFit$bayes

truncFit <- estimateTruncParamModel(x, threshold)
truncFit$estSD
truncBayes <- truncFit$bayesRule
paramDat$censorBayes[paramDat$selected] <- censorBayes
paramDat$truncBayes[paramDat$selected] <- truncBayes
paramDat$fullStein <- paramDat$z * (1 - 1 / var(paramDat$z))
paramDat$naiveStein <- with(paramDat, z * (1 - 1 / var(z[selected])))
paramDat$mle[paramDat$selected] <- univTruncNormMLE(paramDat$z[paramDat$selected], threshold)

ggplot(paramDat) + geom_point(aes(x = z, y = mu, col = selected), alpha = 0.5) +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(xintercept = c(-threshold, threshold), linetype = 2) +
  xlab("Observed Z-Scores") + ylab("True Mean/Estimates") +
  geom_point(aes(x = z, y = censorBayes), size = 0.2, col = "blue") +
  geom_point(aes(x = z, y = truncBayes), size = 0.2, col = "red") +
  geom_point(data = subset(paramDat, selected), aes(x = z, y = fullStein), size = 0.2, col = "black") +
  geom_point(data = subset(paramDat, selected), aes(x = z, y = naiveStein), size = 0.2, col = "orange") +
  geom_point(data = subset(paramDat, selected), aes(x = z, y = mle), size = 0.2, col = "dark green")
# ggsave("figures/parametric example scatter.pdf", height = 5, width = 8)

# Plotting histograms of errors
library(reshape2)
errors <- with(subset(paramDat, selected),
               data.frame(zsigns = sign(z),
                          mu = mu,
                          naiveStein = naiveStein - mu,
                          fullStein = fullStein - mu,
                          censorBayes = censorBayes - mu,
                          truncBayes = truncBayes - mu,
                          mle = mle - mu,
                          observed = z - mu))
apply(errors[, -1], 2, function(x, s) mean(s * x), errors$zsigns)
errors$sign <- "Positive"
errors$sign[errors$zsigns == -1] <- "Negative"
errors$zsigns <- NULL
errors <- melt(errors, id = c("mu", "sign"))
names(errors)[3:4] <- c("estimate", "error")
ggplot(errors) + geom_histogram(aes(x = error, y =  ..density.., fill = sign, col = sign), alpha = 0.1) +
  facet_wrap(~ estimate) + theme_bw() +
  geom_vline(xintercept = 0)
# ggsave("figures/parametric example error histogram.pdf", height = 5, width = 8)


