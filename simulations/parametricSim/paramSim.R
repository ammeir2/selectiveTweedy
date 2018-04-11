library(ggplot2)
library(selectiveTweedy)
library(truncnorm)
library(magrittr)

set.seed(1)
n <- 400
muSD <- 2
threshold <- 2

runSim <- function(config) {
  n <- config[["n"]]
  reps <- config[["reps"]]
  muSD <- config[["muSD"]]
  threshold <- config[["threshold"]]
  selection <- config[["selection"]]

  if(selection == 1) {
    threshold <- c(-Inf, threshold)
  } else {
    threshold <- c(-abs(threshold), abs(threshold))
  }

  resultList <- vector(reps, mode = "list")
  estimates <- c("censor", "trunc", "mle", "oracle", "stein", "naive")
  meanErrors <- matrix(nrow = reps, ncol = length(estimates)) %>% data.frame()
  names(meanErrors) <- estimates
  medianErrors <- meanErrors
  rmses <- medianErrors
  sdEst <- matrix(nrow = reps, ncol = 2) %>% data.frame()
  names(sdEst) <- c("censor", "trunc")
  for(m in 1:reps) {
    # Generating Data
    sample <- matrix(ncol = 2, nrow = n)
    nselected <- 0
    while(nselected < n) {
      mu <- rnorm(1, sd = muSD)
      z <- rnorm(1, mu, 1)
      if(z > threshold[2] | z < threshold[1]) {
        nselected <- nselected + 1
        sample[nselected, ] <- c(mu, z)
      }
    }
    x <- sample[, 2]
    mu <- sample[, 1]

    # Censoring Model --------
    censorFit <- estimateCensoredParamModel(x, threshold = threshold)
    censorBayes <- censorFit$bayes
    censorSD <- censorFit$estSD

    # Truncation Model -------
    truncFit <- estimateTruncParamModel(x, threshold = threshold,
                                        nIntPoints = 1000)
    truncSD <- truncFit$estSD
    truncBayes <- truncFit$bayesRule

    # Other -------
    mle <- univTruncNormMLE(x, threshold = threshold, xsd = 1)
    oracle <- x - x / (1 + muSD^2)
    stein <- x - x / var(x)
    naive <- x

    # Reporting -------
    meanErrors[m, ] <- c(censor = mean(censorBayes - mu),
                             trunc = mean(truncBayes - mu),
                             mle = mean(mle - mu),
                             oracle = mean(oracle - mu),
                             stein = mean(stein - mu),
                             naive = mean(naive - mu))
    medianErrors[m, ] <- c(censor = median(censorBayes - mu),
                               trunc = median(truncBayes - mu),
                               mle = median(mle - mu),
                               oracle = median(oracle - mu),
                               stein = median(stein - mu),
                               naive = median(naive - mu))
    rmses[m, ] <- c(censor = sqrt(mean((censorBayes - mu)^2)),
                       trunc = sqrt(mean((truncBayes - mu)^2)),
                       mle = sqrt(mean((mle - mu)^2)),
                       oracle = sqrt(mean((oracle - mu)^2)),
                       stein = sqrt(mean((stein - mu)^2)),
                       naive = sqrt(mean((naive - mu)^2)))
    sdEst[m, ] <- c(censor = censorSD, trunc = truncSD)
    print(c(m, colMeans(medianErrors[1:m, , drop = FALSE])))
    print(c(colMeans(sdEst[1:m, , drop = FALSE])))
  }

  results <- list()
  results$meanErrors <- meanErrors
  results$medianErrors <- medianErrors
  results$rmses <- rmses
  results$sdEst <- sdEst
  results$config <- config
  return(results)
}

configurations <- expand.grid(n = 400,
                              reps = 1000,
                              muSD = 2,
                              threshold = 2,
                              selection = 1)
set.seed(1)
# results <- apply(configurations, 1, runSim)
# saveRDS(results, "simulations/parametricSim/results/paramSim_A_seed1_1000reps.rds")

# Processing results ---------
results <- readRDS( "simulations/parametricSim/results/paramSim_A_seed1_1000reps.rds")
results <- results[[1]]

# mean error histograms
meanError <- results$meanErrors
meanError <- melt(meanError)
names(meanError) <- c("Estimate", "Error")
ggplot(meanError) +
  geom_histogram(aes(x = Error, y = ..density..), bins = 50, col = "black", fill = "white") +
  facet_wrap(~ Estimate, scales = "free") +
  theme_bw() +
  geom_vline(xintercept = 0, col = "red") +
  ggtitle("Mean Error")
# ggsave("figures/paramExampleMeanError.pdf", width = 8, height = 5)

# median error histogram
medianErrors <- results$medianErrors
medianErrors <- melt(medianErrors)
names(medianErrors) <- c("Estimate", "Error")
ggplot(medianErrors) +
  geom_histogram(aes(x = Error, y = ..density..), bins = 50, col = "black", fill = "white") +
  facet_wrap(~ Estimate, scales = "free") +
  theme_bw() +
  geom_vline(xintercept = 0, col = "red") +
  ggtitle("Median Error")
# ggsave("figures/paramExampleMedianError.pdf", width = 8, height = 5)

# rmse error histogram
rmses <- results$rmses
rmses <- melt(rmses)
names(rmses) <- c("Estimate", "RMSE")
ggplot(rmses) +
  geom_histogram(aes(x = RMSE), bins = 50, col = "black", fill = "white") +
  facet_wrap(~ Estimate, scales = "free") +
  theme_bw() +
  geom_vline(xintercept = 0, col = "red") +
  ggtitle("RMSE")
# ggsave("figures/paramExampleRMSE.pdf", width = 8, height = 5)

# standard deviation estimates
sdest <- results$sdEst
sdest <- melt(sdest)
names(sdest) <- c("Model", "SD_estimate")
ggplot(sdest) +
  geom_histogram(aes(x = SD_estimate), bins = 50, col = "black", fill = "white") +
  theme_bw() +
  facet_wrap(~ Model) +
  geom_vline(xintercept = 2, col = "red") +
  ggtitle("Estimated Standard Deviation of Mean Distribution")
# ggsave("figures/paramExampleSDests.pdf", width = 8, height = 3)

