run.sim <- function(config) {
  m <- config[["m"]]
  pNull <- config[["pNull"]]
  threshold <- config[["threshold"]]
  muSD <- config[["muSD"]]
  reps <- config[["reps"]]

  muMu <- c(0, 0)
  muSD <- c(0.1, muSD)
  probs <- c(pNull, 1 - pNull)

  k <- 5

  results <- list()
  for(rep in 1:reps) {
    # Genearting data -------------
    mu <- sampNormalMixture(m, muMu, muSD, probs)
    z <- rnorm(m, mu, 1)
    x <- z[(z) > threshold]
    whichSelected <- which((z) > threshold)
    true <- mu[(z) > threshold]

    # Naive naive -----------
    naive <- x
    plot(density(naive - true))
    abline(v = 0)

    # Oracle Tweedy -------
    invisible(capture.output(fullfit <- mixtools::normalmixEM(z, k = k)))
    oracleTweedy <- numeric(length(x))
    slot <- 1
    for(i in whichSelected) {
      adjustment <- sum(sapply(1:k, function(j) -fullfit$posterior[i, j] * (z[i] - fullfit$mu[j]) / fullfit$sigma[j]^2))
      oracleTweedy[slot] <- z[i] + adjustment
      slot <- slot + 1
    }
    lines(density(oracleTweedy - true), col = "red")


    # Naive Tweedy ---------------
    invisible(capture.output(naivefit <- mixtools::normalmixEM(x, k = 2)))
    naiveTweedy <- numeric(length(x))
    for(i in 1:length(oracleTweedy)) {
      adjustment <- sum(sapply(1:2, function(j) -naivefit$posterior[i, j] * (x[i] - naivefit$mu[j]) / naivefit$sigma[j]^2))
      naiveTweedy[i] <- x[i] + adjustment
    }
    lines(density(naiveTweedy - true), col = "blue")


    # Selective Tweedy ------------
    # plot(density(z))
    selectiveFit <- estimateTweedy(x, k = k, M = m,
                                   lower = -Inf, upper = threshold,
                                   iterations = 300, nPost = 1,
                                   nEstimates = 10, keepEach = 10,
                                   mixIters = 10)
    selective <- predict(selectiveFit)
    lines(density(selective - true), col = "green")


    print(c(mean((naive - true)^2),
          mean((oracleTweedy - true)^2),
          mean((naiveTweedy - true)^2),
          mean((selective - true)^2)))

    iterResult <- list()
    iterResult$param <- config
    iterResult$true <- true
    iterResult$naive <- naive
    iterResult$oracle <- oracleTweedy
    iterResult$nTweedy <- naiveTweedy
    iterResult$sTweedy <- selective
    results[[rep]] <- iterResult
  }

  return(results)
}

configurations <- expand.grid(m = c(500, 5000),
                              pNull = c(0.95, 0.8),
                              threshold = 2,
                              muSD = 2,
                              reps = 2)

results <- apply(configurations, 1, run.sim)
