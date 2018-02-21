sampToM <- function(m, threshold, method) {
  if(method == "gmm") {
    k <- 2
    mu <- rnorm(k, 0, 0.5)
    mu <- c(-abs(mu), abs(mu))
    sd <- runif(k, min = 0.9, max = 1)
    sd <- rep(sd, 2)
    probs <- runif(k)
    probs <- rep(probs, 2) / 2
    k <- k * 2
  } else if(method == "slackslab") {
    pNull <- 0.9
  } else if(method == "exp") {
    rate <- 1
  } else if(method == "unif") {
    range <- c(-3, 3)
  }

  samp <- numeric(m)
  thetas <- numeric(m)
  censored <- rep(NA, m * 1000)
  naccept <- 0
  rejects <- 0
  while(naccept < m) {
    if(method == "exp") {
      theta <- rexp(rate)
    } else if(method == "slackslab") {
      if(runif(1) < pNull) {
        theta <- rnorm(1, 0, 0.25)
      } else {
        theta <- rnorm(1, 0, 2)
      }
    } else if(method == "unif") {
      theta <- runif(1, min = range[1], max = range[2])
    } else if(method == "gmm") {
      assign <- sample.int(k, 1, FALSE, probs)
      theta <- rnorm(1, mu[assign], sd[assign])
    }
    z <- rnorm(1, theta)
    if(z > threshold[2] | z < threshold[1]) {
      naccept <- naccept + 1
      samp[naccept] <- z
      thetas[naccept] <- theta
    } else {
      rejects <- rejects + 1
      censored[rejects] <- z
    }
  }

  return(list(samp = samp, thetas = thetas, censored = na.omit(censored)))
}

run.sim <- function(config) {
  m <- config[["m"]]
  pNull <- config[["pNull"]]
  threshold <- config[["threshold"]]
  reps <- config[["reps"]]
  method <- config[["method"]]
  threshold <- c(-abs(threshold), abs(threshold))

  results <- list()
  for(rep in 1:reps) {
    # Genearting data -------------
    threshold <- c(-2, 2)
    data <- sampToM(m, threshold, method)

    # Oracle -----------
    all <- c(data$samp, data$censored)
    samp <- data$samp
    thetas <- data$thetas
    gmmfit <- truncDeconvolution(all, threshold = c(0, 0), method = "exp-splines",
                                 df = 16,
                                 full_prior = FALSE,
                                 npost = 1,
                                 iterations = 50, delay = 20, save_densities = TRUE,
                                 progressBar = TRUE, arbvar = FALSE, maxit = 20)
    oracle <- gmmfit$estimates[1:m]
    mean((oracle - thetas) * sign(samp))
    mean((samp - thetas) * sign(samp))
    plot(density(thetas))
    lines(density(oracle), col = "red")


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

configurations <- expand.grid(m = c(100, 200),
                              pNull = c(0.95, 0.8),
                              threshold = c(1.5, 3),
                              reps = 2,
                              dataGeneration = c("exp", "slackslab", "unif", "gmm"))

results <- apply(configurations, 1, run.sim)
