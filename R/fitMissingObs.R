sampToNumber <- function(naccept, mixfit, lower, upper, maxsamp = NULL) {
  if(is.null(maxsamp)) {
    maxsamp <- naccept * 10^3
  }
  samp <- rep(NA, maxsamp)
  accepted <- 0
  for(i in 1:length(samp)) {
    cluster <- sample.int(length(mixfit$lambda), 1, prob = mixfit$lambda)
    candidate <- rnorm(1, mixfit$mu[cluster], mixfit$sigma[cluster])
    if(candidate < lower | candidate > upper) {
      accepted <- accepted + 1
      if(accepted == naccept) {
        break
      }
    } else {
      samp[i] <- candidate
    }
  }

  return(samp[!is.na(samp)])
}

sampRejections <- function(nsamp, mixfit, lower, upper) {
  samp <- numeric(nsamp)
  slot <- 1
  while(slot <= nsamp) {
    candidate <- sampNormalMixture(1, mixfit$mu, mixfit$sigma, mixfit$lambda)
    if(candidate > lower & candidate < upper) {
      samp[slot] <- candidate
      slot <- slot + 1
    }
  }
  return(samp)
}

estimateTweedy <- function(x, k = 10,
                           M = NULL, lower = NULL, upper = NULL,
                           iterations = 400, nPost = 1,
                           nEstimates = 10, keepEach = 10,
                           mixIters = 10, verbose = FALSE,
                           lbound = 0,
                           maxsamp = NULL, ...) {
  keepIters <- seq(from = iterations, length.out = nEstimates, by = -keepEach)
  keepIters <- unique(round(keepIters))
  keepIters <- keepIters[keepIters > 0]

  cloneX <- rep(x, nPost)
  if(is.null(lower)) {
    lower <- max(x[x < 0])
  }
  if(is.null(upper)) {
    upper <- min(x[x > 0])
  }

  # print(plot(density(x)))

  selectionProp <- numeric(iterations)
  slot <- 1
  fitList <- list()

  for(i in 1:iterations) {
    if(i == 1) {
      fit <- list()
      fit$lambda <- 1
      fit$mu <- mean(x)
      fit$sigma <- sd(x)
    }

    if(is.null(M)) {
      augsamp <- sampToNumber(length(cloneX), fit, lower, upper, maxsamp = maxsamp)
    } else {
      augsamp <- sampRejections(nPost * (M - length(x)), fit, lower, upper)
    }

    augment <- c(cloneX, augsamp)
    selectionProp[i] <- length(cloneX) / length(augment)
    if(i %% 20 == 0) {
      if(verbose) print(c(i, mean(selectionProp[max(1, i - 20):i])))
    }

    if(k > 1) {
      invisible(capture.output(
        fit <- mixtools::normalmixEM(augment, k = k, maxit = mixIters, ...)
      ))
    } else {
      fit <- list()
      fit$lambda <- 1
      fit$mu <- mean(augment)
      fit$sigma <- sd(augment)
      fit$posterior <- matrix(1, nrow = length(augment))
    }

    ############# Removing low sigma
    if(lbound > 0) {
      keep <- fit$sigma >= lbound
      if(sum(keep) < 2) {
        keep <- 1:k %in% order(fit$sigma, decreasing = TRUE)[1:2]
      }
      fit$mu <- fit$mu[keep]
      fit$sigma <- fit$sigma[keep]
      fit$lambda <- fit$lambda[keep]
      fit$lambda <- fit$lambda / sum(fit$lambda)
      fit$posterior <- fit$posterior[, keep]
      fit$posterior <- t(apply(fit$posterior, 1, function(x) x / sum(x)))
    }
    #################

    if(i %in% keepIters) {
      fitList[[slot]] <- fit
      slot <- slot + 1
    }
  }

  results <- list()
  results$fit <- fitList
  results$selectionProp <- selectionProp
  results$x <- x
  class(results) <- "tweedyEstimate"
  return(results)
}

predict.tweedyEstimate <- function(obj, newX = NULL) {
  if(is.null(newX)) {
    x <- obj$x
  } else {
    x <- newX
  }

  k <- length(obj$fit[[1]]$mu)
  adjusted <- matrix(nrow = length(x), ncol = length(obj$fit))
  for(j in 1:length(obj$fit)) {
    fit <- obj$fit[[j]]
    for(i in 1:nrow(adjusted)) {
      if(!is.null(newX)) {
        post <- dnorm(x[i], fit$mu, fit$sigma) * fit$lambda
        post <- post / sum(post)
      } else {
        post <- fit$posterior[i, ]
      }
      adjustment <- -sum(post * (x[i] - fit$mu) / fit$sigma^2)
        # sapply(1:k, function(l) -post[l] * (x[i] - fit$mu[l]) / fit$sigma[l]^2))
      adjusted[i, j] <- x[i] + adjustment
    }
  }

  # return(rowMeans(adjusted))
  return(apply(adjusted, 1, median))
}

censoredJS <- function(x, M = NULL, lower = NULL, upper = NULL,
                       iterations = 200) {
  # Setting threshold
  if(is.null(lower) & is.null(upper)) {
    upper <- min(abs(x))
    lower <- -upper
  } else if(is.null(lower) & !is.null(upper)) {
    lower <- -abs(upper)
    if(lower == upper) {
      stop("invalid threshold values!")
    }
  } else if(is.null(upper) & !is.null(lower)) {
    upper <- abs(lower)
    if(lower == upper) {
      stop("invalid threshold values!")
    }
  }

  sdEst <- sqrt(mean(x^2))
  if(!is.null(M)) {
    sampM <- M
  }

  estimates <- numeric(iterations)
  Mvec <- numeric(iterations)
  for(i in 1:iterations) {
    if(is.null(M)) {
      pselect <- pnorm(upper, sd = sdEst, lower.tail = FALSE) + pnorm(lower, sd = sdEst)
      sampM <- rnbinom(1, length(x), prob = pselect)
    }

    augment <- truncnorm::rtruncnorm(sampM, a = lower, b = upper, mean = 0, sd = sdEst)
    sdEst <- sqrt(mean(c(x^2, augment^2)))
    estimates[i] <- sdEst
    Mvec[i] <- sampM
  }

  result <- list(sdEst = mean(estimates), estimates = estimates, sampM = Mvec,
                 x = x)
  class(result) <- "censoredJS"
  return(result)
}

predict.censoredJS <- function(obj, newX = NULL) {
  if(!is.null(newX)) {
    x <- newX
  } else {
    x <- obj$x
  }

  pred <- rep(0, length(x))
  ests <- obj$estimates
  Mvec <- obj$sampM
  nEstimates <- ceiling(length(ests) - length(ests) / 2)
  start <- length(ests) - nEstimates + 1
  for(i in start:length(ests)) {
    M <- Mvec[i]
    pred <- pred + nEstimates^-1 * x * (1 - (Mvec[i] - 2) / (Mvec[i] * ests[i]^2))
  }
  return(pred)
}






