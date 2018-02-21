# Helper functions -----------------
computeGmmDens <- function(thetas, df, symmetric = TRUE, ...) {
  if(df == 1) {
    if(symmetric) {
      sd <- sqrt(mean(thetas^2))
      mu <- 0
    } else {
      mu <- mean(thetas)
      sd <- sqrt(mean((thetas - mu)^2))
    }
    fit <- list()
    fit$mu <- mu
    fit$sigma <- sd
    fit$lambda <- 1
    return(fit)
  }

  # If df > 1
  if(symmetric) {
    thetas <- c(-abs(thetas), thetas)
  }
  fit <- normalmixEM(x = thetas, k = df, ...)
  if(symmetric) {
    fit$sigma <- rep(fit$sigma, 2)
    fit$lambda <- rep(fit$lambda, 2) / 2
    fit$mu <- c(-abs(fit$mu), abs(fit$mu))
  }
  return(fit)
}
computeSplineDens <- function(thetas, support, df, symmetric = TRUE) {
  if(symmetric) {
    tab <- table(abs(thetas))
  } else {
    tab <- table((thetas))
  }
  x <- as.numeric(names(tab))
  y <- log(tab / sum(tab))
  if(is.null(df)) {
    sfit <- smooth.spline(x, y)
  } else {
    sfit <- smooth.spline(x, y, df = df)
  }
  if(symmetric) {
    dens <- predict(sfit, abs(support))$y
  } else {
    dens <- predict(sfit, (support))$y
  }
  dens <- exp(dens - max(dens))
  dens <- dens / sum(dens)
  return(dens)
}
computePolyCoefs <- function(thetas, coefs,  support, symmetric = TRUE) {
  thetas <- c(-abs(thetas), abs(thetas))
  coefs <- optim(coefs, fn = optimFuncPoly, gr = optimGradPoly, method = "BFGS",
                 support = support, thetas = thetas)$par
  return(coefs)
}
computePolyDens <- function(coefs, support) {
  df <- length(coefs)
  supp <- sapply(1:df, function(i) support^i)
  dens <- supp %*% coefs
  dens <- exp(dens - max(dens))
  dens <- dens / sum(dens)
  return(dens)
}
dtrunc <- function(z, thetas, thershold, log = FALSE) {
  probs <- log(pnorm(threshold[1], thetas) + pnorm(threshold[2], thetas, lower.tail = FALSE))
  densities <- dnorm(z, thetas, log = TRUE)
  if(log) {
    return(densities - probs)
  } else{
    return(exp(densities - probs))
  }
}
sampGMM <- function(nsamp, fit) {
  probs <- fit$lambda
  sigma <- fit$sigma
  mu <- fit$mu
  k <- sample.int(length(probs), nsamp, replace = TRUE, prob = probs)
  samp <- rnorm(nsamp, mean = mu[k], sd = sigma[k])
  return(samp)
}
gmmdens <- function(thetas, fit) {
  k <- length(fit$mu)
  dens <- matrix(nrow = length(thetas), ncol = k)
  for(i in 1:k) {
    dens[, i] <- fit$lambda[k] * dnorm(thetas, mean = fit$mu[i], sd = fit$sigma[i])
  }
  dens <- rowSums(dens)
  return(log(dens))
}
sampFromDensGMM <- function(z, fit, nsamp, threshold) {
  # fit$sigma <- fit$sigma * 2
  thetas <- sort(sampGMM(nsamp, fit))
  # impdens <- gmmdens(thetas, fit)

  # fit$sigma <- fit$sigma / 2
  # gdens <- gmmdens(thetas, fit)
  # impweights <- exp(gdens - impdens)

  dens <- dtrunc(z, thetas, threshold)# * impweights
  dens <- dens / sum(dens)
  cdf <- cumsum(dens)
  samp <- thetas[min(which(runif(1) < cdf))]
  return(samp)
}
sampFromDens <- function(z, support, thetadens, threshold) {
  zdens <- dtrunc(z, support, threshold)
  dens <- thetadens * zdens
  dens <- dens / sum(dens)
  cdf <- cumsum(dens)
  samp <- support[min(which(runif(1) < cdf))]
  return(samp)
}
sampTnorm <- function(z, theta, threshold) {
  ppos <- pnorm(threshold[2], mean = theta, lower.tail = FALSE, log.p = TRUE)
  pneg <- pnorm(threshold[1], mean = theta, lower.tail = TRUE, log.p = TRUE)
  ppos <- 1 / (1 + exp(pneg - ppos))
  if(runif(1) < ppos) {
    samp <- rtruncnorm(1, mean = theta, a = threshold[2], b = Inf)
  } else {
    samp <- rtruncnorm(1, mean = theta, a = -Inf, b = threshold[1])
  }
  return(samp)
}
sampToSelections <- function(naccept, thetadens, threshold, support, method) {
  pselect <- pnorm(threshold[1]) + pnorm(threshold[2], lower.tail = FALSE)
  forsamp <- rep(NA, 10 / pselect)
  accepted <- 0
  sampInd <- 1
  if(method != "GMM") {
    thetacdf <- cumsum(thetadens)
  }
  while(accepted < naccept) {
    if(method == "GMM") {
      theta <- sampGMM(1, thetadens)
    } else {
      theta <- support[min(which(runif(1) < thetacdf))]
    }
    z <- rnorm(1, mean = theta)
    if(z < threshold[1] | z > threshold[2]) {
      accepted <- accepted + 1
    } else {
      forsamp[sampInd] <- z
      sampInd <- sampInd + 1
    }
  }
  return(as.numeric(na.omit(forsamp)))
}
computeConditionalEst <- function(z, threshold) {
  dtr <- function(mu, z, threshold) {
    dtrunc(z, mu, threshold, log = TRUE)
  }
  cond <- optimize(dtr, interval = c(-z, z), z = z, threshold = threshold,
                   maximum = TRUE)$maximum
  return(cond)
}
sampFromDensSplines <- function(z, support, thetadens, threshold) {
  zdens <- dtrunc(z, support, threshold)
  dens <- thetadens * zdens
  dens <- dens / sum(dens)
  cdf <- cumsum(dens)
  samp <- support[min(which(runif(1) < cdf))]
  return(samp)
}
sampCondThetaPoly <- function(z, coefs, support, thetadens, threshold = 2) {
  zdens <- dtrunc(z, support, threshold)
  zdens[is.infinite(zdens)] <- 0
  dens <- zdens * thetadens
  cdf <- cumsum(dens) / sum(dens)
  return(support[min(which(runif(1) < cdf))])
}
compExpPoly <- function(coefs, support, thetas) {
  df <- length(coefs)
  supp <- sapply(1:df, function(i) support^i)
  dens <- supp %*% coefs
  dens <- exp(dens - max(dens))
  suffExp <- as.numeric(t(supp) %*% dens / sum(dens))
  return(suffExp)
}
optimFuncPoly <- function(coefs, support, thetas) {
  df <- length(coefs)
  dens <- sapply(1:df, function(i) support^i) %*% coefs
  maxdens <- max(dens)
  dens <- exp(dens - maxdens)
  Z <- sum(dens)
  thetadens <- (sapply(1:df, function(i) thetas^i) %*% coefs - maxdens) - log(Z)
  return(-mean(thetadens))
}
optimGradPoly <- function(coefs, support, thetas) {
  df <- length(coefs)
  suffStat <- colMeans(sapply(1:df, function(i) thetas^i))
  suffExp <- compExpPoly(coefs, support, thetas)
  # print(suffExp - suffStat)
  return(suffExp - suffStat)
}


truncDeconvolution <- function(samp, threshold,
                               method = c("GMM", "exp-family", "exp-splines"),
                               df = 3, support = NULL, symmetric = TRUE,
                               iterations = 100,
                               delay = 20,
                               npost = 1,
                               full_prior = FALSE,
                               save_densities = FALSE,
                               progressBar = TRUE, ...) {
  method <- method[1]
  if(!(method %in% c("GMM", "exp-family", "exp-splines"))) {
    stop("method not supported!")
  }

  if(is.null(support)) {
    if(method == "GMM") {
      support <- NULL
    } else {
      support <- seq(from = -max(abs(samp)), to = max(abs(samp)), by = 0.1)
    }
  }

  if(method == "exp-family") {
    coefs <- rep(0, df)
  } else {
    coefs <- NULL
  }

  if(full_prior) {
    misProp <- 0
  } else {
    misProp <- NULL
  }

  if(save_densities) {
    priorDensities <- list()
  } else {
    priorDensities <- NULL
  }

  if(length(threshold) != 2 |
     threshold[1] > threshold[2]) {
    stop("threshold hold must be a vector of length 2, truncation is assumed to be z<threhsold[1] or z>threshold[2]! ")
  }

  thetas <- samp
  estimates <- rep(0, length(samp))
  thetasamp <- samp
  if(progressBar) pb <- progress_bar$new(total = iterations)
  for(i in 1:iterations) {
    # Estimating Density of Prior ----------------
    if(method == "GMM") {
      try(invisible(capture.output(thetadens <- computeGmmDens(thetasamp, df, symmetric = symmetric, ...))))
    } else if(method == "exp-splines") {
      thetadens <- computeSplineDens(thetas, support, df, symmetric = symmetric)
    } else if(method == "exp-family") {
      coefs <- computePolyCoefs(thetas, coefs,  support, symmetric = symmetric)
      thetadens <- computePolyDens(coefs, support)
    }

    # Sampling from Conditional Distribution ---------------
    if(!full_prior) {
      if(method == "GMM") {
        current <- sapply(rep(samp, npost), sampFromDensGMM, thetadens, 10^3, threshold)
      } else {
        current <- sapply(rep(samp, npost), sampFromDens, support, thetadens, threshold)
      }
      thetasamp <- current
    } else {
      if(method == "GMM") {
        current <- sapply(rep(samp, npost), sampFromDensGMM, thetadens, 10^3, threshold = c(0, 0))
      } else {
        current <- sapply(rep(samp, npost), sampFromDens, support, thetadens, threshold = c(0, 0))
      }
      censored <- sampToSelections(length(samp) * npost, thetadens, threshold, support, method = method)
      thetasamp <- c(censored, current)
    }

    # Aggregating Results ------------
    if(i >= delay) {
      w <- 1 / max(i - delay, 1)
      currentEst <- colMeans(do.call("rbind", split(current, f = sapply(1:npost, function(x) rep(x, length(samp))))))
      estimates <- (1 - w) * estimates + currentEst * w
      if(full_prior) {
        misProp <- (1 - w) * misProp + length(censored) / (length(samp) * npost + length(censored)) * w
      }
    }

    # Saving density estimates --------------
    if(i > delay & save_densities) {
      priorDensities[[i - delay]] <- thetadens
    }
    if(progressBar) pb$tick()
    # print(c(i, thetadens$mu))
    # print(thetadens$sigma)
    # print(thetadens$lambda)
  }

  # Returning results --------------
  results <- list()
  results$estimates <- estimates
  results$priorDensities <- priorDensities
  results$estMisProp <- misProp
  results$expfamCoefs <- coefs
  results$support <- support
  results$call <- match.call()
  return(results)
}








