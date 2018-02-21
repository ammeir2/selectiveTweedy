normParetoMixEM <- function(x, threshold, barrier = 0.005, emIterations = 100, tol = 10^-5) {
  normfit <- c(mean(x), sd(x))
  paretofit <- paretoML(x, location = threshold, barrier = barrier)$par
  probs <- c(0.5, 0.5)
  loglik <- -Inf

  for(i in 1:emIterations) {
    # E-step
    normDens <- dnorm(x, normfit[1], normfit[2], log = TRUE) + log(probs[1])
    paretoDens <- paretoDens(x, threshold, paretofit[1], paretofit[2]) + log(probs[2])
    normWeights <- 1 / (1 + exp(paretoDens - normDens))
    paretoWeights <- 1 - normWeights

    # oldLoglik <- loglik
    # loglik <- mean(log(exp(normDens) + exp(paretoDens)))
    # print(loglik)
    # if(loglik - oldLoglik < tol) break

    # M-step
    wmean <- weighted.mean(x, normWeights)
    wvar <- weighted.mean((x - wmean)^2, normWeights)
    normfit <- c(wmean, sqrt(wvar))

    paretofit <- paretoML(x, location = threshold, weights = paretoWeights, barrier = barrier)$par

    probs[1] <- mean(normWeights)
    probs[2] <- 1 - probs[1]
  }

  return(list(x = x,
              normfit = normfit, paretofit = paretofit,
              probs = probs,
              posteriors = cbind(normWeights, paretoWeights),
              threshold = threshold,
              loglik = loglik, iterations = iterations))
}

normParetoMixDens <- function(x, mixfit) {
  normdens <- dnorm(x, mixfit$normfit[1], mixfit$normfit[2]) * mixfit$probs[1]
  paretodens <- paretoDens(x, mixfit$threshold, mixfit$paretofit[1], mixfit$paretofit[2], FALSE) * mixfit$probs[2]
  return(normdens + paretodens)
}

normParetoTweedy <- function(mixfit) {
  normTweed <- (mixfit$x - (mixfit$x - mixfit$normfit[1]) / mixfit$normfit[2]^2) * mixfit$posteriors[, 1]
  paretoTweed <- paretoTweedy(mixfit$x, location = mixfit$threshold, scale = mixfit$paretofit[1], shape = mixfit$paretofit[2])
  paretoTweed <- paretoTweed * mixfit$posteriors[, 2]
  return(normTweed + paretoTweed)
}
