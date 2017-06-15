optimExpFamily <- function(zsamp, threshold, degree,
                           stepCoef = 0.15, stepRate = 0.7,
                           delay = 100, maxiter = 5000,
                           nsamp = 50) {
  if(round(degree / 2) * 2 != degree) {
    warning("degree of ploynomial must be even, using `degree + 1' !")
    J <- degree + 1
  } else {
    J <- degree
  }

  if(length(threshold) == 1) {
    lthreshold <- -abs(threshold)
    uthreshold <- abs(threshold)
  } else if(length(threshold) == 2) {
    lthreshold <- threshold[1]
    uthreshold <- threshold[2]
  } else {
    stop("Length of threshold must be either 1 or 2!")
  }

  beta <- rep(0, J)
  beta[J] <- - 10^-4
  beta[1] <-  0
  suffStat <- sapply(1:length(beta), function(x) mean(zsamp ^ x))

  t <- 1
  betapath <- matrix(nrow = maxiter, ncol = length(beta))
  tries <- 1
  diffMat <- matrix(nrow = maxiter, ncol = J)
  sampMat <- numeric(maxiter)
  while(t <= maxiter) {
    suffSamp <- sampleSuffStat(nsamp, beta, lthreshold, uthreshold, zsamp)
    sampMat[t] <- suffSamp[1]
    step <-  1 / max(1, t - delay)^stepRate * (suffStat - suffSamp) * stepCoef^(2 * 1:J)
    step <- sign(step) * pmin(abs(step), 0.01)
    newbeta <- beta + step
    #newbeta[J] <- min(newbeta[J], -10^-6)
    if(all(!is.nan(newbeta))) {
      beta <- newbeta
      betapath[t, ] <- beta
      diff <- (suffSamp - suffStat) / (1 + suffSamp)
      print(suffSamp - suffStat)
      diffMat[t, ] <- diff
      t <- t + 1
    }

    #print(c(t, round(colMeans(diffMat[max(1, t - 1000):(t-1), , drop = FALSE]), 2)))

    tries <- tries + 1
    if(tries > maxiter * 10) {
      print("estimation routine unstable!")
      break
    }
  }

  beta <- colMeans(betapath[(maxiter - 200):maxiter, ])

  return(list(coef = beta, optimPath = betapath))
}
