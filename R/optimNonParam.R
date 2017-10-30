optimExpFamily <- function(zsamp, threshold, degree,
                           stepCoef = 0.15, stepRate = 0.7,
                           delay = 100, maxiter = 5000,
                           nsamp = 50) {
  # if(round(degree / 2) * 2 != degree) {
  #   warning("degree of ploynomial must be even, using `degree + 1' !")
  #   J <- degree + 1
  # } else {
  #   J <- degree
  # }

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
  beta[1] <- 0
  beta[2] <- -1
  suffStat <- sapply(1:length(beta), function(j) mean(zsamp^j) / max(abs(zsamp))^j)

  t <- 1
  betapath <- matrix(nrow = maxiter, ncol = length(beta))
  tries <- 1
  diffMat <- matrix(nrow = maxiter, ncol = J)
  init <- rexp(1)
  sampsd <- 1
  while(t <= maxiter) {
    mhtries <- rep(0, 2)
    init <- sample(zsamp, 1)
    samp <-  mhSampler(init, beta, lthreshold, uthreshold, sampsd,
                       100, 10, 500, mhtries)
    suffSamp <- sapply(1:J, function(j) mean(samp^j) / max(abs(zsamp))^j)

    step <-  1 / max(1, t - delay)^stepRate * (suffStat - suffSamp) * stepCoef
    step <- sign(step) * pmin(abs(step), 0.01)
    newbeta <- beta + step
    # newbeta[J] <- min(newbeta[J], -10^-3)
    beta <- newbeta
    betapath[t, ] <- beta
    diff <- (suffSamp - suffStat)
    print(round(suffSamp - suffStat, 3))
    diffMat[t, ] <- diff
    t <- t + 1

    # if(all(!is.nan(newbeta))) {
    # }

    # Updating MH coef
    mhrate <- mhtries[2] / mhtries[1]
    if(mhrate > 0.234) {
      sampsd <- sampsd * 1.04
    } else {
      sampsd <- sampsd * 0.95
    }
    sampsd <- max(sd(zsamp) * 1.5, min(sampsd, sd(zsamp) / 10))

    #print(c(t, round(colMeans(diffMat[max(1, t - 1000):(t-1), , drop = FALSE]), 2)))

    tries <- tries + 1
    if(tries > maxiter * 10) {
      print("estimation routine unstable!")
      break
    }
  }

  colMeans(diffMat[(maxiter - 500):maxiter, ])
  beta <- colMeans(betapath[(maxiter - 200):maxiter, ])

  return(list(coef = beta, optimPath = betapath))
}
