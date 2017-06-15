rejectSamp <- function(nsamp, beta, lthreshold, uthreshold, zsamp) {
  # polynom <- c(beta[1] + theta, beta[2:length(beta)]) * 1:length(beta)
  # root <- polyroot(polynom)
  # root <- Re(root[abs(Im(root)) < epsilon])
  # root <- root[root > 0]
  # values <- computeFz(root, beta, TRUE) - (0.5 * dexp(abs(root), theta, TRUE))
  # posvalues <- exp(values)
  #
  # polynom <- c(beta[1] - theta, beta[2:length(beta)]) * 1:length(beta)
  # root <- polyroot(polynom)
  # root <- Re(root[abs(Im(root)) < epsilon])
  # root <- root[root <= 0]
  # values <- computeFz(root, beta, TRUE) - (0.5 * dexp(abs(root), theta, TRUE))
  # negvalues <- exp(values)
  #
  # values <- c(posvalues, negvalues)
  # M <- max(values)

  theta <- 1 / (quantile(abs(zsamp), c(0.05, 0.25, 0.5, 0.75, 0.95)))
  range <- seq(from = -10, to = 10, by = 0.01)
  M <- numeric(length(theta))
  for(i in 1:length(M)) {
    M[i] <- exp(max(computeFz(range, beta, log = TRUE) - log(0.5) - log(theta[i]) + abs(range) * theta[i]))
  }

  return(rejectSampCpp(nsamp, nsamp * 1000, lthreshold, uthreshold, M, theta, beta))
}

sampleSuffStat <- function(nsamp, beta, lthreshold, uthreshold, zsamp) {
  samp <- rejectSamp(nsamp, beta, lthreshold, uthreshold, zsamp)
  samp <- samp[(samp != 0) & !is.na(samp)]

  samp <- sapply(1:length(beta), function(x) mean(samp ^ x))
  return(samp)
}
