integrate <- function(reps, beta, theta) {
  samp <- rexp(reps, theta) * (1 - 2 * rbinom(reps, 1, 0.5))
  dirdens <- dexp(abs(samp), theta, TRUE) + log(0.5)
  fzdens <- computeFz(samp, beta, TRUE)
  w <- exp(fzdens - dirdens)
  return(mean(w))
}

n <- 10^4
ind <- rbinom(n, 1, 0.5)
mu <- rnorm(n, 0.75)
#mu <- rexp(n, 0.5) * (1 - 2 * ind)
#mu <- runif(n, min = -3, max = 3)
#mu <- ind * rnorm(n, -3) + (1 - ind) * rnorm(n, 2)
osamp <- rnorm(n, mu)
threshold <- 0.5
zsamp <- osamp[abs(osamp) > threshold]
J <- 4

beta <- rep(0, J)
beta[J] <- - 10^-4
beta[1] <-  0
suffStat <- sapply(1:length(beta), function(x) mean(zsamp ^ x))
t <- 1
maxiter <- 4000
betapath <- matrix(nrow = maxiter, ncol = length(beta))
tries <- 1
diffMat <- matrix(nrow = maxiter, ncol = J)
sampMat <- numeric(maxiter)
while(t <= maxiter) {
  suffSamp <- sampleSuffStat(1, beta, threshold, zsamp)
  sampMat[t] <- suffSamp[1]
  step <-  1 / max(1, t - 200)^0.7 * (suffStat - suffSamp) * 0.15^(2 * 1:J)
  step <- sign(step) * pmin(abs(step), 0.05)
  newbeta <- beta + step
  newbeta[J] <- min(newbeta[J], -10^-J)
  if(all(!is.nan(newbeta))) {
    beta <- newbeta
    betapath[t, ] <- beta
    diff <- (suffSamp - suffStat) / (1 + suffSamp)
    diffMat[t, ] <- diff
    t <- t + 1
  }

  print(c(t, round(colMeans(diffMat[max(1, t - 1000):(t-1), , drop = FALSE]), 2)))

  tries <- tries + 1
  if(tries > maxiter * 10) {
    print("estimation routine unstable!")
    break
  }
}

forplot <- melt(betapath)
names(forplot) <- c("iter", "degree", "est")
ggplot(forplot) + geom_line(aes(x = iter, y = est, col = factor(degree)))

beta <- colMeans(betapath[(maxiter - 100):maxiter, ])
range <- seq(from = -10, to = 10, length.out = 1000)
c <- integrate(10^5, beta, 0.1)
dens <- exp(as.numeric(sapply(1:length(beta), function(x) range ^ x) %*% beta))
plot(range, dens * 1 / c, type = "l", ylim = c(0, 0.3))
lines(density(osamp), type = "l", col = "red")
lines(density(zsamp), type = "l", col = "blue")
lines(density(sampMat[(maxiter - 1000):maxiter]), col = "green")

