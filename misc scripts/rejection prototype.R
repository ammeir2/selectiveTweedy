fz <- function(z, beta) {
  result <- numeric(length(z))
  for(i in 1:length(z)) {
    zi <- z[i]^(1:length(beta))
    result[i] <- exp(sum(beta * zi))
  }
  return(result)
}

rdir <- function(n, theta) {
  sign <- 1 - 2*rbinom(n, 1, 0.5)
  samp <- rexp(n, theta)
  return(samp * sign)
}

densRatio <- function(z, beta, theta) {
  fz(z, beta) / (0.5 * dexp(z, theta))
}

dlaplace <- function(z, theta, log = FALSE) {
  val <- 0.5 * theta * exp(- abs(z) * theta)
  if(log) val <- log(val)
  return(val)
}



threshold <- 2
beta <- c(-0.5, 0, 1, -0.5)
theta <- 0.5
epsilon <- 10^-5

# Finding maximum of ratio -----------
optimTheta <- getBestTheta(beta, nvals = 4)
M <- optimTheta[2]
theta <- optimTheta[1]
n <- 10^4
range <- seq(from = -10, to = 10, length.out = 10^3)
plot(range, sapply(range, fz, beta), type = "l" , col = "red")
lines(range, dlaplace(range, theta))
ratios <- log(fz(range, beta)) - dlaplace(range, theta, log = TRUE)
plot(range, ratios)
M <- exp(max(ratios))

nToSamp <- min(max(n * M * 2, n), 10^5)
samp <- rdir(nToSamp, theta)
ratios <- fz(samp, beta) / (0.5 * dexp(abs(samp), 1 / theta) * M)
keep <- runif(nToSamp) < ratios
samp <- samp[keep]
plot(density(samp))
lines(range, fz(range, beta) * 0.25, type = "l" , col = "red")


