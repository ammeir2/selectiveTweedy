# Binomial Model ------------------
set.seed(10)
sampSize <- 1000
threshold <- 2
musd <- 1
mu <- rnorm(sampSize, 0, musd)
x <- rnorm(sampSize, mu)
z <- x[x > threshold]
selectMu <- mu[x > threshold]
plot(density(x), ylim = c(0, 0.5))
lines(density(z), col = "red")
iterations <- 100
estimate <- numeric(iterations)
sdTrue <- sqrt(1 + musd^2)
sdEst <- sd(z)
naiveSD <- sqrt(mean(z^2))
censored <- numeric(sampSize - length(z))
for(i in 1:iterations) {
  for(j in 1:length(censored)) {
    samp <- Inf
    while(samp > threshold) {
      samp <- rnorm(1, 0, sdEst)
    }
    censored[j] <- samp
  }
  # sdEst <- sd(c(z, censored))
  sdEst <- sqrt(mean(c(z, censored)^2))
  estimate[i] <- sdEst
  print(c(sdEst, mean(c(z, censored))))
}
sdEst <- mean(estimate)
lines(density(c(z, censored)), col = "blue")

par(mfrow = c(1, 2), mar = c(3, 4.5, 2, 2))
range <- seq(from = -8, to = 8, by = 0.01)
plot(range, dnorm(range, mean = 0, sd = sdEst), type = "l", ylab = "Estimated Marginal Density", xlab = "",
     col = "blue")
lines(range, dnorm(range, mean = 0, sd = naiveSD), type = "l", col = "red")
lines(range, dnorm(range, mean = 0, sd = sqrt(mean(x^2))), type = "l", col = "black")

stein <- z - (sampSize - 2) / sampSize * z / sdTrue^2
naiveStein <- z - (sampSize - 2) / sampSize * z / naiveSD^2
condStein <- z - (sampSize - 2) / sampSize * z / sdEst^2
plot(density(stein - selectMu), ylim = c(0, 0.8),
     ylab = "Error Density", xlab = "", main = "")
lines(density(naiveStein - selectMu), type = "l", col = "red")
lines(density(condStein - selectMu), type = "l", col = "blue")
abline(v = 0)


# Negative Binomial Model -------------
sampleUntil <- function(m, sd, threshold){
  samp <- numeric(m * 10^3)
  rejected <- 0
  slot <- 1
  while(rejected < m) {
    samp[slot] <- rnorm(1, sd = sd)
    if(samp[slot] > threshold) {
      rejected <- rejected + 1
    } else {
      slot <- slot + 1
    }
  }
  return(samp[1:(slot - 1)])
}
set.seed(10)
m <- 100
musd <- 1
threshold <- 2
mu <- rnorm(m * 10^3, sd = musd)
z <- rnorm(m * 10^3, mean = mu)
selected <- cumsum(z > threshold)
keep <- min(which(selected == m))
z <- z[1:keep]
selected <- z > threshold
x <- z[selected]
selectMu <- mu[1:keep][selected]
iterations <- 100
sdest <- sqrt(mean(x^2))
estimates <- numeric(iterations)
for(i in 1:iterations) {
  augment <- c(x, sampleUntil(m, sdest, threshold))
  sdest <- sqrt(mean(augment^2))
  estimates[i] <- sdest
}
sdEst <- mean(estimates)
naiveSD <- sqrt(mean(x^2))
oracleSD <- sqrt(mean(z^2))

# pdf("figures/nbSteinExample.pdf",pagecentre=T, width=8,height=3.5 ,paper = "special")
par(mfrow = c(1, 2), mar = c(3, 4.5, 2, 2))
range <- seq(from = -8, to = 8, by = 0.01)
plot(range, dnorm(range, mean = 0, sd = oracleSD), type = "l",
     ylab = "Estimated Marginal Density", xlab = "",
     col = "black")
lines(range, dnorm(range, mean = 0, sd = naiveSD), type = "l", col = "red")
lines(range, dnorm(range, mean = 0, sd = sdEst, type = "l", col = "blue")

stein <- x - x / oracleSD^2
naiveStein <- x - x / naiveSD^2
condStein <- x - x / sdEst^2
plot(density(stein - selectMu), ylim = c(0, 0.8),
     ylab = "Error Density", xlab = "", main = "")
lines(density(naiveStein - selectMu), type = "l", col = "red")
lines(density(condStein - selectMu), type = "l", col = "blue")
abline(v = 0)
# dev.off()






