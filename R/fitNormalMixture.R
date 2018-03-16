#' @import progress
fitNormalMixture <- function(x, k = 2, iter = 10, init = NULL) {
  if(is.null(init)) {
    kmFit <- kmeans(x, k)
    assign <- kmFit$cluster
    mu <- as.numeric(kmFit$centers)
    sd <- sapply(1:k, function(a) sd(x[assign == a]))
    probs <- sapply(1:k, function(a) mean(assign == a))
  } else {
    mu <- init$mu
    sd <- init$sd
    probs <- init$probs
  }

  for(i in 1:iter) {
    # E-Step ---------------------
    post <- sapply(1:k, function(a) dnorm(x, mean = mu[a], sd = sd[a]) * probs[a])
    post <- t(apply(post, 1, function(p) p / sum(p)))

    # M-Step ---------
    probs <- colMeans(post)
    mu <- apply(post, 2, function(w) weighted.mean(x, w))
    sd <- sqrt(sapply(1:k, function(a) weighted.mean((x - mu[a])^2, post[, a])))
  }

  return(list(mu = mu, sd = sd, probs = probs, post = post))
}

sampNormalMixture <- function(n, mu, sd, probs) {
  k <- length(mu)
  assign <- sample.int(k, size = n, replace = TRUE, prob = probs)
  samp <- rnorm(n, mean = mu[assign], sd = sd[assign])
  return(samp)
}

sampGammaMixture <- function(n, gammfit) {
  cluster <- sample.int(3, n, replace = TRUE, prob = gammfit$lambda)
  samp <- rgamma(n,
                shape = gammafit$gamma.pars[1, cluster],
                scale = gammafit$gamma.pars[2, cluster])
  return(samp)
}


