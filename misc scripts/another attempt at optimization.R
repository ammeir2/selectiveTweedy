tries <- rep(0, 2)
init <- 3
lthreshold <- -3
uthreshold <- -1
sampsd <- 3

sample <- mhSampler(init, beta, lthreshold, uthreshold, sampsd,
                    1000, 100, 200, tries)
plot(density(sample))
tries
