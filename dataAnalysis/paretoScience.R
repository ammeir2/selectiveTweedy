convertToZ <- function(t, type, df1, df2) {
  z <- numeric(length(t))
  for(i in 1:length(t)) {
    print(i)
    if(type[i] == "F") {
      z[i] <- qnorm(pf(t[i], df1[i], df2[i]))
    } else if(type[i] == "t") {
      z[i] <- qnorm(pt(t[i], df1[i]))
    } else if(type[i] %in% c("Chi2", "Chi")) {
      z[i] <- qnorm(pchisq(t[i], df1[i]))
    } else if(type[i] == "z") {
      z[i] <- t[i]
    } else if(type[i] == "r") {
      z[i] <- atan(t[i]) * sqrt(df1[1] - 3)
    }
    print(z[i])
    a <- 0
  }

  return(z)
}

# Loading data ------------------------
RPPdata <- readRDS(file = "dataAnalysis/rppclean.rds")

sortPvals <- sort(RPPdata$T.pval.recalc.O)
sum(sortPvals[1:(length(sortPvals)-1)] == sortPvals[2:length(sortPvals)])
which(sortPvals[1:(length(sortPvals)-1)] == sortPvals[2:length(sortPvals)])

# Select the completed replication studies
RPPdata <- dplyr::filter(RPPdata, !is.na(T.pval.USE.O), !is.na(T.pval.USE.R))

# We have 99 studies for which p-values and effect sizes could be calculated
nrow(RPPdata)
# We have 97 studies for which p-values of the original effect were significantly below .05
idOK <- complete.cases(RPPdata$T.r.O,RPPdata$T.r.R)
sum(idOK)

# Starting Analysis -------------------
RPPdata$Test.statistic.O
RPPdata$T.Test.value.O
stats <- RPPdata$T.Test.value.O
statType <- RPPdata$T.Test.Statistic.O
df1 <- RPPdata$T.df1.O
df2 <- RPPdata$T.df2.O
df1[statType == "r"] <- RPPdata$T.N.O

originalZ <- qnorm(1 - pmax(RPPdata$T.pval.recalc.O, 10^-15) / 2)
repZ <- qnorm(1 - pmax(RPPdata$T.pval.recalc.R, 10^-15) / 2) * (1 - 2 * (RPPdata$Direction.R == "opposite"))
threshold <- qnorm(0.975)
whichNA <- which(is.na(originalZ) | originalZ < threshold | is.na(repZ))
originalZ <- originalZ[-whichNA]
repZ <- repZ[-whichNA]

paretofit <- paretoTweedyCorrection(originalZ, threshold)
mixfit <- normParetoMixEM(originalZ, threshold, emIterations = 500, barrier = 0)
conditional <- sapply(originalZ, truncNormMLE, threshold = threshold, sd = 1)
npTweed <- normParetoTweedy(mixfit)
bw_adj <- 1
deriv <- (bw_adj * hns(originalZ, deriv.order=1)) %>% kdde(x, h = ., deriv.order=1)
dens <- (bw_adj * hns(originalZ)) %>% kdde(x, h = .)
nonparamAdj <- pmin(predict(deriv, x = originalZ) /predict(dens, x = originalZ), 0)
sorted <- sort(originalZ)
sortAdj <- sorted + sort(nonparamAdj)
nonparamTweed <- sortAdj[match(originalZ, sorted)] %>% pmax(0)

par(mfrow = c(2, 2), mar = rep(4, 4))
plot(density(originalZ), main = "Original Study + Estimated Pareto Density")
lines(sort(originalZ), paretoDens(sort(originalZ), threshold, paretoParam[1], paretoParam[2], FALSE),
      col = "red")
lines(sort(originalZ), normParetoMixDens(sort(originalZ), mixfit), col = "green")
legend("topright", col = c("black", "red", "green"), legend = c("observed", "pareto-ML", "npmix"), lty = 1)

plot(originalZ, repZ, main = "Observed + Estimators")
lines(sort(originalZ), paretofit$estimate[order(originalZ)], col = 'red', lwd = 2)
lines(sort(originalZ), conditional[order(originalZ)], col = 'blue', lwd = 2)
lines(sort(originalZ), npTweed[order(originalZ)], col = 'green', lwd = 2)
lines(sort(originalZ), nonparamTweed[order(originalZ)], col = 'orange', lwd = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
legend("bottomright", lty = 1, col = c("red", "blue", "green", "orange"), legend = c("tweedy", "cond", "npmix", "nonparam"))

plot(density(repZ), main = "Density of Replication + Adjusted Estimates", ylim = c(0, 0.4))
lines(density(paretofit$estimate), col = "red")
lines(density(conditional), col = "blue")
lines(density(npTweed), col = "green")
legend("topright", col = c("black", "red", "blue", "green"), legend = c("rep", "tweed", "cond", "npmix"), lty = 1)
mean(originalZ - repZ)
mean(paretofit$estimate - repZ)
mean(conditional - repZ)
mean(npTweed - repZ)
mean(nonparamTweed - repZ)

median(paretofit$estimate - repZ)
median(conditional - repZ)
median(npTweed - repZ)


