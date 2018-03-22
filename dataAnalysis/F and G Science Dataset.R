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

# plot(density(tweedy))
# lines(density(originalZ))
# lines(density(repZ))
conditional <- sapply(originalZ, truncNormMLE, threshold = threshold, sd = 1)
truncDeconvFit <- truncDeconv(originalZ, threshold = qnorm(.975),
                              meanValues = NULL,
                              twoSided = FALSE,
                              splineDegree = 3,
                              binWidth = 0.1)
empBayes <- predict(truncDeconvFit)$bayes
ptnFit <- paretoTruncNormMix(originalZ, threshold = qnorm(.975), normComps = 3, paretoComp = FALSE)
tweedy <- predict(ptnFit)

par(mfrow = c(1, 2), mar = rep(4, 4))
plot(density(repZ), main = "Original Study + Estimated Pareto Density", col = "blue", ylim = c(0, 0.5))
lines(density(originalZ), col = "black")
lines(density(tweedy), col = "red")
lines(density(conditional), col = "orange")
lines(density(empBayes), col = "green")
legend("topright", lty = 1,
       col = c("black", "blue", "red", "orange", "green"),
       legend = c("original", "replication", "F Model", "Conditional", "G Model"))


plot(originalZ, repZ, main = "Observed + Estimators")
lines(sort(originalZ), predict(ptnFit)[order(originalZ)], col = 'red', lwd = 2)
lines(sort(originalZ), conditional[order(originalZ)], col = 'blue', lwd = 2)
lines(sort(originalZ), empBayes[order(originalZ)], col = 'green', lwd = 2)
abline(a = 0, b = 1, lty = 2, lwd = 2)
legend("bottomright", lty = 1, col = c("red", "blue", "green"), legend = c("F-Model", "cond", "G-model"))

mean(originalZ - repZ)
mean(tweedy - repZ)
mean(conditional - repZ)
mean(empBayes - repZ)

mean((originalZ - repZ)^2)
mean((tweedy - repZ)^2)
mean((conditional - repZ)^2)
mean((empBayes - repZ)^2)





