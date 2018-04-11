permutationTest <- function(x, y, reps = 1000) {
  diffs <- x - y
  obsdiff <- mean(diffs)
  n <- length(diffs)
  nullSamp <- numeric(reps)
  for(i in 1:reps) {
    signs <- rbinom(n, 1, 0.5)
    nullSamp[i] <- mean(diffs * signs)
  }

  return(list(pval = mean(nullSamp > obsdiff),
         diff = obsdiff,
         samp = nullSamp))
}

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
library(ggplot2)
origSampSize <- RPPdata$T.N.O
Nratio <- sqrt(RPPdata$T.N.R / RPPdata$T.N.O)[-whichNA]
naiveRepMu <- originalZ * Nratio
naivePower <- pnorm(-threshold, mean = naiveRepMu) + 1 - pnorm(threshold, mean = naiveRepMu)
replicated <- RPPdata$Replicate.R[-whichNA]
naiveplot <- data.frame(Original = originalZ,
                        Replication = repZ,
                        Power = naivePower,
                        Replicated = replicated)
naiveplot$Replicated <- tolower(naiveplot$Replicated)
ggplot(naiveplot) +
  geom_point(aes(x = Original, y = Replication, size = Power, col = Replicated), alpha = 0.75) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
naiveperm <- permutationTest(originalZ, repZ)
naiveperm$pval

conditional <- sapply(originalZ, truncNormMLE, threshold = threshold, sd = 1)
condRepMu <- conditional * Nratio
condPower <- pnorm(-threshold, mean = condRepMu) + 1 - pnorm(threshold, mean = condRepMu)
condplot <- data.frame(Conditional = conditional,
                        Replication = repZ,
                        Power = condPower,
                        Replicated = replicated)
condplot$Replicated <- tolower(condplot$Replicated)
ggplot(condplot) +
  geom_point(aes(x = Conditional, y = Replication, size = Power, col = Replicated), alpha = 0.75) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
condperm <- permutationTest(conditional, repZ, reps = 10^4)
condperm$pval

log(conditional) - log(repZ)


