results <- list()
for(i in 1:100) {
  # try(results[[i]] <- readRDS(file = paste("simulations/results/slackSlab 1 seed ", i, ".rds",sep = "")))
  try(results[[i]] <- readRDS(file = paste("simulations/results/slackSlab 2 laplace seed ", i, ".rds",sep = "")))
}

results <- do.call("c", do.call("c", results))
forplot <- list()
for(i in 1:length(results)) {
  res <- results[[i]]
  dat <- data.frame(t(res$param),
                    naive = mean((res$naive - res$true)),
                    oracle = mean((res$oracle - res$true)),
                    nTweedy = mean((res$nTweedy - res$true)),
                    selective = mean((res$sTweedy - res$true)))
  forplot[[i]] <- dat
}

forplot <- data.table::rbindlist(forplot)

library(dplyr)
library(reshape2)
library(ggplot2)
forplot$reps <- NULL
forplot <- melt(forplot, id = c("m", "pNull", "threshold", "muSD"))
names(forplot)[5:6] <- c("method", "error")


ggplot(forplot, aes(x = method, y = error)) +
  geom_boxplot() + theme_bw() +
  facet_grid(muSD ~ m, labeller = "label_both") +
  geom_hline(yintercept = 0) + ylim(-2.5, 2.5)
