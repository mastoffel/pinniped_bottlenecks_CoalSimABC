#### plot the calculated posteriors
load("abc_estimates/abc_full_simcoal.RData")

abc_full[[1]]

## load model probabilities
model_probs <- read.table(paste0("results/model_probs/", "sims_simcoal1000k_model_selection.txt"))
model_bot <- model_probs$bot > 0.5

abc_bot <- list(abc_full[[1]][model_bot, ], abc_full[[2]][model_bot])
abc_bot[[1]]

par(mfrow = c(4, 4), mar=c(4,4,1,1))

for (i in 97:112) {
  hist(abc_bot[[2]][[i]]$adj.values, 
       main = abc_bot[[1]][i, 2], 
       xlab = abc_bot[[1]][i, 3],
       #xlim = c(0, 0.005),
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}



length(abc_bot[[2]])


