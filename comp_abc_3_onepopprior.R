#### plot the calculated posteriors for all bottlenecked species
load("abc_estimates/abc_sims_10000k_neut.RData")

abc_full[[1]]
summary(abc_full[[2]][[11]])
nrow(abc_full[[1]])

par(mfrow = c(2, 5), mar=c(4,4,1,1))

for (i in 61:70) {
  hist(abc_full[[2]][[i]]$adj.values, 
       main = abc_full[[1]][i, 2], 
       xlab = abc_full[[1]][i, 3],
       #xlim = c(0, 100000),
       breaks = 1000) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}

hist(abc_full[[2]][[85]]$adj.values, 
     main = abc_full[[1]][i, 2], 
     xlab = abc_full[[1]][i, 3],
     #xlim = c(0, 0.005),
     breaks = 1000, col = "red")
density(rgamma(10000, 1.5, rate = 1500), add = TRUE)







