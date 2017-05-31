#### plot the calculated posteriors for all bottlenecked species
load("abc_estimates/abc_sims_5000k.RData")

abc_full[[1]]
nrow(abc_full[[1]])

par(mfrow = c(3, 5), mar=c(4,4,1,1))

for (i in 91:105) {
  hist(abc_full[[2]][[i]]$adj.values, 
       main = abc_full[[1]][i, 2], 
       xlab = abc_full[[1]][i, 3],
       #xlim = c(0, 0.005),
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}

hist(abc_full[[2]][[85]]$adj.values, 
     main = abc_full[[1]][i, 2], 
     xlab = abc_full[[1]][i, 3],
     #xlim = c(0, 0.005),
     breaks = 1000, col = "red")
density(rgamma(10000, 1.5, rate = 1500), add = TRUE)



plot(density(abc_full[[2]][[86]]$adj.values))
lines(density(rgamma(10000, 1.5, rate = 1500)), add = TRUE)


plot(rgamma(10000, 1.5, rate = 1500), add = TRUE)




