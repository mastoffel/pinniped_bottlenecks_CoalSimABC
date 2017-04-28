#### plot the calculated posteriors
load("abc_estimates/abc_full.RData")

abc_full[[1]]

par(mfrow = c(4, 7), mar=c(4,4,1,1))

for (i in 197:224) {
  hist(abc_full[[2]][[i]]$adj.values, 
       main = abc_full[[1]][i, 2], 
       xlab = abc_full[[1]][i, 3],
       #xlim = c(0, 0.005),
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}



abc_full[[1]]


