#### plot the calculated posteriors
load("abc_estimates/abc_5k.RData")
load("abc_estimates/abc_50k.RData")
load("abc_estimates/abc_500k.RData")


abc_5k[[1]]

par(mfrow = c(8, 9), mar=c(4,4,1,1))

for (i in 1:72) {
  hist(abc_5k[[2]][[i]]$adj.values, 
       main = abc_5k[[1]][i, 2], 
       xlab = abc_5k[[1]][i, 3],
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}

abc_5k[[1]]



par(mfrow = c(8, 10), mar=c(4,4,1,1))
for (i in 1:80) {
  hist(abc_50k[[2]][[i]]$adj.values, 
       main = abc_50k[[1]][i, 2], 
       xlab = abc_50k[[1]][i, 3],
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}


par(mfrow = c(8, 9), mar=c(4,4,1,1))
for (i in 1:72) {
  hist(abc_500k[[2]][[i]]$adj.values, 
       main = abc_500k[[1]][i, 2], 
       xlab = abc_500k[[1]][i, 3],
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}