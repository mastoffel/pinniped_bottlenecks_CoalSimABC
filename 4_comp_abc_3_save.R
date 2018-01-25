# save abc estimates to RData file
# replace filename neut with bot when working on bottleneck model estimates.

load("abc_estimates/abc_sims_10000k_neut.RData")
abc_full[[1]]
length(abc_full[[2]])

# aim: save everything in one data.frame
abc_complete <- abc_full[[1]]

# save adjusted and unadjusted values to list
all_vals_adj <- list()
all_vals_unadj <- list()
for (i in 1:length(abc_full[[2]])) {
  all_vals_adj <- c(all_vals_adj, list(abc_full[[2]][[i]]$adj.values))
  all_vals_unadj <- c(all_vals_unadj, list(abc_full[[2]][[i]]$unadj.values))
}
  
# add values to abc_complete
abc_complete$adj_vals <- all_vals_adj
abc_complete$unadj_vals <- all_vals_unadj

save(abc_complete, file = paste0("abc_estimates/abc_10000k_neut_complete.RData"))
