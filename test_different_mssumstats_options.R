# seal descriptive data #
seal_descriptives <- read_excel("../data/all_data_seals.xlsx")

###### prepare data #######

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28] # 
# (1) ning NULL, nloc NULL
# (2) nind 30, nloc NULL
# (3) nind 30, nloc = 10
# (4) nind NULL, nloc = 10



##### calculate summary statistics #####
#all_sumstats_full <- lapply(all_seals_full, function(x) mssumstats(x, datatype = "microsats", by_pop = "cluster", 
#start_geno = 4, mratio = "loose"))

cl <- parallel::makeCluster(getOption("cl.cores", detectCores()-20))
clusterEvalQ(cl, c(library("sealABC")))
all_sumstats_full <- parallel::parLapply(cl, all_seals_full, 
                                         function(x) mssumstats(x, by_pop = NULL, start_geno = 4, mratio = "loose",
                                                                rarefaction = TRUE, nresamp = 1000, nind = 40, nloc = NULL)) # 
stopCluster(cl)


sum_per_clust <- function(mssumstats_output) {
  if (nrow(mssumstats_output) > 1) {
    out <- as.data.frame(t(apply(mssumstats_output, 2, mean)))
  } else
    out <- mssumstats_output
  out
}

all_sumstats_full <- lapply(all_sumstats_full, sum_per_clust)

# as data.frame
all_sumstats_full <- do.call(rbind, all_sumstats_full)
all1 <- all_sumstats_full

par(mfrow = c(1, 4), mar=c(4,4,1,1))
boxplot(all1$num_alleles_mean, ylim = c(0, 20))
boxplot(all2$num_alleles_mean, ylim = c(0, 20))
boxplot(all3$num_alleles_mean, ylim = c(0, 20))
boxplot(all4$num_alleles_mean, ylim = c(0, 20))

boxplot(all1$mean_allele_range)
boxplot(all2$mean_allele_range)
boxplot(all3$mean_allele_range)
boxplot(all4$mean_allele_range)

boxplot(all1$mratio_mean, ylim = c(0, 1))
boxplot(all2$mratio_mean, ylim = c(0, 1))
boxplot(all3$mratio_mean, ylim = c(0, 1))
boxplot(all4$mratio_mean, ylim = c(0, 1))

boxplot(all1$prop_low_afs_mean, ylim = c(0, 1))
boxplot(all2$prop_low_afs_mean, ylim = c(0, 1))
boxplot(all3$prop_low_afs_mean, ylim = c(0, 1))
boxplot(all4$prop_low_afs_mean, ylim = c(0, 1))


boxplot(all1$num_alleles_sd)
boxplot(all2$num_alleles_sd)
boxplot(all3$num_alleles_sd)
boxplot(all4$num_alleles_sd)
