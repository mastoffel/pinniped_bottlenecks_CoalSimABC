# comparative analysis of abc output 1

# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library("abctools")
library(abc)
library(ggplot2)
library(readxl)
library(dplyr)
library(magrittr)

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28] # 

# calculate summary statistics by population 
cl <- parallel::makeCluster(getOption("cl.cores", detectCores()-20))
clusterEvalQ(cl, c(library("sealABC")))
all_sumstats_full <- parallel::parLapply(cl, all_seals_full, 
                                         function(x) mssumstats(x, by_pop = NULL, start_geno = 4, mratio = "loose",
                                                                rarefaction = TRUE, nresamp = 1000, nind = 30, nloc = 5))
stopCluster(cl)

# for clustered populations, a mean between clusters is calculated
sum_per_clust <- function(mssumstats_output) {
  if (nrow(mssumstats_output) > 1) {
    out <- as.data.frame(t(apply(mssumstats_output, 2, mean)))
  } else
    out <- mssumstats_output
  out
}
all_sumstats_full <- lapply(all_sumstats_full, sum_per_clust)
all_sumstats_full <- do.call(rbind, all_sumstats_full)
all_sumstats_full$species <- row.names(all_sumstats_full)

# names(sims)
sumstats <- c("num_alleles_mean", "num_alleles_sd",
              "exp_het_mean",  "exp_het_sd", 
              "mean_allele_size_sd","sd_allele_size_sd",
              "mean_allele_range", "sd_allele_range",
              "mratio_mean", "mratio_sd",
              "prop_low_afs_mean",
              "mean_allele_size_kurtosis", "sd_allele_size_kurtosis",
              "obs_het_mean", "obs_het_sd"
)



# plot histograms for all summary statistics for all pop_size simulations relative to the empirical data

# load simulations
path_to_sims <- paste0("sims_2000k_optimal.txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

sims_bot <- sims[sims$model == "bot", ]
sims_neut <- sims[sims$model == "neut", ]

#
plot_emp_to_sim <- function(sumstat){

pdf(file = paste0("plots_simsstats/simcoal.", sumstat, ".pdf"), width = 10, height = 20)
  
  par(mfrow = c(10, 2), mar=c(2,2,1, 1))
  for (i in 1:length(all_sumstats_full$species)){
    title <- all_sumstats_full$species[i]
    # bottleneck
    hist(sims_bot[[sumstat]], main = title, xlab = sumstat)
    abline(v=all_sumstats_full[i, sumstat], col="red", lwd = 3)
    abline(v=mean(sims_bot[[sumstat]], na.rm = TRUE), col="blue", lwd = 3)
    # neutral
    hist(sims_neut[[sumstat]], main = title, xlab = sumstat)
    abline(v=all_sumstats_full[i, sumstat], col="red", lwd = 3)
    abline(v=mean(sims_neut[[sumstat]], na.rm = TRUE), col="blue", lwd = 3)
  }
  
dev.off()
}

sapply(sumstats, plot_emp_to_sim)
  
#}




