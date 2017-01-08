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

seal_descriptives <- read_excel("../data/seal_data_complete.xlsx")
seal_descriptives %<>% mutate(abund_level = ifelse(Abundance < 20000, "5k", ifelse(Abundance < 300000, "50k", "500k")))

sumstats <- c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
              "prop_low_afs_mean", "prop_low_afs_sd", "exp_het_mean", "exp_het_sd", "obs_het_mean", "obs_het_sd",
              "mean_allele_size_var", "sd_allele_size_var", "het_excess")



# plot histograms for all summary statistics for all pop_size simulations relative to the empirical data


for (pop_size in c("500k", "5k", "50k" )){ # 

# load simulations
path_to_sims <- paste0("sims_pop", pop_size, "_sim200k.txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

sims_bot <- sims[sims$model == "bot", ]
sims_neut <- sims[sims$model == "neut", ]

seal_descriptives_sub <- seal_descriptives[seal_descriptives$abund_level == pop_size, ]
# sumstat <- "num_alleles_mean"

plot_emp_to_sim <- function(sumstat){

pdf(file = paste0("plots_simsstats/", pop_size, "_", sumstat, ".pdf"), width = 10, height = 20)
  
  par(mfrow = c(10, 2), mar=c(2,2,1, 1))
  for (i in 1:length(seal_descriptives_sub$species)){
    title <- seal_descriptives_sub$species[i]
    # bottleneck
    hist(sims_bot[[sumstat]], main = title, xlab = sumstat)
    abline(v=seal_descriptives_sub[i, sumstat], col="red", lwd = 3)
    abline(v=mean(sims_bot[[sumstat]], na.rm = TRUE), col="blue", lwd = 3)
    # neutral
    hist(sims_neut[[sumstat]], main = title, xlab = sumstat)
    abline(v=seal_descriptives_sub[i, sumstat], col="red", lwd = 3)
    abline(v=mean(sims_neut[[sumstat]], na.rm = TRUE), col="blue", lwd = 3)
  }
  
dev.off()
}

sapply(sumstats, plot_emp_to_sim)
  
}




