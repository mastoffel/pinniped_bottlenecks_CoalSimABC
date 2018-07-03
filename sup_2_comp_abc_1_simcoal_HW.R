## This scripts only purpose is the repetition of ABC analyses based on individuals
# from the largest genetic cluster.

### ABC analysis part 1: Model selection and evaluation for largest clusters.
## This script will take the simulations from 1_sim_msats_simcoal as input and
# (1) plot the different summary statistics as boxplots for all models
# (2) calculate the probabilities of the models for all species
# (3) calculate the fit of the empirical data to the model
# (4) output all plots / txt files under plots and results

# files needed:
# (1) all genotypes from ../data/seal_data_largest_clust_and_pop_29.xlsx
# (2) simulation results sims_10000k.txt

# packages
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library(abc)
library(ggplot2)
library(readxl)
library(dplyr)
library(magrittr)
library(parallel)
library(readr)
library(dplyr)

######  preparation ######

# how many cores should be left free?
cores_not_to_use <- 20

###### load genetic data #######

# load all_seals data for the 29 full datasets
all_seals_full <- sealABC::read_excel_sheets("data/seal_data_largest_clust_and_pop_all_hw_30.xlsx")[1:30] # 

shortcut_save <- "_HW"
##### calculate summary statistics #####

if (!exists(paste0("data/all_sumstats_40ind_30", shortcut_save, ".txt"))) {
  cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - cores_not_to_use ))
  clusterEvalQ(cl, c(library("sealABC")))
  all_sumstats_full <- parallel::parLapply(cl, all_seals_full, 
                                           function(x) mssumstats(x, by_pop = NULL, start_geno = 4, mratio = "loose",
                                                                  rarefaction = TRUE, nresamp = 1000, nind = 40, nloc = NULL)) # 
  stopCluster(cl)
  
  # calculate the mean across clusters
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
  # write rownames as column
  all_sumstats_full <- dplyr::add_rownames(all_sumstats_full, var = "species")
  # write to txt file
  write_delim(all_sumstats_full, paste0("data/all_sumstats_40ind_30", shortcut_save, ".txt"))
} 

# get summary stats data
all_sumstats_full <- read_delim(paste0("data/all_sumstats_40ind_30", shortcut_save, ".txt"), delim = " ")

#########################################


#### select summary statistics ######
sumstats <- c("num_alleles_mean", 
              "prop_low_afs_mean",
              "mean_allele_range",  
              "mratio_mean",  
              "exp_het_mean")

all_sumstats_full <- all_sumstats_full[sumstats]


####### run abc step 1 ########

sim_name <- "sims_10000kbot500"

### load simulations, stored in main folder atm ###
path_to_sims <- paste0(sim_name, ".txt")
sims <- fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

### subsetting and definition of model selection parameters

# subset summary statistics for species of a given population size, pop_size
all_sumstats <- all_sumstats_full

# subset seal descriptive and summary data for species of a given population size, pop_size
all_seals <- all_seals_full

# parameter columns in simulation data.frame
param_start <- which(names(sims) == "sample_size")
param_end <- which(names(sims) == "range_constraint")
params <- c(param_start:param_end)
# create a character vector with models
models <- sims$model
# tolerance rate
tol <- 0.0005
# cross-validation replicates / number of replicates used to estimate the null distribution of the goodness-of-fit statistic
cv_rep <- 100
# method for model selection with approximate bayesian computation, see ?postpr
method <- 'mnlogistic'
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]

model_selection <- TRUE

if (model_selection) {
  
  dir_modselection <- "results/model_probs/"
  
  if (!dir.exists(dir_modselection)) dir.create(dir_modselection)
  
  #check probabilites for all species
  cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - cores_not_to_use))
  clusterEvalQ(cl, c(library("sealABC"), library("abc")))
  all_probs <- parallel::parApply(cl, all_sumstats, 1, abc::postpr, index = models, 
                                  sumstat = sims_stats, tol = tol, method = method)
  stopCluster(cl)
  
  save(all_probs, file = paste0("all_probs_10000k_500bot", shortcut_save, ".RData"))
  
  all_probs_df <- do.call(rbind, lapply(all_probs, function(x) {
    if (!is.numeric(x$pred)) return(NA)
    out <- round(x$pred, 3)
    out
  }))
  row.names(all_probs_df) <- names(all_seals) # new
  write.table(all_probs_df, file = paste0(dir_modselection, sim_name, shortcut_save, "_model_selection_30.txt"))
  
}


