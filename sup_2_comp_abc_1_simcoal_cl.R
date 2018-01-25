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
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop_29.xlsx") # 


# set to TRUE is largest clusters instead of full datasets
calc_on_cluster <- TRUE
if (calc_on_cluster){
  # for calculation on clustered data (just for retrieving the names)
  bottleneck_cl <- read_delim("../data/bottleneck_results_29_cl.txt", col_names = TRUE, delim = " ")
  all_seals_full <- all_seals_full[bottleneck_cl$id]
  shortcut_save <- "_cl"
} else {
  all_seals_full  <- all_seals_full[1:29]
  shortcut_save <- ""
}

##### calculate summary statistics #####

if (!exists(paste0("../data/all_sumstats_40ind_29", shortcut_save, ".txt"))) {
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
  write_delim(all_sumstats_full, paste0("../data/all_sumstats_40ind_29", shortcut_save, ".txt"))
} 

# get summary stats data
all_sumstats_full <- read_delim(paste0("../data/all_sumstats_40ind_29", shortcut_save, ".txt"), delim = " ")

#########################################


#### select summary statistics ######
sumstats <- c("num_alleles_mean", 
              "prop_low_afs_mean",
              "mean_allele_range",  
              "mratio_mean",  
              "exp_het_mean")

all_sumstats_full <- all_sumstats_full[sumstats]


####### run abc step 1 ########

sim_name <- "sims_10000k"

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
  
# just calculate the goodness of fit for each model?=
just_goodness_of_fit <- FALSE

if (!just_goodness_of_fit) {
  
# ##### (1) first visual checks  -------------------------------------------------
# dir_check1 <- "model_evaluation/check1_sumstats/"
# if (!dir.exists(dir_check1)) dir.create(dir_check1)
#   
# # check whether sumstats are different across models
# pdf(paste0(dir_check1, sim_name, ".pdf"))
# par(mfcol = c(3, 7), mar = c(4,4,1,1))
#   for (i in sumstats) {
#     boxplot(sims[[i]] ~ models, main = i)
#   }
# dev.off() #only 129kb in size
#   
#   
# ### (2) can abc at all distinguish between the 2 models ? ----------------------
# # leave-one-out cv and confusion matrix
# dir_check2 <- "model_evaluation/check2_models/"
# if (!dir.exists(dir_check2)) dir.create(dir_check2)
# 
# # model cv, rejection method as it throws an error otherwise
# cv.modsel <- cv4postpr(models, sims_stats, nval = cv_rep, tol = tol, method = "rejection")
# # summary
# post_probs_summary <- summary(cv.modsel)
# # model_missclassification probabilites
# write_delim(x =  round(as.data.frame(post_probs_summary$probs),3), 
#             path = paste0(dir_check2, sim_name, "_model_missclass_probs.txt"))
# # model_missclassification frequencies
# write_delim(x =  as.data.frame(post_probs_summary$conf.matrix), 
#             path = paste0(dir_check2, sim_name, "_model_missclass_freq.txt"))
# 
# # plot
# png(paste0(dir_check2, sim_name, "_confusion_mat.png"), width = 4, height = 4, units = "in", res = 300)
# plot(cv.modsel, names.arg = model_names)
# dev.off() #only 129kb in size
  
  
### (3) model selection --------------------------------------------------------
dir_modselection <- "results/model_probs/"

if (!dir.exists(dir_modselection)) dir.create(dir_modselection)
  
#check probabilites for all species
cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - cores_not_to_use))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_probs <- parallel::parApply(cl, all_sumstats, 1, abc::postpr, index = models, 
                                sumstat = sims_stats, tol = tol, method = method)
stopCluster(cl)
  
all_probs_df <- do.call(rbind, lapply(all_probs, function(x) {
    out <- round(x$pred, 3)
    out
}))
write.table(all_probs_df, file = paste0(dir_modselection, sim_name, shortcut_save, "_model_selection.txt"))

}

# save them with name
mod_sel <- read.table("../abc_analysis/results/model_probs/sims_10000k_cl_model_selection.txt", sep = " ")
all_probs_df_names <- cbind(names(all_seals_full), mod_sel)
names(all_probs_df_names)[1] <- "species"

write.table(all_probs_df_names, file = paste0(dir_modselection, sim_name, 
                                              shortcut_save, "_model_selection.txt"), row.names = FALSE)


# ### (4) Does the preferred model provide a good fit to the data?----------------
# 
# dir_modeval <- "model_evaluation/check3_modeval/"
# 
# if (!dir.exists(dir_modeval)) dir.create(dir_modeval)
# # calculate all fits
# cl <- makeCluster(getOption("cl.cores", detectCores() - cores_not_to_use))
# clusterEvalQ(cl, c(library("sealABC"), library("abc")))
# all_fits_bot <- parApply(cl, all_sumstats, 1, abc::gfit, sumstat = sims_stats, 
#                          nb.replicate = cv_rep, tol = tol, subset = models == "bot")
# 
# all_fits_neut <- parApply(cl, all_sumstats,  1, abc::gfit, sumstat = sims_stats, 
#                           nb.replicate = cv_rep, tol = tol, subset = models == "neut")
# stopCluster(cl)
#   
# all_names <- names(all_seals)
#   
# # goodness of fit plots
# pdf(file = paste0(dir_modeval, sim_name, "_goodnessoffit.pdf"), width = 18, height = 14)
# par(mfrow = c(5, 6), mar=c(4,4,1,1))
# sapply(1:length(all_names), function(x) {
#     # check whether there are NAs, and if yes, exchange with mean
#     if (any(is.na(all_fits_bot[[x]]$dist.sim))) {
#       location <- which(is.na(all_fits_bot[[x]]$dist.sim))
#       all_fits_bot[[x]]$dist.sim[location] <- mean(all_fits_bot[[x]]$dist.sim, na.rm = TRUE)
#     }
#     plot(all_fits_bot[[x]], main = paste0(all_names[x], "_bot"))
#     
#     if (any(is.na(all_fits_neut[[x]]$dist.sim))) {
#       location <- which(is.na(all_fits_neut[[x]]$dist.sim))
#       all_fits_neut[[x]]$dist.sim[location] <- mean(all_fits_neut[[x]]$dist.sim, na.rm = TRUE)
#     }
#     plot(all_fits_neut[[x]], main = paste0(all_names[x], "_neut"))
#   })
#   dev.off()
#   
#   
# # save the p_values and distances
# #if (!dir.exists("results/goodnessoffit_p")) dir.create("results/goodnessoffit_p")
#   
# all_p_bot <- unlist(lapply(all_fits_bot, function(x) summary(x)$pvalue))
# all_p_neut <-  unlist(lapply(all_fits_neut, function(x) summary(x)$pvalue))
#   
# p_df <- data.frame(species = all_names, bot_p = all_p_bot, neut_p = all_p_neut)
# write.table(p_df, file = paste0(dir_modeval, sim_name, "_p_vals_fit.txt"), row.names = FALSE)
  

  
  
  
 