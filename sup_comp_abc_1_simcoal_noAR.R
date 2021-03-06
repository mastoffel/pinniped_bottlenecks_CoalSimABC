### COMPARATIVE ANALYSIS WITH ABC BASED ON SIMULATIONS AND EMPIRICAL DATA

## This script will use all summary statistics except for AR

## This script will
# (1) plot the different summary statistics as boxplots for all models
# (2) calculate the probabilities of the models for all species
# (3) calculate the fit of the empirical data to the model
# (4) output all plots / txt files under plots and results


# analyse simulated microsatellites under three different models000
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

######  load data ######

# seal descriptive data #
# seal_descriptives <- read_excel("../data/all_data_seals.xlsx")

# how many cores should be left free?
cores_not_to_use <- 10

###### prepare data #######

# load all_seals data for the 29 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop_29.xlsx")[1:29] # 


##### calculate summary statistics #####

if (!exists("../data/all_sumstats_40ind_29.txt")) {
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
  write_delim(all_sumstats_full, "../data/all_sumstats_40ind_29.txt")
} 

# get summary stats data
all_sumstats_full <- read_delim("../data/all_sumstats_40ind_29.txt", delim = " ")

#########################################


#### select summary statistics ######

sumstats <- c("prop_low_afs_mean",
              "mean_allele_range",  
              "mratio_mean",  
              "exp_het_mean")

# "mean_allele_size_sd", "mean_allele_range", 
#"mean_allele_size_kurtosis", "sd_allele_size_kurtosis"
# ) # , 
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
cv_rep <- 2
# method for model selection with approximate bayesian computation, see ?postpr
method <- 'mnlogistic'
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]


just_goodness_of_fit <- FALSE

if (!just_goodness_of_fit) {
  
  ##### (1) first visual checks  -------------------------------------------------
  dir_check1 <- "model_evaluation/check1_sumstats/"
  if (!dir.exists(dir_check1)) dir.create(dir_check1)
  
  # check whether sumstats are different across models
  pdf(paste0(dir_check1, sim_name, "noAR.pdf"))
  par(mfcol = c(3, 7), mar = c(4,4,1,1))
  for (i in sumstats) {
    boxplot(sims[[i]] ~ models, main = i)
  }
  dev.off() #only 129kb in size
  
  
  ### (2) can abc at all distinguish between the 2 models ? ----------------------
  # leave-one-out cv and confusion matrix
  dir_check2 <- "model_evaluation/check2_models/"
  if (!dir.exists(dir_check2)) dir.create(dir_check2)
  
  # model cv, rejection method as it throws an error otherwise
  cv.modsel <- cv4postpr(models, sims_stats, nval = cv_rep, tol = tol, method = "rejection")
  # summary
  post_probs_summary <- summary(cv.modsel)
  # model_missclassification probabilites
  write_delim(x =  round(as.data.frame(post_probs_summary$probs),3), 
              path = paste0(dir_check2, sim_name, "_model_missclass_probs_noAR.txt"))
  # model_missclassification frequencies
  write_delim(x =  as.data.frame(post_probs_summary$conf.matrix), 
              path = paste0(dir_check2, sim_name, "_model_missclass_freq_noAR.txt"))
  
  # plot
  png(paste0(dir_check2, sim_name, "_confusion_mat_noAR.png"), width = 4, height = 4, units = "in", res = 300)
  plot(cv.modsel, names.arg = model_names)
  dev.off() #only 129kb in size
  
  
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
  write.table(all_probs_df, file = paste0(dir_modselection, sim_name, "_model_selection_noAR.txt"))
  
}

### (4) Does the preferred model provide a good fit to the data?----------------

dir_modeval <- "model_evaluation/check3_modeval/"

if (!dir.exists(dir_modeval)) dir.create(dir_modeval)
# calculate all fits
cl <- makeCluster(getOption("cl.cores", detectCores() - cores_not_to_use))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_fits_bot <- parApply(cl, all_sumstats, 1, abc::gfit, sumstat = sims_stats, 
                         nb.replicate = cv_rep, tol = tol, subset = models == "bot")

all_fits_neut <- parApply(cl, all_sumstats,  1, abc::gfit, sumstat = sims_stats, 
                          nb.replicate = cv_rep, tol = tol, subset = models == "neut")
stopCluster(cl)

all_names <- names(all_seals)

# goodness of fit plots
pdf(file = paste0(dir_modeval, sim_name, "_goodnessoffit_noAR.pdf"), width = 18, height = 14)
par(mfrow = c(5, 6), mar=c(4,4,1,1))
sapply(1:length(all_names), function(x) {
  # check whether there are NAs, and if yes, exchange with mean
  if (any(is.na(all_fits_bot[[x]]$dist.sim))) {
    location <- which(is.na(all_fits_bot[[x]]$dist.sim))
    all_fits_bot[[x]]$dist.sim[location] <- mean(all_fits_bot[[x]]$dist.sim, na.rm = TRUE)
  }
  plot(all_fits_bot[[x]], main = paste0(all_names[x], "_bot"))
  
  if (any(is.na(all_fits_neut[[x]]$dist.sim))) {
    location <- which(is.na(all_fits_neut[[x]]$dist.sim))
    all_fits_neut[[x]]$dist.sim[location] <- mean(all_fits_neut[[x]]$dist.sim, na.rm = TRUE)
  }
  plot(all_fits_neut[[x]], main = paste0(all_names[x], "_neut"))
})
dev.off()


# save the p_values and distances
#if (!dir.exists("results/goodnessoffit_p")) dir.create("results/goodnessoffit_p")

all_p_bot <- unlist(lapply(all_fits_bot, function(x) summary(x)$pvalue))
all_p_neut <-  unlist(lapply(all_fits_neut, function(x) summary(x)$pvalue))

p_df <- data.frame(species = all_names, bot_p = all_p_bot, neut_p = all_p_neut)
write.table(p_df, file = paste0(dir_modeval, sim_name, "_p_vals_fit_noAR.txt"), row.names = FALSE)


# PCAs
# ?gfitpca
# sims_minus_NA <- sims_stats[-which(rowSums(is.na(as.matrix(sims_stats))) > 0), ]
# models_minus_NA <- models[-which(rowSums(is.na(as.matrix(sims_stats))) > 0)]

# 
# pdf(file = paste0("plots/", pop_size, "_pca.pdf"), width = 18, height = 14)
# par(mfrow = c(5, 2), mar=c(3,3,1,1))
# sapply(1:length(all_names), function(x) {
#   gfitpca(all_sumstats[x, ], sims_minus_NA, index = models_minus_NA, cprob = 0.05, main = all_names[x])
# })
# dev.off()



