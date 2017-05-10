### COMPARATIVE ANALYSIS WITH ABC BASED ON SIMULATIONS AND EMPIRICAL DATA
## This script will
# (1) plot the different summary statistics as boxplots for all models
# (2) calculate the probabilities of the models for all species
# (3) calculate the fit of the empirical data to the model
# (4) output all plots / txt files under plots and results

# load packages ------------
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


# preparations --------

# name of the simulation (without the .txt), defines the names of all output files too
sim_name <- "onepopprior_500k_gamma_varsamp"

# parameter definition model selection (and later abc?)
# tolerance rate
tol <- 0.001
# cross-validation replicates / number of replicates used to estimate the null distribution of the goodness-of-fit statistic
cv_rep <- 2
# reps for confusion matrix (see ?cv4postpr)
nval_conf_matrix <- 5
# method for model selection with approximate bayesian computation, see ?postpr
# "rejection", "mnlogistic", "neuralnet"
method <- "neuralnet"
# leave how many cores free?
cores_not_to_use <- 10 

# load empirical microsatellites and calc sumstats ----------- 

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28] # 

# calculate summary statistics by population 
all_sumstats_full <- lapply(all_seals_full, function(x) mssumstats(x, datatype = "microsats", by_pop = "cluster",
                                                                   start_geno = 4, mratio = "loose"))

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

# load reference table, i.e. the simulations ----------
# sim_name defined above name of all the saved files during the analysis
path_to_sims <- paste0(sim_name, ".txt")
sims <- as.data.frame(fread(path_to_sims, stringsAsFactors = FALSE))

# choose summary statistics for abc ---------

names(all_sumstats_full)
sumstats <- c("num_alleles_mean", "num_alleles_sd",
              "prop_low_afs_mean", "prop_low_afs_sd",
               "exp_het_mean", "mean_allele_range", "mratio_mean") # , 

all_sumstats <- all_sumstats_full[sumstats]

# prepare simulation data for abc ---------

# parameter columns in simulation data.frame
first_ss <- which(names(sims) == "num_alleles_mean")
params <- c(1:(first_ss-1))
# create a character vector with models
models <- sims$model
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]

# (1) first visual checks -----------
if (!dir.exists("plots/model_prob_plots")) dir.create("plots/model_prob_plots")

# check whether sumstats are different across models
pdf(paste0("plots/model_prob_plots/", sim_name, ".pdf"))
par(mfcol = c(3, 7), mar=c(4,4,1,1))
for (i in sumstats){
  boxplot(sims[[i]] ~ models, main = i)
}
dev.off() #only 129kb in size

# (2) can abc at all distinguish between the 4 models ? ----------
cv.modsel <- cv4postpr(models, sims_stats, nval=nval_conf_matrix, tol=tol, method=method)
s <- summary(cv.modsel)
png(paste0("plots/model_prob_plots/", sim_name, "confusion_mat.png"), width=4, height=4, units="in", res=300)
plot(cv.modsel, names.arg= model_names)
dev.off() #only 129kb in size

# (3) model selection with postpr() ------------

if (!dir.exists("results/model_probs")) dir.create("results/model_probs")
# posterior model probablities for all species
cl <- parallel::makeCluster(getOption("cl.cores", detectCores()-cores_not_to_use))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_probs <- parallel::parApply(cl, all_sumstats, 1, abc::postpr, index = models, 
                                sumstat = sims_stats, tol = tol, method = method)
stopCluster(cl)

all_probs_df <- do.call(rbind, lapply(all_probs, function(x) {
  out <- round(x$pred, 3)
  out
}))
write.table(all_probs_df, file = paste0("results/model_probs/",sim_name, "_model_selection.txt"))

# (4) Does the preferred model provide a good fit to the data?------------

# extremely time intense!
# for all species, check the fit for
if (!dir.exists("plots/goodnessoffit")) dir.create("plots/goodnessoffit")
# calculate all fits
cl <- makeCluster(getOption("cl.cores", detectCores()-cores_not_to_use))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_fits_bot <- parApply(cl, all_sumstats, 1, abc::gfit, sumstat = sims_stats, nb.replicate = cv_rep, tol = tol, subset = models == "bot")
all_fits_neut <- parApply(cl, all_sumstats,  1, abc::gfit, sumstat = sims_stats, nb.replicate = cv_rep, tol = tol, subset = models == "neut")
stopCluster(cl)

all_names <- names(all_seals_full)

# goodness of fit plots
pdf(file = paste0("plots/goodnessoffit/", sim_name, "_goodnessoffit.pdf"), width = 18, height = 14)
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
if (!dir.exists("results/goodnessoffit_p")) dir.create("results/goodnessoffit_p")

all_p_bot <- unlist(lapply(all_fits_bot, function(x) summary(x)$pvalue))
all_p_neut <-  unlist(lapply(all_fits_neut, function(x) summary(x)$pvalue))

p_df <- data.frame(species = all_names, bot_p = all_p_bot, neut_p = all_p_neut)
write.table(p_df, file = paste0("results/goodnessoffit_p/",sim_name, "_p_vals_fit.txt"), row.names = FALSE)

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

# (5) abc analysis preparation  ---------

# just use bottleneck model
mod <- "bot"

# subset parameters and summary statistics from simulations under the bottleneck model
stat_mod <- subset(sims_stats, subset=models==mod)
par_mod <- subset(sims_param, subset=models==mod)


# check whether a parameter can be estimated at all by looking at the error (optional here as very time intense)

# before inference, we see whether a parameter can be estimated at all
# parameter <- "N_bot"
#  
# cv_nbot <- function(method, nval = 50, tols = tol){
#    cv_res_rej <- cv4abc(data.frame(Na=par_mod[, parameter]), stat_mod, nval=nval,
#                         tols=tols, method=method, statistic = "mean")
#  }
# 
# all_cv_nbot <- lapply(c("ridge"), cv_nbot, nval = 5)

# (6) abc analysis  ---------

## transform model parameters to absolute values
par_mod <- par_mod %>% 
  mutate(N_bot = N0 * N_bot) %>%
  mutate(N_hist_bot = N0 * N_hist_bot) %>%
  mutate(start_bot = 4 * N0 * start_bot) %>%
  mutate(end_bot = 4 * N0 * end_bot) 

## abc method choice, all three possible
all_methods <- c("neuralnet", "ridge", "loclinear") # , "neuralnet"

# extract species names
all_species <- row.names(all_sumstats)
# parameters to estimate posteriors
all_parameters <- c("N_bot", "N0", "mu", "start_bot", "end_bot", "N_hist_bot", "sigma2_g", "p_single") 

# get all combinations of method, species and parameters
all_args <- expand.grid(all_methods,all_species, all_parameters)
all_args <- data.frame(apply(all_args, 2, as.character), stringsAsFactors = FALSE)
names(all_args) <- c("methods", "species", "pars")

## run abc
cl <- parallel::makeCluster(getOption("cl.cores", detectCores()-cores_not_to_use))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_probs <- parallel::parApply(cl, all_args, 1, abc::abc, target = all_sumstats[x["species"], ], 
                                param = par_mod[, x["pars"]], sumstat = stat_mod, tol = tol, method=x["methods"])
stopCluster(cl)


## run abc
# abc_est <- apply(all_args, 1, function(x) {
#   abc(target = all_sumstats[x["species"], ], param = par_mod[, x["pars"]], 
#       sumstat = stat_mod, tol = tol, method=x["methods"])
# })

# create a list of 3. 
# The first element ist the parameter data.frame for the abc. 
# The second element are the corresponding abc objects.
# the third element are the prior values

abc_full <- list(all_args, abc_est, par_mod)
save(abc_full, file = "abc_estimates/abc_full_gamma_all_methods.RData")







