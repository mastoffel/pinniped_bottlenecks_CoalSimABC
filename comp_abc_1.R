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

seal_descriptives <- read_excel("../data/all_data_seals.xlsx")
# add factor for abundance --> effective population size considered to be a maximum of one fifth of the abundance
seal_descriptives %<>% mutate(abund_level = ifelse(Abundance < 20000, "5k", ifelse(Abundance < 300000, "50k", "500k")))


###### prepare data

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28]

# calculate summary statistics
all_sumstats_full <- lapply(all_seals_full, function(x) mssumstats(x[4:ncol(x)], type = "microsats", data_type = "empirical"))
# as data.frame
all_sumstats_full <- do.call(rbind, all_sumstats_full)
# add factor for abundance level
# all_sumstats_full$abund_level <- seal_descriptives$abund_level

# which ss to use
# names(sims)
sumstats <- c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
              "prop_low_afs_mean", "prop_low_afs_sd")

all_sumstats_full <- all_sumstats_full[sumstats]

# subset all under 50k (effective population size maximum of 10k)
# pop_size <- "10k" # alternative 100k, 1000k

# pop_size <- "100k"
# run loop for all simulations

for (pop_size in c("5k", "50k", "500k")){ #

# load simulations
path_to_sims <- paste0("sims_simple_pop", pop_size, "_sim100k.txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# subset empirical summary statistics with the relevant populations
all_sumstats <- all_sumstats_full[seal_descriptives$abund_level == pop_size, ]

# subset all_seals with the relevant populations
all_seals <- all_seals_full[seal_descriptives$abund_level == pop_size]

# parameter columns
params <- c(1:12)
# character vector with models
models <- sims$model
# tolerance rate
tol <- 0.005
# cross-validation replicates
cv_rep <- 30
# method
method <- "mnlogistic"


# some processing
# extract names of all models
all_models <- names(table(models))

# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]


##### (1) first visual checks 
# check whether sumstats are different across models
pdf(paste0("plots/", pop_size, "_boxplots.pdf"))
par(mfcol = c(3, 3), mar=c(4,4,1,1))
for (i in sumstats){
  boxplot(sims[[i]] ~ models, main = i)
}
dev.off() #only 129kb in size


### (2) can abc at all distinguish between the 4 models ?
cv.modsel <- cv4postpr(models, sims_stats, nval=cv_rep, tol=tol, method=method)
s <- summary(cv.modsel)
png(paste0("plots/", pop_size, "_confusion_matrix.png"), width=4, height=4, units="in", res=300)
plot(cv.modsel, names.arg= all_models)
dev.off() #only 129kb in size



### (3) model selection

# abc_mod_probs <- function(obs_stats, models, sims_stats, tol, method, sumstats){
#   # model probabilities
#   mod_prob <- abc::postpr(obs_stats, models, sims_stats, tol = tol, method = method)
#   # sum_prob <- summary(mod_prob)
# }

#check probabilites for all species
cl <- makeCluster(getOption("cl.cores", detectCores()-5))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_probs <- parallel::parApply(cl, all_sumstats, 1, abc::postpr, index = models, sumstat = sims_stats, tol = tol, method = method)
stopCluster(cl)

all_probs_df <- do.call(rbind, lapply(all_probs, function(x) {
  out <- round(x$pred, 3)
  out
} ))
write.table(all_probs_df, file = paste0("output/", pop_size, "_model_selection.txt"))



#### Does the preferred model provide a good fit to the data?  
# calc_fit <- function(obs_stats, model, sumstats, all_sumstats, sims_stats) {
#   res_gfit_bot <- abc::gfit(target = all_sumstats[species_name, sumstats], sumstat = sims_stats[models == model, ], 
#                        nb.replicate = cv_rep, tol = tol)
# }

# calculate all fits
cl <- makeCluster(getOption("cl.cores", detectCores()-5))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_fits_bot <- parApply(cl, all_sumstats, 1, abc::gfit, sumstat = sims_stats, nb.replicate = cv_rep, tol = tol, subset = models == "bot")
all_fits_neut <- parApply(cl, all_sumstats,  1, abc::gfit, sumstat = sims_stats, nb.replicate = cv_rep, tol = tol, subset = models == "neut")
stopCluster(cl)

all_names <- names(all_seals)

# plots
pdf(file = paste0("plots/", pop_size, "_goodnessoffit.pdf"), width = 18, height = 14)
par(mfrow = c(5, 4), mar=c(4,4,1,1))
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
# summary(res_gfit_bot)$pvalue

all_p_bot <- unlist(lapply(all_fits_bot, function(x) summary(x)$pvalue))
all_p_neut <-  unlist(lapply(all_fits_neut, function(x) summary(x)$pvalue))

p_df <- data.frame(species = all_names, bot_p = all_p_bot, neut_p = all_p_neut)
write.table(p_df, file = paste0("output/", pop_size, "_p_vals_fit.txt"), row.names = FALSE)


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



} # end loop