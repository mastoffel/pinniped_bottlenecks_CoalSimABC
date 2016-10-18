# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC", force = TRUE)
library(sealABC)
library(data.table)
library(reshape2)
library("abctools")
library(abc)
library(ggplot2)

sims <- fread("sims_simple_pop100k_sim10k.txt", stringsAsFactors = FALSE)
# sims <- as.data.frame(sims)

sims <- as.data.frame(sims)

# which stats to use
sumstats <- c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
              "prop_low_afs_mean", "prop_low_afs_mean")
models <- sims$model
##### first visual checks 
# check whether sumstats are different across models
par(mfcol = c(3, 4), mar=c(5,5,1,1))
for (i in sumstats){
  boxplot(sims[[i]] ~ models, main = i)
}



# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

# hawaiian monk seal
genotypes <- all_seals[[12]]
genotypes <- genotypes[4:ncol(genotypes)]

# calculate summary statistics
obs_stats <- mssumstats(genotypes, type = "microsats", data_type = "empirical")

# divide stats and parameters
sims_stats <- sims[13:(ncol(sims)-1)] # last column is model factor
sims_param <- sims[1:12]
models <- sims[["model"]]

?AS.select
# ASchoice <- AS.select(obs_stats, sims_param, sims_stats)
# ASchoice$best

# select best summary statistics ("mean_allele_size_var", "mean_allele_size_var",)
sims_stats <- sims_stats[c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
                           "prop_low_afs_mean", "prop_low_afs_mean")] # "mean_allele_size_var"
obs_stats <- obs_stats[c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
                         "prop_low_afs_mean", "prop_low_afs_mean")] # "mratio_mean",


# can abc distinguish between the 4 models ?
# sims_stats <- as.matrix(sims_stats)
# sims_stats[is.infinite(sims_stats)] <- NA
cv.modsel <- cv4postpr(models, sims_stats, nval=10, tol=.01, method="rejection")
s <- summary(cv.modsel)
plot(cv.modsel, names.arg=c("bot", "neut"))


# model probabilities
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.001, method = "mnlogistic")
summary(mod_prob)
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.001, method = "rejection")
summary(mod_prob)
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.05, method = "neuralnet")
summary(mod_prob)


# goodness of fit
sims_stats <- as.matrix(sims_stats)
# sims_stats <- sims_stats[-(rowSums(is.infinite(sims_stats))>0), ]

res_gfit_bott <- gfit(target=obs_stats, sumstat=sims_stats[models == "bot",], tol = 0.001, statistic=mean, nb.replicate=100)
plot(res_gfit_bott, main="Histogram under H0")
summary(res_gfit_bott)

res_gfit_neut <- gfit(target=obs_stats, sumstat= sims_stats[models == "neut",], statistic=mean,tol = 0.01,  nb.replicate=10)
plot(res_gfit_neut, main="Histogram under H0")
summary(res_gfit_neut)


# or PCA
sims_stats_new <- apply(sims_stats, 2, function(x){
  x[!is.finite(x)] <- mean(x[is.finite(x)], na.rm =TRUE)
  x
} )
gfitpca(target=obs_stats, sumstat=sims_stats, index=models, cprob=0.01)


# posterior predictive checks (Gelman) // critics: using the data twice. 
mylabels <- c("num_alleles_mean", "mratio_mean", "mratio_sd", "prop_low_afs_mean")
# mylabels <- c("num_alleles_mean", "mean_allele_range", "mratio_mean")
par(mfrow = c(2,2), mar=c(5,2,4,0))
for (i in mylabels){
  hist(sims_stats[models == "bot",i],breaks=40, xlab=i, main="")
  abline(v = obs_stats[i], col = 2)
}



# parameter inference under the specified model
mod <- "bot"

stat_bot <- subset(sims_stats, subset=models==mod)
head(stat_bot)

par_bot <- subset(sims_param, subset=models==mod)
head(par_bot)

# before inference, we see whether a parameter can be estimated at all
cv_res_rej <- cv4abc(data.frame(Na=par_bot[,"N_bot"]), stat_bot, nval=100,
                     tols=c(0.01), method="loclinear", statistic = "median")
summary(cv_res_rej) #should be as low as possible
plot(cv_res_rej, caption = "rejection")

# parameter inference
# some transformations to par_bot
library(dplyr)
par_bot <- par_bot %>% 
  mutate(N_bot = N0 * N_bot) %>%
  mutate(N_hist_bot = N0 * N_hist_bot) %>%
  mutate(start_bot = 4 * N0 * start_bot) %>%
  mutate(end_bot = 4 * N0 * end_bot)  


res <- abc(target = unlist(obs_stats), param = par_bot[, "end_bot"], 
           sumstat = stat_bot, tol = 0.01, method="loclinear")

summary(res)
hist(res, breaks = 100)

par(cex=.8)
plot(res, param=par_bot$N_bot)





