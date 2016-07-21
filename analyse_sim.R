# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)

sims <- fread("sims_full_broad_priors.txt", stringsAsFactors = FALSE)
sims <- as.data.frame(sims)
# sims <- read.table("sims_full_broad_priors.txt", stringsAsFactors = FALSE)
# names(sims) <- unlist(sims[1, ])
# sims <- sims[-1, ]
# transform everything to numeric
# sims <- cbind(do.call(cbind, lapply(sims[-ncol(sims)], as.numeric)), sims[ncol(sims)])

#
par(mfcol = c(3,4))
boxplot(sims$num_alleles_mean ~ sims$model)
boxplot(sims$mean_allele_size_sd ~ sims$model)
boxplot(sims$mean_allele_range ~ sims$model)
boxplot(sims$exp_het_mean ~ sims$model)
boxplot(sims$obs_het_mean ~ sims$model)
boxplot(sims$mratio_mean ~ sims$model)

library(ggplot2)

par(mfcol = c(3,3))
boxplot(sims$num_alleles_sd ~ sims$model)
boxplot(sims$sd_allele_size_sd ~ sims$model)
boxplot(sims$sd_allele_range ~ sims$model)
boxplot(sims$exp_het_sd ~ sims$model)
boxplot(sims$obs_het_sd ~ sims$model)
boxplot(sims$mratio_sd ~ sims$model)


# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

lapply(all_seals, function(x) range(unlist(x[4:ncol(x)]), na.rm = TRUE))

# hawaiian monk seal
genotypes <- all_seals[[33]]
genotypes <- genotypes[4:ncol(genotypes)]

# calculate summary statistics
obs_stats <- mssumstats(genotypes, type = "microsats")

# divide stats and parameters
sims_stats <- sims[14:(ncol(sims)-1)] # last column is model factor
sims_param <- sims[1:13]
models <- sims[["model"]]

library(abc)

# select best summary statistics
sims_stats <- sims_stats[c("num_alleles_mean", "mean_allele_range", "mratio_mean","mean_allele_size_sd")]
obs_stats <- obs_stats[c("num_alleles_mean", "mean_allele_range", "mratio_mean", "mean_allele_size_sd")]


# can abc distinguish between the 4 models ?
cv.modsel <- cv4postpr(models, sims_stats, nval=5, tol=.01, method="rejection")
s <- summary(cv.modsel)
plot(cv.modsel, names.arg=c("bot", "neut", "decl", "exp"))


# model probabilities
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.001, method = "mnlogistic")
summary(mod_prob)
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.1, method = "rejection")
summary(mod_prob)
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.01, method = "neuralnet")
summary(mod_prob)


# goodness of fit
sims_stats <- as.matrix(sims_stats)

res_gfit_bott <- gfit(target=obs_stats, sumstat= sims_stats[models == "bot",], tol = 0.01, statistic=median, nb.replicate=1000)
plot(res_gfit_bott, main="Histogram under H0")
summary(res_gfit_bott)

res_gfit_neut <- gfit(target=obs_stats, sumstat= sims_stats[models == "neut",], statistic=mean, nb.replicate=10)
plot(res_gfit_neut, main="Histogram under H0")
summary(res_gfit_neut)

res_gfit_decl <- gfit(target=obs_stats, sumstat= sims_stats[models == "decl",], statistic=mean, nb.replicate=100)
plot(res_gfit_decl, main="Histogram under H0")
summary(res_gfit_decl)

res_gfit_exp <- gfit(target=obs_stats, sumstat= sims_stats[models == "exp",], statistic=mean, nb.replicate=100)
plot(res_gfit_exp, main="Histogram under H0")
summary(res_gfit_exp)

# or PCA
sims_stats_new <- apply(sims_stats, 2, function(x){
  x[!is.finite(x)] <- mean(x[is.finite(x)], na.rm =TRUE)
  x
} )
gfitpca(target=obs_stats, sumstat=sims_stats_new, index=models, cprob=.01)


# posterior predictive checks (Gelman) // critics: using the data twice. 
mylabels <- c("num_alleles_mean","mean_allele_size_sd", "mean_allele_range", "exp_het_mean", "obs_het_mean",
              "mratio_mean")
mylabels <- c("num_alleles_mean", "mean_allele_range", "mratio_mean")
par(mfrow = c(2,3), mar=c(5,2,4,0))
for (i in mylabels){
  hist(sims_stats[models == "exp",i],breaks=40, xlab=i, main="")
  abline(v = obs_stats[i], col = 2)
}



# parameter inference under the bottleneck model

stat_bot <- subset(sims_stats, subset=models=="exp")
head(stat_bot)

par_bot <- subset(sims_param, subset=models=="exp")
head(par_bot)

# before inference, we see whether a parameter can be estimated at all

cv_res_rej <- cv4abc(data.frame(Na=par_bot[,"N_bot"]), stat_bot, nval=100,
                     tols=c(.01,.05, .3), method="rejection", statistic = "mean")
summary(cv_res_rej) #should be as low as possible
plot(cv_res_rej, caption = "rejection")

# parameter inference
# some transformations to par_bot
library(dplyr)
par_bot <- par_bot %>% 
  mutate(N_bot = N0 * N_bot) %>%
  mutate(N_hist = N0 * N_hist) %>%
  mutate(N_hist_decl = N0 * N_hist_decl) %>%
  mutate(N_hist_exp = N0 * N_hist_exp) %>%
  mutate(start_bot = 4 * N0 * start_bot) %>%
  mutate(end_bot = 4 * N0 * end_bot)  %>%
  mutate(start_decl = 4 * N0 * start_decl) %>%
  mutate(start_exp = 4 * N0 * start_exp)

res <- abc(target = unlist(obs_stats), param = par_bot[, "start_exp"], 
           sumstat = stat_bot, tol = 0.01, method="rejection")

summary(res)
hist(res, breaks = 100)

par(cex=.8)
plot(res, param=par_bot$start_bot)





