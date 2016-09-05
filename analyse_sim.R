# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library("abctools")
library(abc)
library(ggplot2)

sims <- fread("sims_50000_6mod.txt", stringsAsFactors = FALSE)
# sims <- as.data.frame(sims)

sims2 <- fread("sims_50000_6mod_nr2.txt", stringsAsFactors = FALSE)
# sims2 <- as.data.frame(sims2)

sims3 <- fread("sims_50000_6mod_nr3.txt", stringsAsFactors = FALSE)
# sims3 <- as.data.frame(sims2)

sims <- rbindlist(list(sims, sims2, sims3), use.names = TRUE, fill = TRUE)

sims <- as.data.frame(sims)

# sims_melt <- data.table::melt(sims, id.vars )

sims_melt <- melt(sims[c(15,17,28)], na.rm = TRUE, id.var = "model")
ggplot(data = sims_melt, aes(x = variable, y = value)) + 
        geom_boxplot(aes(x = model)) +
        facet_wrap(~ variable)
        

par(mfcol = c(4,4))

boxplot(sims$num_alleles_mean ~ sims$model)
boxplot(sims$mean_allele_size_var ~ sims$model)
boxplot(sims$mean_allele_range ~ sims$model)
boxplot(sims$exp_het_mean ~ sims$model)
boxplot(sims$obs_het_mean ~ sims$model)
boxplot(sims$mratio_mean ~ sims$model)
# boxplot(sims$g2 ~ sims$model)
boxplot(sims$het_excess ~ sims$model)


par(mfcol = c(3,3))
boxplot(sims$num_alleles_sd ~ sims$model)
boxplot(sims$sd_allele_size_var ~ sims$model)
boxplot(sims$sd_allele_range ~ sims$model)
boxplot(sims$exp_het_sd ~ sims$model)
boxplot(sims$obs_het_sd ~ sims$model)
boxplot(sims$mratio_sd ~ sims$model)


# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

lapply(all_seals, function(x) range(unlist(x[4:ncol(x)]), na.rm = TRUE))

# hawaiian monk seal
genotypes <- all_seals[[13]]
genotypes <- genotypes[4:ncol(genotypes)]

# calculate summary statistics
obs_stats <- mssumstats(genotypes, type = "microsats")

# divide stats and parameters
sims_stats <- sims[19:(ncol(sims)-1)] # last column is model factor
sims_param <- sims[1:18]
models <- sims[["model"]]

?AS.select
# ASchoice <- AS.select(obs_stats, sims_param, sims_stats)
# ASchoice$best

# select best summary statistics
sims_stats <- sims_stats[c("num_alleles_mean", "mean_allele_range", "exp_het_mean", "mean_allele_size_var", "mratio_mean", "het_excess" )] # "mean_allele_size_var"
obs_stats <- obs_stats[c("num_alleles_mean", "mean_allele_range", "exp_het_mean","mean_allele_size_var", "mratio_mean", "het_excess" )] # "mratio_mean",


# can abc distinguish between the 4 models ?
# sims_stats <- as.matrix(sims_stats)
# sims_stats[is.infinite(sims_stats)] <- NA
cv.modsel <- cv4postpr(models, sims_stats, nval=100, tol=.01, method="rejection", subset = !(is.infinite(sims_stats$mratio_mean)))
s <- summary(cv.modsel)
plot(cv.modsel, names.arg=c("bot", "neut", "decl", "exp"))


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

res_gfit_bott <- gfit(target=obs_stats, sumstat=sims_stats[models == "bot",], tol = 0.01, statistic=mean, nb.replicate=10)
plot(res_gfit_bott, main="Histogram under H0")
summary(res_gfit_bott)

res_gfit_neut <- gfit(target=obs_stats, sumstat= sims_stats[models == "neut",], statistic=mean,tol = 0.01,  nb.replicate=10)
plot(res_gfit_neut, main="Histogram under H0")
summary(res_gfit_neut)

res_gfit_decl <- gfit(target=obs_stats, sumstat= sims_stats[models == "decl_IA",], statistic=mean,tol = 0.01,  nb.replicate=10)
plot(res_gfit_decl, main="Histogram under H0")
summary(res_gfit_decl)

res_gfit_exp <- gfit(target=obs_stats, sumstat= sims_stats[models == "exp_IA",], statistic=mean, tol = 0.01, nb.replicate=100)
plot(res_gfit_exp, main="Histogram under H0")
summary(res_gfit_exp)

# or PCA
sims_stats_new <- apply(sims_stats, 2, function(x){
  x[!is.finite(x)] <- mean(x[is.finite(x)], na.rm =TRUE)
  x
} )
gfitpca(target=obs_stats, sumstat=sims_stats_new, index=models, cprob=0.1)


# posterior predictive checks (Gelman) // critics: using the data twice. 
mylabels <- c("num_alleles_mean", "mean_allele_range", "mean_allele_size_var", "exp_het_mean", "mratio_mean", "het_excess")
# mylabels <- c("num_alleles_mean", "mean_allele_range", "mratio_mean")
par(mfrow = c(2,3), mar=c(5,2,4,0))
for (i in mylabels){
  hist(sims_stats[models == "exp_rec",i],breaks=40, xlab=i, main="")
  abline(v = obs_stats[i], col = 2)
}



# parameter inference under the specified model
mod <- "exp_rec"

stat_bot <- subset(sims_stats, subset=models==mod)
head(stat_bot)

par_bot <- subset(sims_param, subset=models==mod)
head(par_bot)

# before inference, we see whether a parameter can be estimated at all
cv_res_rej <- cv4abc(data.frame(Na=par_bot[,"N_bot"]), stat_bot, nval=100,
                     tols=c(.01), method="loclinear", statistic = "mean")
summary(cv_res_rej) #should be as low as possible
plot(cv_res_rej, caption = "rejection")

# parameter inference
# some transformations to par_bot
library(dplyr)
par_bot <- par_bot %>% 
  mutate(N_bot = N0 * N_bot) %>%
  mutate(N_hist_bot = N0 * N_hist_bot) %>%
  mutate(N_hist_decl = N0 * N_hist_decl) %>%
  mutate(N_hist_exp = N0 * N_hist_exp) %>%
  mutate(start_bot = 4 * N0 * start_bot) %>%
  mutate(end_bot = 4 * N0 * end_bot)  %>%
  mutate(start_decl_IA = 4 * N0 * start_decl_IA) %>%
  mutate(start_decl_recent = 4 * N0 * start_decl_recent) %>%
  mutate(start_exp_IA = 4 * N0 * start_exp_IA) %>%
  mutate(start_exp_recent = 4 * N0 * start_exp_recent)

res <- abc(target = unlist(obs_stats), param = par_bot[, "start_exp_recent"], 
           sumstat = stat_bot, tol = 0.01, method="loclinear")

summary(res)
hist(res, breaks = 100)

par(cex=.8)
plot(res, param=par_bot$N_bot)





