# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library("abctools")
library(abc)
library(ggplot2)

path_to_sims <- "sims/sims_100000_6mod_lowmu.txt"

# sims <- rbindlist(list(sims, sims2, sims3), use.names = TRUE, fill = TRUE)

sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# which stats to use
sumstats <- c("num_alleles_mean", "mean_allele_range", "exp_het_mean", "mratio_mean")
# parameter columns
params <- c(1:18)
# character vector with models
models <- sims$model
# tolerance rate
tol <- 0.01
# cross-validation replicates
cv_rep <- 5
# method
method <- "mnlogistic"


# check whether sumstats are different across models
par(mfcol = c(length(sumstats), 1), mar=c(5,5,3,1))
for (i in sumstats){
  boxplot(sims[[i]] ~ models, main = i)
}

# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")


# seal species
genotypes <- all_seals[[19]]
genotypes <- genotypes[4:ncol(genotypes)]


abc_analysis <- function(genotypes, sims = sims, params = params, sumstats  = sumstats , models = models, 
                         tol = 0.01, cv_rep = 5, method = "mnlogistic", conf_matrix = FALSE) {
  
  # extract names of all models
  all_models <- names(table(models))
  
  # calculate empirical summary statistics
  obs_stats <- sealABC::mssumstats(genotypes, type = "microsats")

  # divide stats and parameters
  sims_stats <- sims[sumstats] 
  obs_stats <- obs_stats[sumstats] 
  sims_param <- sims[params]

  # can abc distinguish between the 4 models ?
  # sims_stats <- as.matrix(sims_stats)
  # sims_stats[is.infinite(sims_stats)] <- NA
  if (conf_matrix) {
    cv.modsel <- cv4postpr(models, sims_stats, nval=cv_rep, tol=tol, method=method)
    s <- summary(cv.modsel)
    png("confusion_matrix.png", width=4, height=4, units="in", res=300)
    plot(cv.modsel, names.arg= all_models)
    dev.off() #only 129kb in size
  }
  
  # model probabilities
  mod_prob <- postpr(obs_stats, models, sims_stats, tol = tol, method = method)
  sum_prob <- summary(mod_prob)
  
  # goodness of fit
  e1 <- environment()
  calc_gfit <- function(model){
    res_gfit <- gfit(target=obs_stats, sumstat=sims_stats[models == model,], tol = tol, statistic=mean, nb.replicate=cv_rep)
    # assign(paste("res_gfit_", model, sep = ""), res_gfit, envir = e1)
  }
  all_gfits <- lapply(all_models, calc_gfit)
  names(all_gfits) <- all_models

# or PCA
# sims_stats_new <- apply(sims_stats, 2, function(x){
#   x[is.na(x)] <- mean(x, na.rm =TRUE)
#   x
# } )
# gfit_pca_new(target=obs_stats, sumstat=sims_stats_new, index=models, cprob=0.01)


# posterior predictive checks (Gelman) // critics: using the data twice. 
mylabels <- sumstats

png("ppc.png", width=8, height=8, units="in", res=600)
par(mfrow = c(6,4), mar=c(5,5,2,1))

for (k in all_models){
for (i in mylabels){
  hist(sims_stats[models == k, i],breaks=40, xlab=i, main=k)
  abline(v = obs_stats[i], col = 2)
}
}
dev.off()



# parameter inference under the specified model
mod <- "bot"

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

}



