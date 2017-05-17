# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library("abctools")
library(abc)
library(ggplot2)

path_to_sims <- "sims_simple_pop100k_sim500k.txt"

# load data 
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# which stats to use
names(sims)
sumstats <- c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
              "prop_low_afs_mean", "prop_low_afs_sd")
# parameter columns
params <- c(1:12)
# character vector with models
models <- sims$model
# tolerance rate
tol <- 0.001
# cross-validation replicates
cv_rep <- 10
# method
method <- "mnlogistic"

# some processing

# extract names of all models
all_models <- names(table(models))

# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]


##### first visual checks 
# check whether sumstats are different across models
par(mfcol = c(3, 3), mar=c(5,5,1,1))
for (i in sumstats){
  boxplot(sims[[i]] ~ models, main = i)
}

### can abc at all distinguish between the 4 models ?
cv.modsel <- cv4postpr(models, sims_stats, nval=cv_rep, tol=tol, method=method)
s <- summary(cv.modsel)
# png("confusion_matrix.png", width=4, height=4, units="in", res=300)
plot(cv.modsel, names.arg= all_models)
# dev.off() #only 129kb in size

# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

# seal species
# genotypes <- all_seals[[19]]
# genotypes <- genotypes[4:ncol(genotypes)]

all_sumstats <- lapply(all_seals, function(x) mssumstats(x[4:ncol(x)], type = "microsats", data_type = "empirical"))
all_sumstats <- do.call(rbind, all_sumstats)



### model selection

abc_mod_probs <- function(genotypes, models, sims_stats, tol, method, sumstats){
  genotypes <- genotypes[4:ncol(genotypes)]
  # calculate empirical summary statistics
  obs_stats <- sealABC::mssumstats(genotypes, type = "microsats", data_type = "empirical")
  obs_stats <- obs_stats[sumstats]
  # model probabilities
  mod_prob <- postpr(obs_stats, models, sims_stats, tol = tol, method = "mnlogistic")
  # sum_prob <- summary(mod_prob)
}
#check probabilites for all species
all_probs <- lapply(all_seals[1:28], abc_mod_probs, models, sims_stats, tol, method, sumstats)
all_probs_df <- do.call(rbind, lapply(all_probs, function(x) {
                            sum_abc <- summary(x)
                            # out <- sum_abc$rejection$Prob
                            out <- sum_abc$mnlogistic$Prob
                           } ))
write.table(all_probs_df, file = "output/model_selection.txt")

#abc_analysis <- function(genotypes, sims = sims, params = params, sumstats  = sumstats , models = models, 
#                         tol = 0.01, cv_rep = 5, method = "mnlogistic", conf_matrix = FALSE) {

#### Does the preferred model provide a good fit to the data?  ######

# TEST 1
calc_fit <- function(species_name, model) {
    res_gfit_bot <- gfit(target = all_sumstats[species_name, sumstats], sumstat = sims_stats[models == model, ], 
                       nb.replicate = 100, tol = 0.01)
}

# get all names
all_names <- names(all_seals)[1:28]
# calculate all fits
all_fits_bot <- lapply(all_names, calc_fit, "bot")
all_fits_neut <- lapply(all_names, calc_fit, "neut")

# plots
pdf(file = "plots/goodnessoffit.pdf")
par(mfrow = c(8, 4), mar=c(3,3,1,1))
sapply(1:16, function(x) {
  plot(all_fits_bot[[x]], main = all_names[x])
  plot(all_fits_neut[[x]], main = all_names[x])
})
dev.off()

pdf(file = "plots/goodnessoffit2.pdf")
par(mfrow = c(6, 4), mar=c(3,3,1,1))
sapply(17:28, function(x) {
  plot(all_fits_bot[[x]], main = all_names[x])
  plot(all_fits_neut[[x]], main = all_names[x])
})
dev.off()


# save the p_values and distances
summary(res_gfit_bot)$pvalue

all_p_bot <- unlist(lapply(all_fits_bot, function(x) summary(x)$pvalue))
all_p_neut <-  unlist(lapply(all_fits_neut, function(x) summary(x)$pvalue))

p_df <- data.frame(species = all_names, bot_p = all_p_bot, neut_p = all_p_neut)
write.table(p_df, file = "output/p_vals_fit.txt", row.names = FALSE)


# TEST 2 PCA
obs_stats <- all_sumstats[1:28, sumstats]
# get NA rows
which_NA <- which(rowSums(is.na(as.matrix(sims_stats)))>0)
# for final run
# gfitpca(target=obs_stats["crabeater_seal", ], sumstat=sims_stats[-which_NA, ], index=models[-which_NA], cprob=0.3)






# parameter inference under the specified model
mod <- "bot"

stat_bot <- subset(sims_stats, subset=models==mod)
head(stat_bot)

par_bot <- subset(sims_param, subset=models==mod)
head(par_bot)

# before inference, we see whether a parameter can be estimated at all
cv_res_rej <- cv4abc(data.frame(Na=par_bot[,"N_bot"]), stat_bot, nval=100,
                     tols=c(.005,.01, 0.05), method="loclinear")
summary(cv_res_rej) #should be as low as possible
plot(cv_res_rej, caption = "rejection")

# parameter inference
# some transformations to par_bot
# hawaiian monk seal
genotypes <- all_seals[[21]]
genotypes <- genotypes[4:ncol(genotypes)]

# calculate summary statistics
obs_stats <- mssumstats(genotypes, type = "microsats", data_type = "empirical")
obs_stats <- obs_stats[sumstats]
library(dplyr)
par_bot <- par_bot %>% 
  mutate(N_bot = N0 * N_bot) %>%
  mutate(N_hist_bot = N0 * N_hist_bot) %>%
  mutate(start_bot = 4 * N0 * start_bot) %>%
  mutate(end_bot = 4 * N0 * end_bot) 

res <- abc(target = unlist(obs_stats), param = par_bot[, "N0"], 
           sumstat = stat_bot, tol = 0.01, method="loclinear")

summary(res)
hist(res, breaks = 100)

par(cex=.8)
plot(res, param=par_bot$N_bot)





