# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library(abc)
library(ggplot2)

path_to_sims <- "sims_pop500k_500k.txt"

# load data 
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# which stats to use
names(sims)
sumstats <- c("num_alleles_mean", "num_alleles_sd", "mratio_mean", "mratio_sd",
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
genotypes <- all_seals[[1]]
genotypes <- genotypes[4:ncol(genotypes)]
# genotypes <- genotypes[4:ncol(genotypes)]

# calculate summary statistics
all_sumstats <- mssumstats(genotypes, type = "microsats", data_type = "empirical")
obs_stats <- all_sumstats[sumstats]

### model selection
mod_prob <- abc::postpr(obs_stats, models, sims_stats, tol = 0.001, method = "mnlogistic")
summary(mod_prob)


#### Does the preferred model provide a good fit to the data?  ######
model <- "bot"
# TEST 1
res_gfit_bot <- gfit(target = obs_stats, sumstat = sims_stats[models == model, ], nb.replicate = 10, tol = 0.01)
summary(res_gfit_bot)
plot(res_gfit_bot)

# TEST 2 PCA
# get NA rows
which_NA <- which(rowSums(is.na(as.matrix(sims_stats)))>0)
# for final run
gfitpca(target=obs_stats, sumstat=sims_stats[-which_NA, ], index=models[-which_NA], cprob=0.1)



# parameter inference under the specified model
mod <- "bot"

stat_bot <- subset(sims_stats, subset=models==mod)
head(stat_bot)

par_bot <- subset(sims_param, subset=models==mod)
head(par_bot)

# before inference, we see whether a parameter can be estimated at all
cv_res_rej <- cv4abc(data.frame(Na=par_bot[,"N_bot"]), stat_bot, nval=40,
                     tols=c(.001), method="neuralnet")
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

res <- abc(target = unlist(obs_stats), param = par_bot[, "p_single"], 
           sumstat = stat_bot, tol = 0.005, method="ridge")
options(scipen=999)

summary(res)
hist(res, breaks = 1000)

par(cex=.8)
plot(res, param=par_bot$N_bot)





