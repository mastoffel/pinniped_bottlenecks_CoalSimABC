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
seal_descriptives %<>% mutate(abund_level = ifelse(Abundance < 50000, "10k", ifelse(Abundance < 500000, "100k", "1000k")))


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


pop_size <- "100k"
# run loop for all simulations
#for (pop_size in c("10k", "100k", "1000k")){
  
  
  # load simulations
  path_to_sims <- paste0("sims_simple_pop", pop_size, "_sim500k.txt")
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
  tol <- 0.01
  # cross-validation replicates
  cv_rep <- 5
  # method
  method <- "mnlogistic"
  
  
  # some processing
  # extract names of all models
  all_models <- names(table(models))
  
  # divide stats and parameters
  sims_stats <- sims[sumstats] 
  sims_param <- sims[params]
  
  
  # parameter inference under the specified model
  mod <- "bot"
  
  stat_bot <- subset(sims_stats, subset=models==mod)
  # head(stat_bot)
  
  par_bot <- subset(sims_param, subset=models==mod)
  # head(par_bot)
  
  # before inference, we see whether a parameter can be estimated at all
  cv_res_rej <- cv4abc(data.frame(Na=par_bot[,"N_bot"]), stat_bot, nval=100,
                       tols=c(.005,.01, 0.05), method="loclinear")
  summary(cv_res_rej) #should be as low as possible
  plot(cv_res_rej, caption = "rejection")
  
  library(dplyr)
  par_bot <- par_bot %>% 
    mutate(N_bot = N0 * N_bot) %>%
    mutate(N_hist_bot = N0 * N_hist_bot) %>%
    mutate(start_bot = 4 * N0 * start_bot) %>%
    mutate(end_bot = 4 * N0 * end_bot) 
  
  # parameter inference
  # some transformations to par_bot
  # hawaiian monk seal
  genotypes <- all_seals[[1]]
  genotypes <- genotypes[4:ncol(genotypes)]
  
  # calculate summary statistics
  obs_stats <- mssumstats(genotypes, type = "microsats", data_type = "empirical")
  obs_stats <- obs_stats[sumstats]

  
  res <- abc(target = unlist(obs_stats), param = par_bot[, "N_hist_bot"], 
             sumstat = stat_bot, tol = 0.01, method="rejection")
  
  summary(res)
  hist(res, breaks = 100)
  
  par(cex=.8)
  plot(res, param=par_bot$N_bot)
  
  
  
  