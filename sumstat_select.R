# summary statistics selection with abctools

# comparative analysis of abc output 1

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
library(abctools)

# load all seal data
seal_descriptives <- read_excel("../data/all_data_seals.xlsx")
# add factor for abundance --> effective population size considered to be a maximum of one fifth of the abundance
seal_descriptives %<>% mutate(abund_level = ifelse(Abundance < 20000, "5k", ifelse(Abundance < 300000, "50k", "500k")))


###### prepare data #######

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28]

# minimum number of loci?
min(unlist(lapply(all_seals_full, function(x) (ncol(x)-3)/2)))
min(unlist(lapply(all_seals_full, function(x) nrow(x))))

# calculate summary statistics
all_sumstats_full <- lapply(all_seals_full, function(x) mssumstats(x[4:ncol(x)], type = "microsats", data_type = "empirical"))

# as data.frame
all_sumstats_full <- do.call(rbind, all_sumstats_full)


?AS.select


out <- AS.select(all_sumstats_full[1, ])

pop_size <- "500k"

# load simulations
path_to_sims <- paste0("sims_pop", pop_size, "_sim300k_smallsamp.txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# subset empirical summary statistics with the relevant populations
all_sumstats <- all_sumstats_full[seal_descriptives$abund_level == pop_size, ]

sumstat_afs <- all_sumstats_full[1, ]

# parameter columns
params <- c(1:12)

sims_params <- sims[1:300000, params]
sims_stats <- sims[1:300000, -c(params, ncol(sims))]

out <- AS.select(as.matrix(all_sumstats), as.matrix(sims_params), as.matrix(sims_stats))


for (pop_size in c( "5k", "50k", "500k")){ #
  
  # load simulations
  path_to_sims <- paste0("sims_pop", pop_size, "_sim300k_smallsamp.txt")
  
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
  tol <- 0.001
  # cross-validation replicates
  cv_rep <- 2
  # method
  method <- "neuralnet"
  
  
  # some processing
  # extract names of all models
  all_models <- names(table(models))
  
  # divide stats and parameters
  sims_stats <- sims[sumstats] 
  sims_param <- sims[params]


