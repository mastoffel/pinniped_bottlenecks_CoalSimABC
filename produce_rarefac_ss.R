# the script runs an abc analysis with one set of simulations on 
# multiple emprical microsatellite dataframes

# (1) definition of parameter values and wether to run a cross-validation
# on each parameter via pseudo-observed datasets taken from the prior distribution
# (2) calculates summary statistics with the sealABC package and transforms
# them into a clean data.frame
# (3) user selects summary statistics
# (4) simulations are read from a text file in split up into parameters
# and summary statistics
# (5) user chooses one of the models to do abc on (will be generalized to all models soon)
# (6) cross-validation with cv4abc if specified in (1)
# (7) abc for all species where a given model had the highest probability as
# saved under results/model_probs. abc can be done with different methods
# on all species and paramaters and is saved under abc_estimates/

# load packages 
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


# parameter definition ---------------------------------------------------------
# path to simulation
sims_name <- "sims_5000k_large_bot"

# do a cross validation for the abc parameters? ?cv4abc
calc_cv_for_abc <- FALSE
# if yes, which method?
method_cv <- "loclinear"
# how many pseudo-observed datasets should be evaluated?
nval_cv <- 100
# which tolerance level(s)?
tols_cv <- c(0.0002)

# do abc?
abc_analysis <- TRUE
# tolerance level for abc
tol_abc <- 0.001
## abc method choice, all three possible
all_methods <- c("loclinear") # "ridge", "loclinear", "neuralnet"

# prepare empirical data -------------------------------------------------------

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28]

# calculate summary statistics
# cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - 20))
# clusterEvalQ(cl, c(library("sealABC")))
# all_sumstats_full <- 
#   parallel::parLapply(cl, all_seals_full, mssumstats, by_pop = NULL, start_geno = 4, mratio = "loose",
#                       rarefaction = TRUE, nresamp = 500, nind = 10, nloc = NULL)
# stopCluster(cl)

all_sumstats_full <- lapply(all_seals_full, mssumstats, by_pop = NULL, start_geno = 4, mratio = "loose",
          rarefaction = TRUE, nresamp = 100, nind = 10, nloc = NULL)


sum_per_clust <- function(mssumstats_output) {
  if (nrow(mssumstats_output) > 1) {
    out <- as.data.frame(t(apply(mssumstats_output, 2, mean)))
  } else
    out <- mssumstats_output
  out
}
all_sumstats_full <- lapply(all_sumstats_full, sum_per_clust)

# transform summary statistics to data frame
all_sumstats_full <- do.call(rbind, all_sumstats_full)

save(all_sumstats_full, file = "seal_ss_rarefaction16ind.RData")
