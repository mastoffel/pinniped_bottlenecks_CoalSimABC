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
sims_name <- "sims_2000k_optimal"

# do a cross validation for the abc parameters? ?cv4abc
calc_cv_for_abc <- TRUE
# if yes, which method?
method_cv <- "loclinear"
# how many pseudo-observed datasets should be evaluated?
nval_cv <- 50
# which tolerance level(s)?
tols_cv <- c(0.005)

# tolerance level for abc
tol_abc <- 0.005
## abc method choice, all three possible
all_methods <- c("ridge", "loclinear", "neuralnet") # "ridge", "loclinear", "neuralnet"

# prepare empirical data -------------------------------------------------------

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28]

# calculate summary statistics
cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - 20))
clusterEvalQ(cl, c(library("sealABC")))
all_sumstats_full <- 
  parallel::parLapply(cl, all_seals_full, mssumstats, by_pop = NULL, start_geno = 4, mratio = "loose",
                                              rarefaction = TRUE, nresamp = 1000, nind = 30, nloc = 5)
stopCluster(cl)

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


# select summary statistics for posteriors. ------------------------------------

sumstats <- c("num_alleles_mean", "num_alleles_sd",
              "exp_het_mean", "mratio_mean", "prop_low_afs_mean")
all_sumstats_full <- all_sumstats_full[sumstats]


# load simulations -------------------------------------------------------------
path_to_sims <- paste0(sims_name,".txt")
sims <- fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# create dataframes with parameters and summary statistics ---------------------

# parameter columns in simulation data.frame
param_start <- which(names(sims) == "sample_size")
param_end <- which(names(sims) == "range_constraint")
params <- c(param_start:param_end)
# create a character vector with models
models <- sims$model
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]
  

# just use bottleneck model
mod <- "bot"

# subset sims_stats and sims_param for the given model
stat_mod <- subset(sims_stats, subset = models == mod)
par_mod <- subset(sims_param, subset = models == mod)

# check whether a parameter can be estimated at all ----------------------------
# (optional here as very time intense)

if (calc_cv_for_abc == TRUE) {
  # before inference, we see whether a parameter can be estimated at all
  cv_nbot <- function(method, pars, nval = 5, tols = c(0.001, 0.005, 0.0005)){ #
    cv_res_rej <- cv4abc(data.frame(par_mod)[pars], stat_mod, nval = nval,
                         tols = tols, method = method)
  }
  pars <- c("nbot", "nhist", "tbotend", "tbotstart", "mut_rate", "gsm_param")
  
  cv_res <- cv_nbot(method = method_cv, pars = pars, nval = nval_cv,  tols = tols_cv)
  out <- paste0("abc_estimates/cv_params_", sims_name, ".txt")
  write.table(cv_res, file = out, row.names = FALSE)
}

# run the actual abc analysis --------------------------------------------------

## load model probabilities
model_probs <- read.table(paste0("results/model_probs/", 
                                 sims_name, "_model_selection.txt"))
model_bot <- model_probs$bot > 0.5

# extract species names for species that have been bottlenecked according
# to the model selection
all_species <- row.names(all_sumstats_full)[model_bot]

# parameters to estimate posteriors
all_parameters <- c("nbot", "pop_size", "mut_rate", "tbotstart", 
                    "tbotend", "nhist", "gsm_param") #

# get all combinations of method, species and parameters
all_args <- expand.grid(all_methods, all_species, all_parameters)
all_args <- data.frame(apply(all_args, 2, as.character), stringsAsFactors = FALSE)
names(all_args) <- c("methods", "species", "pars")

## run abc
abc_est <- apply(all_args, 1, function(x) {
  abc(target = all_sumstats[x["species"], ], param = par_mod[, x["pars"]], 
      sumstat = stat_mod, tol = tol_abc, method = x["methods"])
})

# create a list of 2. The first element ist the parameter data.frame for the abc. 
# The second element are the corresponding abc objects.

abc_full <- list(all_args, abc_est)
save(abc_full, file = "abc_estimates/abc_", sims_name, ".RData")

