# comparative analysis of abc output 1

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

#### load all seal descriptive data and add a maximum effect population size variable ####
seal_descriptives <- read_excel("../data/all_data_seals.xlsx")

###### prepare data ########
# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28]

cl <- parallel::makeCluster(getOption("cl.cores", detectCores()-20))
clusterEvalQ(cl, c(library("sealABC")))
all_sumstats_full <- parallel::parLapply(cl, all_seals_full, 
                                         function(x) mssumstats(x, by_pop = NULL, start_geno = 4, mratio = "loose",
                                                                rarefaction = TRUE, nresamp = 1000, nind = 20, nloc = 5))
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

# select summary statistics for posteriors. Should be the same as in comp_abc_1.R
# sumstats <- c("num_alleles_mean", "prop_low_afs_mean",   
#               "mean_allele_range",  "mean_allele_size_var",
#               "exp_het_mean")

sumstats <- c("num_alleles_mean", "num_alleles_sd",
              #"exp_het_mean", "exp_het_sd",
              "mean_allele_size_sd", "sd_allele_size_sd",
               "mean_allele_range", "sd_allele_range",
              "mratio_mean" , "prop_low_afs_mean") 

all_sumstats_full <- all_sumstats_full[sumstats]


# load sims
path_to_sims <- paste0("sims_simcoal1000k.txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

  
### subsetting and definition of model selection parameters
  
# parameter columns in simulation data.frame
param_start <- which(names(sims) == "sample_size")
param_end <- which(names(sims) == "range_constraint")
params <- c(param_start:param_end)
# create a character vector with models
models <- sims$model
# tolerance rate
tol <- 0.0005
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]
  
  
### run the actual abc analysis
  
  ## for loop for running abc for populations under different models
  #for (mod in c("bot")){
  
  # just use bottleneck model
  mod <- "bot"
  
  # subset parameters and summary statistics from simulations under the bottleneck model
  stat_mod <- subset(sims_stats, subset=models==mod)
  par_mod <- subset(sims_param, subset=models==mod)
  
  ## check whether a parameter can be estimated at all by looking at the error (optional here as very time intense)
  
  # if (mod == "bot"){
  
  # # before inference, we see whether a parameter can be estimated at all
  # cv_nbot <- function(method, nval = 50, tols = c(0.001, 0.0001, 0.00005)){
  #   cv_res_rej <- cv4abc(data.frame(Na=par_mod[,"nbot"]), stat_mod, nval=nval,
  #                        tols=tols, method=method)
  # }
  # 
  # all_cv_nbot <- lapply(c("neuralnet"), cv_nbot, nval = 50)
  # assign(paste0("cv_nbot_", pop_size), all_cv_nbot)
  # 
  
  ## abc method choice, all three possible
  all_methods <- c("ridge") # "ridge", "loclinear", "neuralnet"
  
  # extract species names
  all_species <- row.names(all_sumstats)
  # parameters to estimate posteriors
  all_parameters <- c("nbot", "pop_size", "mut_rate", "tbotstart", "tbotend", "nhist", "gsm_param") # "N0", "mu", "start_bot", "end_bot", "N_hist_bot", "sigma2_g"
  
  # get all combinations of method, species and parameters
  all_args <- expand.grid(all_methods,all_species, all_parameters)
  all_args <- data.frame(apply(all_args, 2, as.character), stringsAsFactors = FALSE)
  names(all_args) <- c("methods", "species", "pars")
  
  ## run abc
  abc_est <- apply(all_args, 1, function(x) {
    abc(target = all_sumstats[x["species"], ], param = par_mod[, x["pars"]], 
        sumstat = stat_mod, tol = tol, method=x["methods"])
  })
  
  # create a list of 2. The first element ist the parameter data.frame for the abc. The second element are the corresponding abc objects.
  # assign(paste0("abc_", pop_size), list(all_args, abc_est))
  # saved_file_name <- paste(paste0("abc_", pop_size),"RData",sep=".")
  # save(list(all_args, abc_est), saved_file_name)
  
 abc_full <- list(all_args, abc_est)
 save(abc_full, file = "abc_estimates/abc_full_simcoal.RData")
  
  