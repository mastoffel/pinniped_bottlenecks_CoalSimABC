# ABC analysis part 2: Parameter estimation

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
library(readr)

# parameter definition ---------------------------------------------------------
# path to simulation
sims_name <- "sims_10000k"

# do a cross validation for the abc parameters? ?cv4abc
calc_cv_for_abc <- FALSE
# if yes, which method?
method_cv <- "rejection"      # "rejection"    # "loclinear"
# how many pseudo-observed datasets should be evaluated?
nval_cv <- 1000
# which tolerance level(s)?
tols_cv <- c(0.0005)
# parallel
run_parallel <- TRUE

# do abc?
abc_analysis <- FALSE
# tolerance level for abc
tol_abc <- 0.0005
## abc method choice, all three possible
all_methods <- c("loclinear") # "ridge", "loclinear", "neuralnet"

# prepare empirical data -------------------------------------------------------

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop_29.xlsx")[1:29]

# get summary stats data
all_sumstats_full <- read_delim("../data/all_sumstats_40ind_29.txt", delim = " ")

# select summary statistics for posteriors. ------------------------------------

sumstats <- c("num_alleles_mean", 
              "exp_het_mean", "mratio_mean", "prop_low_afs_mean",
              "mean_allele_range")

# sumstats <- c("mean_allele_range", "num_alleles_mean", "mean_allele_size_sd")

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
  

# both models

for (i in c("bot", "neut")) {
  # mod <- "bot"
  mod <- i
  
  # subset sims_stats and sims_param for the given model
  stat_mod <- subset(sims_stats, subset = models == mod)
  par_mod <- subset(sims_param, subset = models == mod)
  
  # check whether a parameter can be estimated at all ----------------------------
  # (optional here as very time intense)
  
  # if ((calc_cv_for_abc == TRUE) & (run_parallel == FALSE)) {
  #   # before inference, we see whether a parameter can be estimated at all
  #   cv_nbot <- function(method, pars, nval = 5, tols = c(0.0005, 0.001, 0.005)){ #
  #     cv_res_rej <- cv4abc(data.frame(par_mod)[pars], stat_mod, nval = nval,
  #                          tols = tols, method = method)
  #   }
  #   #pars <- c("gsm_param")
  #   pars <- c("nbot", "nhist", "tbotend", "tbotstart", "mut_rate", "gsm_param")
  #   
  #   cv_res <- cv_nbot(method = method_cv, pars = pars, nval = nval_cv,  tols = tols_cv)
  #   out <- paste0("abc_estimates/cv_param_it1000_parallel", sims_name, ".RData")
  #   # write.table(cv_res, file = out, row.names = FALSE)
  #   save(cv_res, file = out)
  # }
  
  
  # in parallel ------------------------------------------------------------------
  if ((calc_cv_for_abc == TRUE) & (run_parallel == TRUE)) {
    # before inference, we see whether a parameter can be estimated at all
    
    cv_nbot <- function(iter, par_mod,stat_mod, method, pars, tols = c(0.0005)){ #
      cv_res_rej <- cv4abc(data.frame(par_mod)[pars], stat_mod, nval = 5,
                           tols = tols, method = method)
    }
    pars <- c("pop_size", "nbot", "nhist", "tbotend", "tbotstart", "mut_rate", "gsm_param")
    # pars <- c("nbot", "nhist", "mut_rate", "gsm_param")
    # tols_cv <- c(0.0005, 0.0001)
    
    cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - 10))
    clusterEvalQ(cl, c(library("abc")))
    all_cv_res <- parLapply(cl, 1:(nval_cv/5), cv_nbot, par_mod, stat_mod, method_cv, pars, tols_cv)
    stopCluster(cl)
    
    # first one
    all_cv <- all_cv_res[[1]]
    
    # extract values for binding
    # true values
    true_vals <- data.frame(do.call(rbind, lapply(all_cv_res, function(x) x$true)))
    cv_samples <- as.numeric(unlist(data.frame(do.call(rbind, lapply(all_cv_res, function(x) data.frame(x$cvsamples))))))
    # estimated values
    estim_vals <- list()
    for (i in 1:length(all_cv$estim)) {
      estim_vals[[i]] <- data.frame(do.call(rbind, lapply(all_cv_res, function(x) x$estim[[i]])))
    }
    names(estim_vals) <- names(all_cv$estim)
    # cv_samples
    all_cv$cvsamples <- cv_samples
    all_cv$true <- true_vals
    all_cv$estim <- estim_vals
    
    out <- paste0("model_evaluation/check4_params/cv_param_it1000_parallel_loclin_", sims_name,"_",mod,"_29", ".RData")
    # write.table(cv_res, file = out, row.names = FALSE)
    save(all_cv, file = out)
  }

}




if (abc_analysis == TRUE) {
  # run the actual abc analysis --------------------------------------------------
  for (i in c("bot", "neut")) {
    mod <- i
    
    ## load model probabilities
    model_probs <- read.table(paste0("results/model_probs/", 
                                     sims_name, "_model_selection.txt"))
    model_bot <- model_probs$bot > 0.5
    model_neut <- model_probs$bot <= 0.5
    # extract species names for species that have been bottlenecked according
    # to the model selection
    if (mod == "bot") {
      all_species <- row.names(all_sumstats_full)[model_bot]
    } else if (mod == "neut") {
      all_species <- row.names(all_sumstats_full)[model_neut]
    }
    
    
    # parameters to estimate posteriors
    all_parameters <- c("nbot", "pop_size", "mut_rate", "tbotstart", 
                        "tbotend", "nhist", "gsm_param") #
    
    # get all combinations of method, species and parameters
    all_args <- expand.grid(all_methods, all_species, all_parameters)
    all_args <- data.frame(apply(all_args, 2, as.character), stringsAsFactors = FALSE)
    names(all_args) <- c("methods", "species", "pars")
    
    ## run abc in parallel
    
    cl <- parallel::makeCluster(getOption("cl.cores", detectCores() - 10))
    clusterEvalQ(cl,library("abc"))
    # clusterExport(cl, varlist = c("all_sumstats", "par_mod", "stat_mod", "tol_abc"), envir = environment())
    
    one_abc <- function(arg_set, all_sumstats_full, par_mod, stat_mod, tol_abc){
      abc(target = all_sumstats_full[arg_set["species"], ], param = par_mod[, arg_set["pars"]], 
          sumstat = stat_mod, tol = tol_abc, method = arg_set["methods"])
    }
    abc_est <- parallel::parApply(cl, all_args, 1, one_abc,  all_sumstats_full, par_mod, stat_mod, tol_abc)
    
    stopCluster(cl)
    
    # create a list of 2. The first element ist the parameter data.frame for the abc. 
    # The second element are the corresponding abc objects.
    
    abc_full <- list(all_args, abc_est)
    save(abc_full, file = paste0("abc_estimates/abc_", sims_name,"_",mod,"_29", ".RData"))
  }
 
}
