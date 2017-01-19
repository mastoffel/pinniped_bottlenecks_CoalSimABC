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

seal_descriptives <- read_excel("../data/all_data_seals.xlsx")
# add factor for abundance --> effective population size considered to be a maximum of one fifth of the abundance
seal_descriptives %<>% mutate(abund_level = ifelse(Abundance < 20000, "5k", ifelse(Abundance < 300000, "50k", "500k")))

###### prepare data ########

# load all_seals data for the 28 full datasets
all_seals_full <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")[1:28]

# calculate summary statistics
all_sumstats_full <- lapply(all_seals_full, function(x) mssumstats(x[4:ncol(x)], type = "microsats", data_type = "empirical"))
# as data.frame
all_sumstats_full <- do.call(rbind, all_sumstats_full)

# which ss to use
# names(sims)
sumstats <- c("num_alleles_mean", "prop_low_afs_mean",   
              "mean_allele_range",  "mean_allele_size_var",
              "exp_het_mean")

all_sumstats_full <- all_sumstats_full[sumstats]

# run loop for all simulations
for (pop_size in c("5k", "50k", "500k")){
  
  # load simulations
  path_to_sims <- paste0("sims_pop", pop_size, "_sim300k_restr.txt")
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
  
  # some processing
  # extract names of all models
  all_models <- names(table(models))
  
  # divide stats and parameters
  sims_stats <- sims[sumstats] 
  sims_param <- sims[params]
  
  # run the actual abc analysis
  
  #for (mod in c("bot")){
  
  # just use bottleneck model
  mod <- "bot"
  
    stat_mod <- subset(sims_stats, subset=models==mod)
    par_mod <- subset(sims_param, subset=models==mod)
    
    ## accuracy of parameter estimate
    # if (mod == "bot"){
    
    # # before inference, we see whether a parameter can be estimated at all
    # cv_nbot <- function(method, nval = 50, tols = c(.001,.01, 0.005)){
    #   cv_res_rej <- cv4abc(data.frame(Na=par_mod[,"N_bot"]), stat_mod, nval=nval,
    #                        tols=tols, method=method)
    # }
    
   # all_cv_nbot <- lapply(c("loclinear", "ridge"), cv_nbot, nval = 30)
   # assign(paste0("cv_nbot_", pop_size), all_cv_nbot)
   # 

    par_mod <- par_mod %>% 
      mutate(N_bot = N0 * N_bot) %>%
      mutate(N_hist_bot = N0 * N_hist_bot) %>%
      mutate(start_bot = 4 * N0 * start_bot) %>%
      mutate(end_bot = 4 * N0 * end_bot) 
    
    # abc method
    all_methods <- c("ridge", "loclinear", "neuralnet") # 
    all_methods <- c("neuralnet")
    # species names
    all_species <- row.names(all_sumstats)
    # parameters to estimate posteriors
    all_parameters <- c("N_bot", "N0", "mu", "start_bot", "end_bot", "N_hist_bot", "sigma2_g", "p_single") # "N0", "mu", "start_bot", "end_bot", "N_hist_bot", "sigma2_g"
    
    # get all combinations of method, species and parameters
    all_args <- expand.grid(all_methods,all_species, all_parameters)
    all_args <- data.frame(apply(all_args, 2, as.character), stringsAsFactors = FALSE)
    names(all_args) <- c("methods", "species", "pars")
    
    ## run abc
    abc_est <- apply(all_args, 1, function(x) {
      abc(target = all_sumstats[x["species"], ], param = par_mod[, x["pars"]], 
          sumstat = stat_mod, tol = 0.001, method=x["methods"])
    })
    
    # create a list of 2. The first element ist the parameter data.frame for the abc. The second element are the corresponding abc objects.
    # assign(paste0("abc_", pop_size), list(all_args, abc_est))
    # saved_file_name <- paste(paste0("abc_", pop_size),"RData",sep=".")
    # save(list(all_args, abc_est), saved_file_name)
    
    if (pop_size == "5k"){
      abc_5k <- list(all_args, abc_est)
      save(abc_5k, file = "abc_estimates/abc_5k.RData")
    }
    if (pop_size == "50k"){
      abc_50k <- list(all_args, abc_est)
      save(abc_50k, file = "abc_estimates/abc_50k.RData")
    }
    if (pop_size == "500k"){
      abc_500k <- list(all_args, abc_est)
      save(abc_500k, file = "abc_estimates/abc_500k.RData")
    }
}  


  
  
  
  