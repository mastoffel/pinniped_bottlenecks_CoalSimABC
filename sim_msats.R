# simulating microsatellite data and calculating summary statistics

library(devtools)
# install_github("mastoffel/sealABC", force = TRUE)
library(sealABC)
# devtools::install_github("andersgs/microsimr", build_vignettes = TRUE, force = TRUE)
library(microsimr)
library(strataG)
library(splitstackshape)
library(stringr)
library(parallel)
library(data.table)
library(readxl)
library(dplyr)
library(truncnorm)
# devtools::install_github('ericarcher/strataG', dependencies = TRUE)

# load all seal data ------------------------------
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")


### select 28 datasets //
# seals <- all_seals[c("antarctic_fur_seal", "galapagos_fur_seal", "stellers_sea_lion_cl_1",
#                      "grey_seal_orkneys", "harbour_seal_waddensee_cl_1", "galapagos_sea_lion",
#                      "south_american_fur_seal_cl_1", "hooded_seal", "mediterranean_monk_seal",
#                      "hawaiian_monk_seal", "bearded_seal_cl_1", "crabeater_seal",
#                      "leopard_seal", "arctic_ringed_seal", "ross_seal",
#                      "weddell_seal_cl_2", "northern_fur_seal_cl_1", "atlantic_walrus_cl_2",
#                      "nes_cl_2", "ses_cl_4", "california_sea_lion", "south_american_sea_lion",
#                      "new_zealand_sea_lion", "saimaa_ringed_seal_cl_1", "lagoda_ringed_seal",
#                      "baltic_ringed_seal", "new_zealand_fur_seal", "australian_fur_seal")]


# modify and simplify seal names
# names(seals) <- str_replace(names(seals), "_cl_1", "")
# names(seals) <- str_replace(names(seals), "_cl_2", "")

# # load abundances
# abundances <- read_excel("../data/abundances.xlsx", sheet = 1)
# 
# # load generation times
# gen_time <- read_excel("../data/seal_data_krueger.xlsx", sheet = 1)[c("dataset_name", "Generation time")]
# names(gen_time) <- str_replace_all(names(gen_time), " ", "_")
# gen_time <- gen_time %>% 
#   filter(!is.na(dataset_name)) %>% 
#   # generation time of arctic ringed seal for other ringed seals
#   mutate(Generation_time = ifelse(is.na(Generation_time), 6804, Generation_time)) %>%
#   mutate(gen_time = Generation_time / 356)
# # mean(gen_time$gen_time)



# which species? -----------------------------------------

# species_name <- "antarctic_fur_seal"

# parameters for species
# N_pop <- as.numeric(abundances[abundances$species == species_name, "Abundance"])

# N_pop <- 100000

# generation time in years
# gen_time <- as.numeric(gen_time[gen_time$dataset_name == species_name, "gen_time"])

# sampling 



gen_time <- 13.07 # mean generation time across species


run_sim <- function(niter, N_pop, model, gen_time) {
  
  N_samp <- round(runif(1, min = 20, max = 100), 0)
  N_loc <-  round(runif(1, min = 5, max = 30), 0)
  # N_loc <- 20
  
  if (niter%%5000 == 0) {
    if (!file.exists("num_iter/iterations.txt")){
      iter_df <- data.frame(model, niter)
      write.table(iter_df, file = "num_iter/iterations.txt", row.names = FALSE, col.names = TRUE)
    } else {
      if (file.exists("num_iter/iterations.txt")){
      all_iter <- read.table("num_iter/iterations.txt", header = 1)
      iter_df <- do.call(rbind, list(all_iter, data.frame(model, niter)))
      write.table(iter_df, file = "num_iter/iterations.txt", row.names = FALSE, col.names = TRUE)
    }
    }
  }
  
  ## diploid effective population size - size 1/10 to 1 of census
  # minium proportion of effective population size to census:
  # 1 /
  prop_prior_N <- 10
  
  N0 <- round(runif(1, min = N_pop / prop_prior_N, max = N_pop), 0)
  # to keep the historical population size in the same prior range as the current population
  # size, relative to N_pop
  Ne_prop <-   N_pop / N0  
  
  ## mutation rate
 # mu <- runif(1, min = 5e-05, max = 5e-03)
  mu <- 0.00005
  ## theta
  theta <- 4 * N0 * mu
  
  sequence_correct <- FALSE
  #### bottleneck end ####
  while (sequence_correct == FALSE) {
    # time is always thought in GENERATIONS
    latest_end_bot <- (2016 - 1950) / gen_time # generations ago, 2000 as latest end
    earliest_end_bot <- (2016 - 1800) / gen_time # generations ago, 1800 as earliest end
    # see ms manual. time is expressed in units of 4*N0
    end_bot <- runif(1, min = latest_end_bot / (4*N0), max = earliest_end_bot / (4*N0))
    
    #### bottleneck start ####
    latest_start_bot <- (2016 - 1850) / gen_time # generations ago, latest bottleneck start is 1850
    earliest_start_bot <- (2016 - 1700) / gen_time # generations ago, earliest bottleneck start is 1700
    # see ms manual. time is expressed in units of 4*N0
    start_bot <- runif(1, min = latest_start_bot / (4*N0), max = earliest_start_bot / (4*N0))
    
    # ensure that the end of the bottleneck is later than the start of the bottleneck
    if (end_bot < start_bot) sequence_correct <- TRUE
  }
  
  ## bottleneck population size reduction from 1/1000 to 2 individuals
  max_bot <- 1000 / N0 # 100 individuals
  min_bot <- 2 / N0   # 2 individuals
  N_bot <- round(runif(1, min = min_bot, max = max_bot), 7) 
  
  # hist(rlnorm(1000, meanlog = 4, sdlog = 3), breaks = 1000)
  # rtruncnorm(10000, a=0, b=1000, mean = 50, sd = 300)
  
  # N_bot <- rnorm(1, mean = 0.01, sd = 0.05)
  
  ## historical population size ranges in the same absolute priors (range) as the current population size
  N_hist_bot <- round(runif(1, min = (1 / prop_prior_N) * Ne_prop, max = Ne_prop), 7)
  
  # Ice age ended roughly 12000 years ago
  # earliest_start_decl_IA <- 12000/gen_time
  # latest start of decline roughly 8000 years ago
  # latest_start_decl_IA <- 8000/gen_time 
  
  # Truncated normal distribution including Ice Age End
  start_decl_IA <- (truncnorm::rtruncnorm(n = 1, a=5000, b=19000, mean = 12000, sd = 1500) / gen_time) / (4*N0)
  start_exp_IA <- (truncnorm::rtruncnorm(n = 1, a=5000, b=19000, mean = 12000, sd = 1500) / gen_time) / (4*N0)
  
  ## start of expansion and decline have same broad priors
 # start_decl_IA <- runif(1, min =  latest_start_decl_IA / (4*N0), max =  earliest_start_decl_IA / (4*N0))
 #  start_exp_IA <- runif(1, min =  latest_start_decl_IA / (4*N0), max =  earliest_start_decl_IA / (4*N0))
  
  ## recent decline and expansion similar to bottleneck
  start_decl_recent <- runif(1, min = latest_start_bot / (4*N0), max = earliest_start_bot / (4*N0))
  start_exp_recent <- runif(1, min = latest_start_bot / (4*N0), max = earliest_start_bot / (4*N0))
  
  # starts are now like the bottleneck start
  # start_decl <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
  # start_exp <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
  
  # historical pop sizes
  pop_param_decl_exp <- 10 # population was a maximum of 10 times larger or smaller in the past
  
  # N_hist_exp <- round(runif(1, min = 1/pop_param_decl_exp, max = 1), 7)
  # N_hist_decl <- round(runif(1, min = 1, max = pop_param_decl_exp), 7)
  
  # truncated normal distribution
  N_hist_decl <- truncnorm::rtruncnorm(n = 1,a = 1, b = 100, mean = pop_param_decl_exp, sd = 5)
  N_hist_exp <- truncnorm::rtruncnorm(n = 1,a = 0, b = 1, mean = 1/pop_param_decl_exp, sd = 0.05)
  
  # N_hist_exp <- rnorm(n = 1, mean = 1 / pop_param_decl_exp, sd = 0.05)
  
  # exponential population decline or increase
  # care: time is now start of bottleneck
  alpha_decl_IA <-  -(1/ start_decl_IA) * log((N_hist_decl * N0) / N0)
  alpha_exp_IA <-   -(1/ start_exp_IA) * log((N_hist_exp * N0) / N0 )
  
  alpha_decl_recent <-  -(1/  start_decl_recent) * log((N_hist_decl * N0) / N0)
  alpha_exp_recent <-   -(1/ start_exp_recent) * log((N_hist_exp * N0) / N0 )
  # 
  
  if (model == "bottleneck") {
    ms_options <- paste("-eN", end_bot, N_bot, "-eN", start_bot, N_hist_bot, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- paste("-eN", start_bot, N_hist_bot, sep = " ")
  }
  
  if (model == "decline_recent") {
    ms_options <- paste("-G", alpha_decl_recent, "-eN", start_decl_recent, N_hist_decl, sep = " ") 
  }
  
  if (model == "expansion_recent") {
    ms_options <- paste("-G", alpha_exp_recent, "-eN", start_exp_recent, N_hist_exp, sep = " ") 
  }
  
  if (model == "decline_IA") {
   # ms_options <- paste("-G", alpha_decl_IA, "-eN", start_decl_IA, N_hist_decl, sep = " ") 
    ms_options <- paste("-eN", start_decl_IA, N_hist_decl, sep = " ") 
  }
  
  if (model == "expansion_IA") {
   # ms_options <- paste("-G", alpha_exp_IA, "-eN", start_exp_IA, N_hist_exp, sep = " ") 
    ms_options <- paste("-eN", start_exp_IA, N_hist_exp, sep = " ") 
  }
  
  p_single = runif(1, min = 0.7, max = 0.99) # probability of multi-step mutation is 0.2
  sigma2_g = runif(1, min = 30, max = 70) # typical step-size ~7
  
  simd_data <- as.data.frame(microsimr::sim_microsats(theta = theta,
                                                      n_ind = N_samp,
                                                      n_loc = N_loc,
                                                      n_pop = 1,
                                                      p_single = p_single,
                                                      sigma2 = sigma2_g,
                                                      ms_options = ms_options), stringsAsFactors = FALSE)
  
  sum_stats <- sealABC::mssumstats(simd_data)
  
  sim_param <- data.frame(N_pop, N0, mu, theta,
                          start_bot, end_bot,
                          N_bot, N_hist_bot,
                          start_decl_IA, start_decl_recent,
                          start_exp_IA, start_exp_recent,
                          N_hist_decl, N_hist_exp,
                          p_single, sigma2_g,
                          N_loc, N_samp)
  out <- cbind(sim_param, sum_stats)
  
  out
}


### number of all simulations
num_sim <- 50000

######### all simulations for 100000 ###############
# census size
N_pop <- 100000


all_models = c("bottleneck", "neutral", "decline_recent", "expansion_recent", "decline_IA", "expansion_IA")

run_sim_per_mod <- function(model){
  
  cl <- makeCluster(getOption("cl.cores", detectCores() - 2))
  clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
  
  sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, model = model, gen_time = gen_time)
  sims_df <- as.data.frame(data.table::rbindlist(sims))
  
  stopCluster(cl)
  
  sims_df
  
}
# 
# run_sim_per_mod <- function(model){
#   sims <- lapply(1:num_sim, run_sim,  N_pop = N_pop, model = model, gen_time = gen_time)
#   sims_df <- as.data.frame(data.table::rbindlist(sims))
#   sims_df
# }


sims <- do.call(rbind, lapply(all_models, run_sim_per_mod))
sims$model <- c(rep("bot", num_sim), rep("neut", num_sim), rep("decl_IA", num_sim), 
                rep("exp_IA", num_sim), rep("decl_rec", num_sim), rep("exp_rec", num_sim))

write.table(sims, file = "sims_50000_6mod_lowmu_highsigma.txt", row.names = FALSE)

# sims <- rbind(sims_bot, sims_neut, sims_decl_IA, sims_exp_IA, sims_decl_rec, sims_exp_rec) #sims_exp








# 
# ######### all simulations for 10000 ###############
# # census size
# N_pop <- 10000
# 
# cl <- makeCluster(getOption("cl.cores", detectCores() - 1))
# clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck", gen_time = gen_time)
# sims_bot <- as.data.frame(data.table::rbindlist(sims))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="neutral", gen_time = gen_time)
# sims_neut <- as.data.frame(data.table::rbindlist(sims))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="decline", gen_time = gen_time)
# sims_decl <- as.data.frame(data.table::rbindlist(sims))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="expansion", gen_time = gen_time)
# sims_exp <- as.data.frame(data.table::rbindlist(sims))
# 
# stopCluster(cl)
# 
# sims <- rbind(sims_bot, sims_neut, sims_decl, sims_exp) #sims_exp
# sims$model <- c(rep("bot", num_sim), rep("neut", num_sim), rep("decl", num_sim), rep("exp", num_sim))
# 
# write.table(sims, file = "sims_smallN_10000.txt", row.names = FALSE)
# 
# 
# 
# 
# ######### all simulations for 1000000 ###############
# # census size
# N_pop <- 1000000
# 
# cl <- makeCluster(getOption("cl.cores", detectCores() - 1))
# clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck", gen_time = gen_time)
# sims_bot <- as.data.frame(data.table::rbindlist(sims))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="neutral", gen_time = gen_time)
# sims_neut <- as.data.frame(data.table::rbindlist(sims))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="decline", gen_time = gen_time)
# sims_decl <- as.data.frame(data.table::rbindlist(sims))
# 
# sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="expansion", gen_time = gen_time)
# sims_exp <- as.data.frame(data.table::rbindlist(sims))
# 
# stopCluster(cl)
# 
# sims <- rbind(sims_bot, sims_neut, sims_decl, sims_exp) #sims_exp
# sims$model <- c(rep("bot", num_sim), rep("neut", num_sim), rep("decl", num_sim), rep("exp", num_sim))
# 
# write.table(sims, file = "sims_largeN_1000000.txt", row.names = FALSE)
