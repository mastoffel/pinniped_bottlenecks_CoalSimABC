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

gen_time <- 13.07 # mean generation time across species


run_sim <- function(niter, model, gen_time) {
  
  # N_samp <- round(runif(1, min = 20, max = 500), 0)
  #  N_loc <-  round(runif(1, min = 5, max = 30), 0)
  
  N_samp <-  50
  N_loc <- 10
  
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
  
  ## diploid effective population size prior: from 1 to 500000
  # N0 <- round(runif(1, min = 1, max = 500000))
  N0 <- rtruncnorm(1, a=1, b=300000, mean = 10000, sd = 50000)
  
  ## mutation rate
  # original mutation prior (however, posteriors rather around 0.001)
 # mu <- runif(1, min = 0.00001, max = 0.009)
  # prior similar to Pritchard et al 1999
  mu <- rgamma(1, 2, rate = 10000e-1)
  # mu <- rtruncnorm(1, a=0, b=0.001, mean = 0.0001, sd = 0.0002)
  
  ## theta
  theta <- 4 * N0 * mu
  
  sequence_correct <- FALSE
  #### bottleneck end ####
  while (sequence_correct == FALSE) {
    # time is always thought in GENERATIONS
    latest_end_bot <- (2016 - 2000) / gen_time # generations ago, 2000 as latest end
    earliest_end_bot <- (2016 - 1800) / gen_time # generations ago, 1800 as earliest end
    # see ms manual. time is expressed in units of 4*N0
    end_bot <- runif(1, min = latest_end_bot / (4*N0), max = earliest_end_bot / (4*N0))
    
    #### bottleneck start ####
    latest_start_bot <- (2016 - 1850) / gen_time # generations ago, latest bottleneck start is 1850
    earliest_start_bot <- (2016 - 1650) / gen_time # generations ago, earliest bottleneck start is 1700
    # see ms manual. time is expressed in units of 4*N0
    start_bot <- runif(1, min = latest_start_bot / (4*N0), max = earliest_start_bot / (4*N0))
    
    # ensure that the end of the bottleneck is later than the start of the bottleneck
    if (end_bot < start_bot) sequence_correct <- TRUE
  }
  
  ### uniform prior for bottleneck
  # N_bot <- round(rtruncnorm(1, a=1, b=1000, mean = 10, sd = 200), 0)
  N_bot <- N0
  while (!(N_bot < N0)){
    N_bot <- round(runif(1, min = 1, max = 500),0)
  }
  N_bot <- N_bot / N0
  
  ## historical population size ranges in the same absolute priors (range) as the current population size
  # At the moment prop_prior_N is N_pop, which means that the historical population size
  # ranges from 1 to N_pop, i.e. has the same range than the current pop size N0
  #N_hist_bot <- round(runif(1, min = 1, max = 500000)) / N0
  N_hist_bot <- rtruncnorm(1, a=1, b=300000, mean = 10000, sd = 50000)/ N0
  # 
  if (model == "bottleneck") {
    ms_options <- paste("-eN", end_bot, N_bot, "-eN", start_bot, N_hist_bot, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- paste("-eN", start_bot, N_hist_bot, sep = " ")
    # ms_options <- NULL
  }
  
  # p_single <- rtruncnorm(1, a=0.7, b=1, mean = 0.85, sd = 0.07)
  # tryout
  p_single <-  runif(1, min = 0.7, max = 1) # 1 - probability of multi-step mutations
  
  
  # sigma2_g <- runif(1, min = 1, max = 30) 
  # sigma2_g <- runif(1, min = 1, max = 15)
  # new try with a gamma distribution
  sigma2_g <- rgamma(1, 4, 0.9)
  
  simd_data <- as.data.frame(microsimr::sim_microsats(theta = theta,
                                                      n_ind = N_samp,
                                                      n_loc = N_loc,
                                                      n_pop = 1,
                                                      p_single = p_single,
                                                      sigma2 = sigma2_g,
                                                      mutation_model = "tpm",
                                                      ms_options = ms_options), stringsAsFactors = FALSE)
  
  sum_stats <- sealABC::mssumstats(simd_data, datatype = "microsimr")
  
  sim_param <- data.frame(N0, mu, theta,
                          start_bot, end_bot,
                          N_bot, N_hist_bot,
                          p_single, sigma2_g,
                          N_loc, N_samp)
  out <- cbind(sim_param, sum_stats)
  
  out
}


### number of all simulations
num_sim <- 100000
file_ext <- "1mio_sims_gammaprior.txt"

all_models = c("bottleneck", "neutral")

######### all simulations for 1000000 ###############
run_sim_per_mod <- function(model){

  cl <- makeCluster(getOption("cl.cores", detectCores()-3))
  clusterEvalQ(cl, c(library("strataG"), library("splitstackshape"), library("truncnorm")))

  sims <- parLapply(cl, 1:num_sim, run_sim, model = model, gen_time = gen_time)
  sims_df <- as.data.frame(data.table::rbindlist(sims))

  stopCluster(cl)

  sims_df

}

# run_sim_per_mod <- function(model){
#   sims <-lapply(1:num_sim, run_sim, model = model, gen_time = gen_time)
#   sims_df <- as.data.frame(data.table::rbindlist(sims))
#   sims_df
# }

sims <- do.call(rbind, lapply(all_models, run_sim_per_mod))
sims$model <- c(rep("bot", num_sim), rep("neut", num_sim))

write.table(sims, file = paste0("onepopprior_", file_ext), row.names = FALSE)

