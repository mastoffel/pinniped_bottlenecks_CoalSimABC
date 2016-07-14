# simulating microsatellite data and calculating summary statistics

library(devtools)
# install_github("mastoffel/sealABC", dependecies = TRUE)
library(sealABC)
# devtools::install_github("andersgs/microsimr")
library(microsimr)
library(strataG)
library(splitstackshape)
library(stringr)
library(parallel)
library(data.table)
library(readxl)
library(dplyr)

# devtools::install_github('ericarcher/strataG', dependencies = TRUE)

# load all seal data ------------------------------
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")


### select 28 datasets //
seals <- all_seals[c("antarctic_fur_seal", "galapagos_fur_seal", "stellers_sea_lion_cl_2",
                     "grey_seal_orkneys", "harbour_seal_waddensee_cl_2", "galapagos_sea_lion",
                     "south_american_fur_seal_cl_2", "hooded_seal", "mediterranean_monk_seal",
                     "hawaiian_monk_seal", "bearded_seal_cl_2", "crabeater_seal",
                     "leopard_seal", "arctic_ringed_seal", "ross_seal",
                     "weddell_seal_cl_1", "northern_fur_seal_cl_1", "atlantic_walrus_cl_1",
                     "nes_cl_1", "ses_cl_1", "california_sea_lion", "south_american_sea_lion",
                     "new_zealand_sea_lion", "saimaa_ringed_seal_cl_2", "lagoda_ringed_seal",
                     "baltic_ringed_seal", "new_zealand_fur_seal", "australian_fur_seal")]

# modify and simplify seal names
names(seals) <- str_replace(names(seals), "_cl_1", "")
names(seals) <- str_replace(names(seals), "_cl_2", "")

# load abundances
abundances <- read_excel("../data/abundances.xlsx", sheet = 1)

# load generation times
gen_time <- read_excel("../data/seal_data_krueger.xlsx", sheet = 1)[-1]
names(gen_time) <- str_replace(names(gen_time), " ", "_")
gen_time <- gen_time %>% 
  select(dataset_name, Generation_time) %>%
  filter(!is.na(dataset_name)) %>% 
  # generation time of arctic ringed seal for other ringed seals
  mutate(Generation_time = ifelse(is.na(Generation_time), 6804, Generation_time)) %>%
  mutate(gen_time = Generation_time / 356)



# which species? -----------------------------------------

species_name <- "antarctic_fur_seal"

# parameters for species
N_pop <- as.numeric(abundances[abundances$species == species_name, "Abundance"])

# generation time in years
gen_time <- as.numeric(gen_time[gen_time$dataset_name == species_name, "gen_time"])

# sampling 
N_samp <- 150
N_loc <- 50
##



run_sim <- function(niter, N_pop, N_samp, N_loc, model = c("bottleneck", "neutral"), gen_time) {
  
  ## diploid effective population size - size 1/10 to 1 of census
  # minium proportion of effective population size to census:
  # 1 /
  prop_prior_N <- 10
  
  N0 <- round(runif(1, min = N_pop / prop_prior_N, max = N_pop), 2)
  # to keep the historical population size in the same prior range as the current population
  # size, relative to N_pop
  Ne_prop <-   N_pop / N0
  
  ## mutation rate
  mu <- runif(1, min = 0.00005, max = 0.0005)
  ## theta
  theta <- 4 * N0 * mu
  
  sequence_correct <- FALSE

  while (sequence_correct == FALSE) {
    # time is always thought in GENERATIONS
    latest_end <- (2016 - 2000) / gen_time # generations ago, 2000 as latest end
    earliest_end <- (2016 - 1800) / gen_time # generations ago, 1800 as earliest end
    # see ms manual. time is expressed in units of 4*N0
    end_bot <- runif(1, min = latest_end / (4*N0), max = earliest_end / (4*N0))
    
    #### bottleneck start ####
    latest_start <- (2016 - 1900) / gen_time # generations ago, latest bottleneck start is 1900
    earliest_start <- (2016 - 1700) / gen_time # generations ago, earliest bottleneck start is 1700
    # see ms manual. time is expressed in units of 4*N0
    start_bot <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
    
    # ensure that the end of the bottleneck is later than the start of the bottleneck
    if (end_bot < start_bot) sequence_correct <- TRUE
  }
  
  ## bottleneck population size reduction from 1/1000 to 1/1000000
  N_bot <- round(runif(1, min = 0.000001, max = 0.001), 7) 
  
  ## historical population size ranges in the same absolute priors as the current population size
  N_hist <- round(runif(1, min = (1 / prop_prior_N) * Ne_prop, max = Ne_prop), 5)
  #
  
  # exponential population decline
  # alpha = -(1/ end_bot) * log(N0/N_hist)
  
  if (model == "bottleneck") {
    ms_options <- paste("-eN", end_bot, N_bot, "-eN", start_bot, N_hist, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- paste("-eN", start_bot, N_hist, sep = " ")
  }
  
  p_single = runif(1, min = 0.6, max = 0.95) # probability of multi-step mutation is 0.2
  sigma2_g = runif(1, min = 5, max = 70) # typical step-size ~7
  
  simd_data <- as.data.frame(microsimr::sim_microsats(theta = theta,
                                                      n_ind = N_samp,
                                                      n_loc = N_loc,
                                                      n_pop = 1,
                                                      p_single = p_single,
                                                      sigma2 = sigma2_g,
                                                      ms_options = ms_options), stringsAsFactors = FALSE)
  
  sum_stats <- sealABC::mssumstats(simd_data)
  
  sim_param <- data.frame(N0, mu, theta,
                          start_bot, end_bot,
                          N_bot, N_hist,
                          p_single, sigma2_g)
  out <- cbind(sim_param, sum_stats)
  
  out
}


# lineprof(run_sim( N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck"))



# number of simulations
num_sim <- 100000


cl <- makeCluster(getOption("cl.cores", detectCores()))
clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))

sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck", gen_time = gen_time)
sims_bot <- as.data.frame(data.table::rbindlist(sims))

sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="neutral", gen_time = gen_time)
sims_neut <- as.data.frame(data.table::rbindlist(sims))

stopCluster(cl)




sims <- rbind(sims_bot, sims_neut) 
sims$model <- c(rep("bot", num_sim), rep("neut", num_sim))


write.table(sims, file = "sims_2mod_afs.txt", row.names = FALSE)

