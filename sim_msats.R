# simulating microsatellite data and calculating summary statistics
library(microsimr)
library(strataG)
library(splitstackshape)
library(stringr)
library(parallel)
library(data.table)



source("abc_tools.R")
# load all seal data
all_seals <- load_seal_data("../data/seal_data_largest_clust_and_pop.xlsx")


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




# which species?

i <- 9

# parameters for species
N_pop <- as.numeric(unlist(abundances[which(str_detect(abundances$species[i], names(seals))), "Abundance"]))
N_samp <- nrow(seals[[i]])
N_loc <- (ncol(seals[[i]]) - 3) / 2
##




run_sim <- function(niter, N_pop, N_samp, N_loc, model = c("bottleneck", "neutral", "decline")) {
  
  ## diploid pop size
  N0 <- round(runif(1, min = N_pop / 10, max = N_pop), 0)
 #  N0 <- N_pop
  # N0 <- 10000
  ## mutation rate
  # mu <- runif(1, min = 0.0005, max = 0.005)
  mu <- 0.0005
  
  ## theta
  theta <- 4 * N0 * mu
  
  # the 'start' here is actually the temporal end of the bottleneck. I´m just
  # thinking backwards in time here.
  
  #### bottleneck end ####
  # time is always thought in GENERATIONS
  latest_end <- 5 # generations ago
  earliest_end <- 15 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  end_bot <- runif(1, min = latest_end / (4*N0), max = earliest_end / (4*N0))
  
  
  #### bottleneck start ####
  latest_start <- 15 # generations ago
  earliest_start <- 50 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  start_bot <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
  
  ## bottleneck population size 1 - 1000, expressed relative to N0
  N_bot <- round(runif(1, min = 0.00001 * N0, max = 0.001 * N0), 4)
  
  ## historical populaiton size 1 - 100 times as big as current
  N_hist <- round(runif(1, min = N_pop / 10 , max = N_pop), 0)
  
  #
  N_hist_decl <- round(runif(1, min = N0 * 10 , max = N0 * 1000), 0)
  
  # exponential population decline
  # alpha = -(1/ end_bot) * log(N0/N_hist)
  
  if (model == "bottleneck") {
    ms_options <- paste("-eN", end_bot, N_bot, "-eN", start_bot, N_hist, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- NULL
  }
  
  if (model == "decline") {
    ms_options <- paste("-eN", start_bot, N_hist, sep = " ")
  }
  
  simd_data <- as.data.frame(microsimr::sim_microsats(theta = theta,
                                                      n_ind = N_samp,
                                                      n_loc = N_loc,
                                                      n_pop = 1,
                                                      ms_options = ms_options), stringsAsFactors = FALSE)
  
  # reshape a little bit
  simd_data <- simd_data[-c(1:2)]
  genotypes <- as.data.frame(splitstackshape::cSplit(simd_data, names(simd_data), "."))
  
  g_types_geno <- new("gtypes", genotypes, ploidy = 2)
  
  # summary statistics
  # according to DIYabc
  # number of alleles across loci
  num_alleles <- strataG::numAlleles(g_types_geno)
  num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
  num_alleles_sd <- sd(num_alleles, na.rm = TRUE)
  
  # allele size variance (actually, sd) across loci
  allele_size_sd <- unlist(lapply(seq(from = 1, to = ncol(simd_data)-1, by = 2),
                                  function(x) sd(unlist(simd_data[x:(x+1)]), na.rm = TRUE)))
  
  mean_allele_size_sd <- mean(allele_size_sd, na.rm = TRUE)
  sd_allele_size_sd <- sd(allele_size_sd, na.rm = TRUE)
  
  # expected heterozygosity
  exp_het <- exptdHet(g_types_geno)
  exp_het_mean <- mean(exp_het, na.rm = TRUE)
  exp_het_sd<- sd(exp_het, na.rm = TRUE)
  # observed heterozygosity
  obs_het <- obsvdHet(g_types_geno)
  obs_het_mean <- mean(obs_het, na.rm = TRUE)
  obs_het_sd <- sd(obs_het, na.rm = TRUE)
  
  # excess
  het_ratio <- mean(obs_het / exp_het, na.rm = TRUE)
  
  out <- data.frame(N0 = N0, t_bot = end_bot, t_hist = start_bot, n_bot = N_bot, n_hist = N_hist,
                    num_alleles_mean = num_alleles_mean, num_alleles_sd = num_alleles_sd,
                    mean_allele_size_sd = mean_allele_size_sd, sd_allele_size_sd = sd_allele_size_sd,
                    exp_het_mean = exp_het_mean,  exp_het_sd =  exp_het_sd,
                    obs_het_mean = obs_het_mean,  obs_het_sd =  obs_het_sd,
                    het_ratio = het_ratio)
  
  out
}



# number of simulations
num_sim <- 20000

# for (model in c("bottleneck", "neutral", "decline")) {
#   
#   # run on cluster
#   cl <- makeCluster(getOption("cl.cores", 2))
#   clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
#   #sims_bot <- do.call(rbind, parLapply(cl, 1:num_sim, run_sim, "bottleneck"))
#   # sims_bot <- do.call(rbind, parLapply(cl, 1:num_sim, run_sim, "bottleneck"))
#   sims <- parLapply(cl, 1:num_sim, run_sim, N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = model)
#   # system.time(sims_bot <- dplyr::rbind_all(sims))
#   sims_model <- as.data.frame(data.table::rbindlist(sims))
#   stopCluster(cl)
#   write.table(sims_model, file = paste("sims_", model, ".txt", sep = ""), row.names = FALSE)
#   
# }


cl <- makeCluster(getOption("cl.cores", 35))
clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck")
sims_bot <- as.data.frame(data.table::rbindlist(sims))
# sims_neut <- do.call(rbind, parLapply(1:num_sim, run_sim, "neutral"))
stopCluster(cl)

cl <- makeCluster(getOption("cl.cores", 35))
clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="neutral")
sims_neut <- as.data.frame(data.table::rbindlist(sims))
# sims_neut <- do.call(rbind, parLapply(1:num_sim, run_sim, "neutral"))
stopCluster(cl)
# boxplot(sims$exp_het)

cl <- makeCluster(getOption("cl.cores", 35))
clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))
sims <- parLapply(cl, 1:num_sim, run_sim,  N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model ="decline")
sims_decl <- as.data.frame(data.table::rbindlist(sims))
# sims_neut <- do.call(rbind, parLapply(1:num_sim, run_sim, "neutral"))
stopCluster(cl)
# boxplot(sims$exp_het)

# sims_bot vs. sims_neut
sims <- rbind(sims_bot, sims_neut, sims_decl)
sims$model <- c(rep("bot", num_sim), rep("neut", num_sim),  rep("decl", num_sim))
#
write.table(sims, file = "sims_3modhms.txt", row.names = FALSE)


