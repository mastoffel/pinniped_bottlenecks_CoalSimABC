# simulating microsatellite data and calculating summary statistics

library(microsimr)
library(strataG)
library(splitstackshape)


  
  ## diploid pop size
  N0 <- round(runif(1, min = 1000, max = 300000), 0)
  # N0 <- 100000
  
  ## mutation rate
  mu <- runif(1, min = 0.0005, max = 0.005)
  # mu <- 0.0005
  
  ## theta
  theta <- 4 * N0 * mu
  
  # the 'start' here is actually the temporal end of the bottleneck. I´m just 
  # thinking backwards in time here.
  
  
  #### bottleneck end ####
  # time is always thought in GENERATIONS
  latest_end <- 1 # generations ago
  earliest_end <- 30 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  end_bot <- runif(1, min = latest_end / (4*N0), max = earliest_end / (4*N0))
  
  
  #### bottleneck start ####
  latest_start <- 15 # generations ago
  earliest_start <- 100 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  start_bot <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
  
  ## bottleneck population size 1 - 10000
  N_bot <- round(runif(1, min = 0.00001, max = 0.001), 4)
  N_hist <- round(runif(1, min = 1, max = 100), 0)
  
  # exponential population decline
  
  alpha = -(1/ end_bot) * log(N0/N_hist)
  
  if (model == "bottleneck") {
    ms_options <- paste("-eN", end_bot, N_bot, "-eN", start_bot, N_hist, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- NULL
  }
  
  if (model == "exp_decline") {
    ms_options <- paste("-G", alpha,"-eG", start_bot, "0", "-eN", start_bot, N_hist, sep = " ")
  }
  
  simd_data <- as.data.frame(microsimr::sim_microsats(theta = theta,
                                                      n_ind = 100,
                                                      n_loc = 20,
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
  allele_size_sd <- unlist(lapply(seq(from = 1, to = ncol(simd_data), by = 2), 
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
  
  out <- data.frame(N0 = N0, t_bot = end_bot, t_hist = start_bot, n_bot = N_bot, n_hist = N_hist, 
                    num_alleles_mean = num_alleles_mean, num_alleles_sd = num_alleles_sd,
                    mean_allele_size_sd = mean_allele_size_sd, sd_allele_size_sd = sd_allele_size_sd,
                    exp_het_mean = exp_het_mean,  exp_het_sd =  exp_het_sd, 
                    obs_het_mean = obs_het_mean,  obs_het_sd =  obs_het_sd)
  
  out 


