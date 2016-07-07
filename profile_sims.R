# simulating microsatellite data and calculating summary statistics
library(microsimr)
library(strataG)
library(splitstackshape)
library(stringr)
library(parallel)
library(data.table)



run_sim <- function( N_pop, N_samp,  N_loc, model = c("bottleneck", "neutral", "decline")) {
  
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
  latest_start <- earliest_end + 1 # generations ago / ensure that start is before end
  earliest_start <- 50 # generations ago
  # see ms manual. time is expressed in units of 4*N0
  start_bot <- runif(1, min = latest_start / (4*N0), max = earliest_start / (4*N0))
  
  ## bottleneck population size 1 - 1000, expressed relative to N0
  N_bot <- round(runif(1, min = 0.000001, max = 0.0001), 4)
  
  ## historical populaiton size 1 - 100 times as big as current
  N_hist <- round(runif(1, min = N_pop / 10, max = N_pop), 0)
  
  #
  N_hist_decl <- round(runif(1, min = N0, max = N0 * 1000), 0)
  
  # exponential population decline
  # alpha = -(1/ end_bot) * log(N0/N_hist)
  
  if (model == "bottleneck") {
    ms_options <- paste("-eN", end_bot, N_bot, "-eN", start_bot, N_hist, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- paste("-eN", start_bot, N_hist, sep = " ")
  }
  
  if (model == "decline") {
    ms_options <- paste("-eN", start_bot, N_hist_decl, sep = " ")
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
num_sim <- 1
N_samp <- 20
N_loc <- 15



measure_time <- function(N_pop) {
  
  out1 <- system.time(run_sim(N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck"))
  out2 <- system.time(run_sim(N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "neutral"))
  out3 <- system.time(run_sim(N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "decline"))
   
  out <- out1[3] + out2[3] + out3[3]
}
## more than 100 000 gets unrealistic
out <- lapply(seq(from = 10, to = 10000, by = 1000), measure_time)
plot(unlist(out))




library(data.table)
measure_ss <- function(N_pop, model = "bottleneck") {
  
  rep_sim <- function(niter) {
    out <- run_sim(N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = model)
  }
  out <- lapply(1:30, rep_sim)
  sims <- as.data.frame(data.table::rbindlist(out))
}

out <- lapply(seq(from = 100, to = 500000, by = 50000), measure_ss)

all_sims <- as.data.frame(data.table::rbindlist(out))

library(ggplot2)
ggplot(all_sims, aes(x = as.factor(N0), y = exp_het_mean)) +
  geom_boxplot() 
  

plot(unlist(lapply(out, function(x) x$exp_het_mean)))



