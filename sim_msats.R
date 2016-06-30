# simulating microsatellite data and calculating summary statistics

library(microsimr)
library(strataG)
library(splitstackshape)


run_sim <- function(niter, model = c("bottleneck", "neutral")) {
  
  # diploid pop size
  N0 <- 100000
  # mutation rate
  mu <- 0.0005
  # theta
  theta <- 4 * N0 * mu
  
  start_bot <- round(runif(1, min = 1, max = 30), 0)
  end_bot <- round(runif(1, min = 15, max = 100), 0)
  
  N_bot <- round(runif(1, min = 0.0001, max = 0.001), 4)
  N_hist <- round(runif(1, min = 1, max = 100), 0)
  
  if (model == "bottleneck") {
    ms_options <- paste("-eN", start_bot, N_bot, "-eN", end_bot, N_hist, sep = " ")
  }
  
  if (model == "neutral") {
    ms_options <- NULL
  }
  
  simd_data <- as.data.frame(microsimr::sim_microsats(theta = theta,
                                                      n_ind = 50,
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
  
  out <- data.frame(t_bot = start_bot, t_hist = end_bot, n_bot = N_bot, n_hist = N_hist, 
                    num_alleles_mean = num_alleles_mean, num_alleles_sd = num_alleles_sd,
                    mean_allele_size_sd = mean_allele_size_sd, sd_allele_size_sd = sd_allele_size_sd,
                    exp_het_mean = exp_het_mean,  exp_het_sd =  exp_het_sd, 
                    obs_het_mean = obs_het_mean,  obs_het_sd =  obs_het_sd)
  
  out 
}

library(parallel)

# number of simulations
num_sim <- 1000

# run on cluster
cl <- makeCluster(getOption("cl.cores", 20))
clusterEvalQ(cl, c(library("strataG"), library("splitstackshape")))

sims_bot <- do.call(rbind, parLapply(cl, 1:num_sim, run_sim, "bottleneck"))
sims_neut <- do.call(rbind, lapply(1:num_sim, run_sim, "neutral"))
stopCluster(cl)


# boxplot(sims$exp_het)

# sims_bot vs. sims_neut
sims <- rbind(sims_bot, sims_neut)
sims$model <- c(rep("bot", num_sim), rep("neut", num_sim))

write.table(sims, file = "sims.txt", row.names = FALSE)

# 
# par(mfcol = c(2,2))
# boxplot(sims$num_alleles_mean ~ sims$model)
# boxplot(sims$mean_allele_size_sd ~ sims$model)
# boxplot(sims$exp_het_mean ~ sims$model)
# boxplot(sims$obs_het_mean ~ sims$model)
