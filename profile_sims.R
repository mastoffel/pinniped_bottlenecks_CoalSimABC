# simulating microsatellite data and calculating summary statistics

library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(microsimr)
library(strataG)
library(splitstackshape)
library(stringr)
library(parallel)
library(data.table)
library(readxl)


library(lineprof)



# which species?

i <- 1

# parameters for species
# N_pop <- as.numeric(unlist(abundances[which(str_detect(abundances$species[i], names(seals))), "Abundance"]))

  

time_sim_microsats <- function(N_pop){
  
  N_pop <- 100
  N_samp <- 100
  N_loc <- 20
  ##
  N0 <- N_pop
  mu <- 0.0005
  theta <- 4 * N0 * mu
  lineprof(simd_data <- microsimr::sim_microsats(theta = theta, n_ind = N_samp, n_loc = N_loc,  n_pop = 1))
  
}

                    
                                                   
                                            
  
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


lineprof(run_sim( N_pop = N_pop, N_samp = N_samp, N_loc = N_loc, model = "bottleneck"))



# number of simulations
num_sim <- 1000

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
write.table(sims, file = "sims_3modafs.txt", row.names = FALSE)





