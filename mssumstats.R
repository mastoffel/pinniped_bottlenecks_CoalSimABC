#' takes output from microsimr and outputs summary statistics
#' 
#' 


mssumstats <- function(simd_data) {
  
  # reshape a little bit
  simd_data <- simd_data[-c(1:2)]
  # genotypes <- splitstackshape::cSplit_f(simd_data, splitCols = names(simd_data), sep = "|")
  genotypes <- splitstackshape::cSplit(simd_data, splitCols = names(simd_data), sep = ".")
  
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

