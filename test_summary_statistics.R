library(diveRsity)
library(sealABC)
library(strataG)
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
# load all seal data
all_seals <- read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

genotypes <- all_seals[[1]]


##### with strataG

g_types_geno <- new("gtypes", genotypes[4:ncol(genotypes)], ploidy = 2)

# number of alleles across loci
num_alleles <- strataG::numAlleles(g_types_geno)
num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
num_alleles_sd <- sd(num_alleles, na.rm = TRUE)

# allele size variance (actually, sd) across loci
allele_size_sd <- unlist(lapply(seq(from = 1, to = ncol(genotypes)-1, by = 2),
                                function(x) sd(unlist(genotypes[x:(x+1)]), na.rm = TRUE)))

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

out <- data.frame(num_alleles_mean = num_alleles_mean, num_alleles_sd = num_alleles_sd,
                  mean_allele_size_sd = mean_allele_size_sd, sd_allele_size_sd = sd_allele_size_sd,
                  exp_het_mean = exp_het_mean,  exp_het_sd =  exp_het_sd,
                  obs_het_mean = obs_het_mean,  obs_het_sd =  obs_het_sd,
                  het_ratio = het_ratio)



# with diveRsity

?basicStats
data("Test_data")

out <- make_genepop(genotypes)
out[] <- lapply(out, as.character)
out <- data.frame(apply(out, 2, function(x) {
  x[is.na(x)] <- ""
  x
}), stringsAsFactors = FALSE)
out[1,1] <- "BLANK"
out$V1 <- as.factor(out$V1)
basicStats(out, rarefaction = FALSE)
