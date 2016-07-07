# analyse simulated microsatellites under three different models000


sims <- read.table("sims_3modhms.txt", stringsAsFactors = FALSE)
names(sims) <- unlist(sims[1, ])
sims <- sims[-1, ]
# transform everything to numeric
sims <- cbind(do.call(cbind, lapply(sims[-ncol(sims)], as.numeric)), sims[ncol(sims)])

#
par(mfcol = c(2,3))
boxplot(sims$num_alleles_mean ~ sims$model)
boxplot(sims$mean_allele_size_sd ~ sims$model)
boxplot(sims$exp_het_mean ~ sims$model)
boxplot(sims$obs_het_mean ~ sims$model)
boxplot(sims$het_ratio ~ sims$model)

par(mfcol = c(2,2))
boxplot(sims$num_alleles_sd ~ sims$model)
boxplot(sims$sd_allele_size_sd ~ sims$model)
boxplot(sims$exp_het_sd ~ sims$model)
boxplot(sims$obs_het_sd ~ sims$model)



# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

# hawaiian monk seal
genotypes <- all_seals[[1]]
genotypes <- genotypes[4:ncol(genotypes)]





#### calculate empirical summary statistics
g_types_geno <- new("gtypes", genotypes, ploidy = 2)
# summary statistics
# according to DIYabc
# number of alleles across loci
num_alleles <- strataG::numAlleles(g_types_geno)
num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
num_alleles_sd <- sd(num_alleles, na.rm = TRUE)

# allele size variance (actually, sd) across loci
allele_size_sd <- unlist(lapply(seq(from = 1, to = ncol(genotypes), by = 2), 
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


obs_stats <- data.frame(num_alleles_mean, num_alleles_sd, mean_allele_size_sd, sd_allele_size_sd,
                        exp_het_mean, exp_het_sd, obs_het_mean, obs_het_sd, het_ratio)

# divide stats and parameters
sims_stats <- sims[6:(ncol(sims)-1)] # last column is model factor
sims_param <- sims[1:5]
models <- sims[["model"]]

library(abc)

# abc
cv.modsel <- cv4postpr(models, sims_stats, nval=5, tol=.5, method="mnlogistic")
s <- summary(cv.modsel)
plot(cv.modsel, names.arg=c("bot", "neut", "decl"))

mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.1, method = "mnlogistic")
summary(mod_prob)
mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.1, method = "rejection")
summary(mod_prob)


res.gfit.bott <- gfit(target=obs_stats, sumstat= sims_stats[models == "neut",], statistic=mean, nb.replicate=100)
plot(res.gfit.bott, main="Histogram under H0")

summary(res.gfit.bott)

# or PCA

gfitpca(target=obs_stats, sumstat=sims_stats, index=models, cprob=.1)

## 

mylabels <- c("num_alleles_mean","mean_allele_size_sd", "exp_het_mean", "obs_het_mean", "het_ratio")
par(mfrow = c(1,5), mar=c(5,2,4,0))
for (i in mylabels){
  hist(sims_stats[models == "decl",i],breaks=40, xlab=i, main="")
  abline(v = obs_stats[i], col = 2)
}



# posterior predictive checks

stat_bot <- subset(sims_stats, subset=models=="bot")
par_bot <- subset(sims_param, subset=models=="bot")

cv.res.rej <- cv4abc(data.frame(Na=par_bot[,"N0"]), stat_bot, nval=10,
              tols=c(.1,.2, .3), method="rejection")

summary(cv.res.rej)


# parameter inference

res <- abc(target = unlist(obs_stats), param = par_bot, sumstat = stat_bot, tol = 0.1, 
           transf=c("log"), method="loclinear")
summary(res)
hist(res)


par(cex=.8)
plot(res, param=par_bot)



