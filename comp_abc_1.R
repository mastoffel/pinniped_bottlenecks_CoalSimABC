# comparative analysis of abc output 1

# analyse simulated microsatellites under three different models000
library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)
library(data.table)
library(reshape2)
library("abctools")
library(abc)
library(ggplot2)

path_to_sims <- "sims_simple_pop100k_sim100k.txt"

# load data 
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# which stats to use
names(sims)
sumstats <- c("num_alleles_mean",  "num_alleles_sd" , "mratio_mean",  "mratio_sd",
              "prop_low_afs_mean", "prop_low_afs_sd")

# parameter columns
params <- c(1:12)
# character vector with models
models <- sims$model
# tolerance rate
tol <- 0.001
# cross-validation replicates
cv_rep <- 100
# method
method <- "mnlogistic"

# some processing

# extract names of all models
all_models <- names(table(models))

# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]


##### first visual checks 
# check whether sumstats are different across models
par(mfcol = c(4, 3), mar=c(5,5,1,1))
for (i in sumstats){
  boxplot(sims[[i]] ~ models, main = i)
}

### can abc at all distinguish between the 4 models ?
cv.modsel <- cv4postpr(models, sims_stats, nval=cv_rep, tol=tol, method=method)
s <- summary(cv.modsel)
png("plots/confusion_matrix.png", width=4, height=4, units="in", res=300)
plot(cv.modsel, names.arg= all_models)
dev.off() #only 129kb in size

# load all seal data
all_seals <- sealABC::read_excel_sheets("../data/seal_data_largest_clust_and_pop.xlsx")

# seal species
# genotypes <- all_seals[[19]]
# genotypes <- genotypes[4:ncol(genotypes)]

all_sumstats <- lapply(all_seals, function(x) mssumstats(x[4:ncol(x)], type = "microsats", data_type = "empirical"))
all_sumstats <- do.call(rbind, all_sumstats)



### model selection

abc_mod_probs <- function(genotypes, models, sims_stats, tol, method, sumstats){
  genotypes <- genotypes[4:ncol(genotypes)]
  # calculate empirical summary statistics
  obs_stats <- sealABC::mssumstats(genotypes, type = "microsats", data_type = "empirical")
  obs_stats <- obs_stats[sumstats]
  # model probabilities
  mod_prob <- abc::postpr(obs_stats, models, sims_stats, tol = tol, method = method)
  # sum_prob <- summary(mod_prob)
}

#check probabilites for all species
cl <- makeCluster(getOption("cl.cores", detectCores()-5))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_probs <- parallel::parLapply(cl, all_seals[1:28], abc_mod_probs, models, sims_stats, tol, method, sumstats)
stopCluster(cl)

all_probs_df <- do.call(rbind, lapply(all_probs, function(x) {
  # sum_abc <- summary(x)
  # if (is.null(sum_abc$rejection$Prop)){
  #   out <- sum_abc$Prob
  # } else {
  #   out <- sum_abc$rejection$Prop
  # }
  out <- round(x$pred, 3)
  out
} ))
write.table(all_probs_df, file = "output/model_selection.txt")




#### Does the preferred model provide a good fit to the data?  
calc_fit <- function(species_name, model, sumstats, all_sumstats, sims_stats) {
  res_gfit_bot <- abc::gfit(target = all_sumstats[species_name, sumstats], sumstat = sims_stats[models == model, ], 
                       nb.replicate = cv_rep, tol = tol)
}

# get all names
all_names <- names(all_seals)[1:28]

# calculate all fits
cl <- makeCluster(getOption("cl.cores", detectCores()-5))
clusterEvalQ(cl, c(library("sealABC"), library("abc")))
all_fits_bot <- parLapply(cl, all_names, calc_fit, "bot", sumstats, all_sumstats, sims_stats)
all_fits_neut <- parLapply(cl, all_names, calc_fit, "neut", sumstats, all_sumstats, sims_stats)
stopCluster(cl)

# plots
pdf(file = "plots/goodnessoffit.pdf")
par(mfrow = c(8, 4), mar=c(3,3,1,1))
sapply(1:16, function(x) {
  plot(all_fits_bot[[x]], main = all_names[x])
  plot(all_fits_neut[[x]], main = all_names[x])
})
dev.off()

pdf(file = "plots/goodnessoffit2.pdf")
par(mfrow = c(6, 4), mar=c(3,3,1,1))
sapply(17:28, function(x) {
  plot(all_fits_bot[[x]], main = all_names[x])
  plot(all_fits_neut[[x]], main = all_names[x])
})
dev.off()


# save the p_values and distances
summary(res_gfit_bot)$pvalue

all_p_bot <- unlist(lapply(all_fits_bot, function(x) summary(x)$pvalue))
all_p_neut <-  unlist(lapply(all_fits_neut, function(x) summary(x)$pvalue))

p_df <- data.frame(species = all_names, bot_p = all_p_bot, neut_p = all_p_neut)
write.table(p_df, file = "output/p_vals_fit.txt", row.names = FALSE)