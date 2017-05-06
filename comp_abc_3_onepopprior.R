#### plot the calculated posteriors
load("abc_estimates/abc_full.RData")

abc_full[[1]]

## load model probabilities
model_probs <- read.table(paste0("results/model_probs/", "onepopprior_1mio_sims_model_selection.txt"))
model_bot <- model_probs$bot > 0.5

abc_bot <- list(abc_full[[1]][model_bot, ], abc_full[[2]][model_bot])


par(mfrow = c(4, 7), mar=c(4,4,1,1))

for (i in 1:17) {
  hist(abc_bot[[2]][[i]]$adj.values, 
       main = abc_bot[[1]][i, 2], 
       xlab = abc_bot[[1]][i, 3],
       #xlim = c(0, 0.005),
       breaks = 100) #caption = abc_5k[[1]][i, 2], ,  xlab = abc_5k[[1]][i, 3]
}



abc_full[[1]]


# load sims
path_to_sims <- paste0("mu_gammaprior_1mio.txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)


### subsetting and definition of model selection parameters

# all sumstats
all_sumstats <- all_sumstats_full
# seal descriptive and summary data 
all_seals <- all_seals_full
# parameter columns in simulation data.frame
params <- c(1:12)
# create a character vector with models
models <- sims$model
# tolerance rate
tol <- 0.0005
# extract names of all models
model_names <- names(table(models))
# divide stats and parameters
sims_stats <- sims[sumstats] 
sims_param <- sims[params]

# just use bottleneck model
mod <- "bot"

# subset parameters and summary statistics from simulations under the bottleneck model
stat_mod <- subset(sims_stats, subset=models==mod)
par_mod <- subset(sims_param, subset=models==mod)

par_mod <- par_mod %>% 
  mutate(N_bot = N0 * N_bot) %>%
  mutate(N_hist_bot = N0 * N_hist_bot) %>%
  mutate(start_bot = 4 * N0 * start_bot) %>%
  mutate(end_bot = 4 * N0 * end_bot) 

plot(abc_bot[[2]][[13]], par_mod$N_bot)
