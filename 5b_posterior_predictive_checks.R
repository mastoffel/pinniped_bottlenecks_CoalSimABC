# evaluate post.predictive checks

library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
# simulated based on posteriors
all_checks <- read_delim("model_evaluation/check5_postpred/sims_10000kbot500_post_pred_checks2.txt", delim = " ")
head(all_checks)

# compare these summary statistics
sumstats <- c("num_alleles_mean", 
              "exp_het_mean", 
              "mratio_mean", 
              "prop_low_afs_mean",
              "mean_allele_range")

# empirical summary stats
all_sumstats_full <- read_delim("data/all_sumstats_40ind_30.txt", delim = " ") %>% 
  select(species, sumstats)

# long format for plotting
all_checks_long <- all_checks %>% 
                      select(species, sumstats) %>% 
                      gather(sumstat, value, -species) %>% 
                      group_by(species, sumstat) 

# observed sumstats long format
all_sumstats_full_long <- all_sumstats_full %>% 
  gather(sumstat, value, -species)

p <- ggplot(all_checks_long , aes(value)) +
    geom_histogram() +
    geom_vline(aes(xintercept = value), all_sumstats_full_long) +
    facet_grid(species ~ sumstat, scales = "free")

ggsave(filename = "post_pred_checks.jpg", plot = p, width = 5, height = 10)

ggplot(all_checks, aes(value)) +
  geom_histogram() +
  facet_wrap(species ~ sumstat)