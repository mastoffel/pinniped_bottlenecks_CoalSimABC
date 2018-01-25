# plotting summary statistics.

# boxplots of sims
library(dplyr)
library(tidyr)
library(ggplot2)

sims_test <- sims %>% 
        group_by(model) %>% 
        sample_n(size = 1000) %>% 
        select(num_alleles_mean, exp_het_mean, mean_allele_range, prop_low_afs_mean, mratio_mean) %>% 
        gather(-model, key = "statistic", value = "value")


ggplot(data = sims_test, aes(statistic, value)) + 
    geom_point() +
    geom_boxplot()

df <- data.frame(x = c(1,2,3,4), y = c(1,2,3,4))
plot(df$x, df$y)
