# script to concatenate two simulation txt files and save as new one
library(data.table)
library(dplyr)
library(tibble)
sim_name1 <- "sims_1000k_4mods"
sim_name2 <- "sims_1000kbot500_iceage_expand"
output_name <- "sims_1000k_4mods_expand.txt"

### load simulations, stored in main folder atm ###
path_to_sims <- paste0(sim_name1, ".txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

# just take first 2mio
sims <- filter(sims, model %in% c("bottleneck", "neutral")) 

### load simulations, stored in main folder atm ###
path_to_sims <- paste0(sim_name2, ".txt")
sims2 <-fread(path_to_sims, stringsAsFactors = FALSE)
sims2 <- as.data.frame(sims2)

# add growth rate to sims
sims <- sims %>% 
  add_column(growth_rate = sims2$growth_rate, .before = "mut_rate")

sims_all <- rbind(sims, sims2)

write.table(sims_all, file = output_name, row.names = FALSE)

