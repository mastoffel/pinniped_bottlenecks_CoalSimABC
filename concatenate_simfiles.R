# script to concatenate two simulation txt files and save as new one
sim_name1 <- "sims_5000k_large_bot"
sim_name2 <- "sims_5000k_large_bot2"
output_name <- "sims_10000k.txt"

### load simulations, stored in main folder atm ###
path_to_sims <- paste0(sim_name1, ".txt")
sims <-fread(path_to_sims, stringsAsFactors = FALSE)
sims <- as.data.frame(sims)

### load simulations, stored in main folder atm ###
path_to_sims <- paste0(sim_name2, ".txt")
sims2 <-fread(path_to_sims, stringsAsFactors = FALSE)
sims2 <- as.data.frame(sims2)

sims_bot <- rbind(sims[sims$model == "bot", ], sims2[sims2$model == "bot", ])
sims_neut <- rbind(sims[sims$model == "neut", ], sims2[sims2$model == "neut", ])

sims_all <- rbind(sims_bot, sims_neut)

write.table(sims_all, file = output_name, row.names = FALSE)
