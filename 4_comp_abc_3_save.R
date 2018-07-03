# save abc estimates to RData file
# replace filename neut with bot when working on bottleneck model estimates.

load("abc_estimates/abc_sims_1000k_4mods_neutral_30.RData")
abc_full[[1]]
length(abc_full[[2]])

# aim: save everything in one data.frame
abc_complete <- abc_full[[1]]

# save adjusted and unadjusted values to list
all_vals_adj <- list()
all_vals_unadj <- list()
for (i in 1:length(abc_full[[2]])) {
  all_vals_adj <- c(all_vals_adj, list(abc_full[[2]][[i]]$adj.values))
  all_vals_unadj <- c(all_vals_unadj, list(abc_full[[2]][[i]]$unadj.values))
}
  
# add values to abc_complete
abc_complete$adj_vals <- all_vals_adj
abc_complete$unadj_vals <- all_vals_unadj

save(abc_complete, file = paste0("abc_estimates/abc_1000k_4mods_neutral_complete_30.RData"))
# 
# library(tidyr)
# abc_bot <- unnest(abc_complete)
# library(dplyr)
# library(ggplot2)
# library(ggbeeswarm)
# abc_bot %>% filter(pars == "nbot") %>% 
#   group_by(species) %>% 
#   #sample_n(4000) %>% 
#   #filter(!(species == "hawaiian_monk_seal" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal")) %>% 
#   ggplot(aes(y = adj_vals, x = species)) + 
#   # geom_quasirandom(alpha = 0.01, size = 2, color = "#053061", width = 0.45, bandwidth = 1.5) +
#   geom_quasirandom(alpha = 0.05, size = 1, color = "#053061", width = 0.47, bandwidth = 2.5) +
#   # geom_beeswarm(priority='density', alpha = 0.5, cex = 0.2, color = "grey") +
#   # geom_jitter(size = 0.5, alpha = 0.1, width = 0.2, color = "grey") +
#   geom_boxplot(width = 0.4, outlier.shape = NA, color = "white", alpha = 0.5, size = 0.2) +
#   stat_summary(fun.y = "estimate_mode", colour = "black", geom = "point", size = 2, shape = 21, fill = "grey") +
#   #stat_summary(fun.y = "mean", colour = "blue", geom = "point") +
#   #theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
#   # scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
#   #ylim(0, 900) +
#   # scale_x_discrete(labels = species_names_bot_twolines) + 
#   scale_y_continuous(breaks = c(seq(from = 0, to = 800, by = 200)), limits = c(0,900)) +
#   #scale_y_discrete(expand = c(0.01, 0))+
#   xlab("") +
#   ylab(expression(Bottleneck~N[e])) +
#   coord_flip() +
#   theme(plot.margin = unit(c(0.5,1.5,0.5,0), "cm"),
#         axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 11),
#         axis.title.y = element_text(size = 13),
#         axis.text.y = element_text(size = 11))

