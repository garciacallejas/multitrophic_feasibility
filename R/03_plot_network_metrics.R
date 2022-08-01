
# plot network metrics

library(tidyverse)
library(cowplot)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------

# vers <- "_v2"
vers <- ""

metrics.observed <- read.csv2(file = paste("results/community_metrics",vers,".csv",sep=""),
                         header = TRUE,
                         stringsAsFactors = FALSE)

metrics.null <- read.csv2(file = paste("results/community_metrics_null",vers,".csv",sep=""),
                     header = TRUE,
                     stringsAsFactors = FALSE)

# metrics.mf <- read.csv2(file = paste("results/community_metrics_mean_field",vers,".csv",sep=""),
#                               header = TRUE,
#                               stringsAsFactors = FALSE)

# join data ---------------------------------------------------------------

metrics.observed$sd <- 0
metrics.observed$type <- "Observed"

# metrics.mf$sd <- 0
# metrics.mf$type <- "Mean field"

# metrics.null.2 <- metrics.null %>% group_by(year,plot,guilds,metric) %>%
#   summarise(metrics.average = mean(value),
#             sd = sd(value))
# metrics.null.2$type <- "null"
# 
# names(metrics.null.2)[which(names(metrics.null.2) == "metrics.average")] <- "value"

# metrics.all <- bind_rows(metrics.observed,metrics.null.2,metrics.mf)
metrics.all <- metrics.observed

metrics.all$guilds <- factor(metrics.all$guilds,levels = c("plants",
                                                 "floral visitors",
                                                 "herbivores",
                                                 "plants-floral visitors",
                                                 "plants-herbivores","all"))

metrics.all$guilds <- recode(metrics.all$guilds, "plants" = "Plants", 
                             "floral visitors" = "Pollinators",
                             "herbivores" = "Herbivores",
                             "plants-floral visitors" = "Plants-Pollinators",
                             "plants-herbivores" = "Plants-Herbivores",
                             "all" = "All")

# metrics.all$type <- recode(metrics.all$type, "observed" = "Observed","null" = "Randomized",
#                            "mean field" = "Mean field")
# metrics.all$type <- factor(metrics.all$type,levels = c("Observed","Randomized","Mean field"))

# subset metrics ----------------------------------------------------------

metrics.all <- subset(metrics.all,metric %in% c("richness","connectance",
                                                "diagonally_dominant_sp",
                                                "avg_diagonal_dominance",
                                                "degree_distribution",
                                                # "kurtosis",
                                                "skewness",
                                                "complexity",
                                                "qual_modularity",
                                                "quant_modularity"))

metrics.all$metric <- recode(metrics.all$metric, "richness" = "Richness", 
                             "connectance" = "Connectance",
                             "intra_inter_ratio" = "Intra/Intersp. ratio",
                             "diagonally_dominant_sp" = "Diagonally dominant sp.",
                             "avg_diagonal_dominance" = "Mean diagonal dominance",
                             "degree_distribution" = "Degree dist. het.",
                             "skewness" = "Skewness",
                             "complexity" = "Complexity",
                             "qual_modularity" = "Modularity (Qualitative)",
                             "quant_modularity" = "Modularity (Quantitative)")

# plot --------------------------------------------------------------------

p2 <- ggplot(metrics.all, aes(x = plot, y = value, group = interaction(year,intraguild.type))) + 
  geom_errorbar(aes(ymin = value-sd, ymax = value+sd,
                    color = intraguild.type),
                position = position_dodge(.2),
                width = 0.2)+
  # geom_col(aes(fill = guilds)) + 
  geom_line(aes(color = intraguild.type),position = position_dodge(.2)) +
  geom_point(aes(fill = intraguild.type),position = position_dodge(.2),shape = 21) +
  facet_grid(metric~guilds,scales = "free_y") +
  scale_x_continuous(breaks = 1:9) +
  labs(y = "metric value") +
  scale_color_OkabeIto() + 
  scale_fill_OkabeIto(darken = 0.2) +
  # theme_cowplot() + 
  theme_bw()+
  theme(strip.background = element_blank())+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  NULL
# p2

metric.means <- metrics.all %>% 
  group_by(guilds, intraguild.type, metric) %>%
  summarise(mean.value = mean(value),
            sd.value = sd(value),
            se.value = sd(value)/sqrt(n()))

metric.means.v2 <- metric.means
# metric.means.v2$mean.value[metric.means.v2$metric == "Intra/Intersp. ratio" &
#                              metric.means.v2$guilds == "All"] <- NA
# metric.means.v2$sd.value[metric.means.v2$metric == "Intra/Intersp. ratio" &
#                              metric.means.v2$guilds == "All"] <- NA

pd <- 0.5
p3 <- ggplot(metric.means) + 

  geom_jitter(data = metrics.all, aes(x = intraguild.type, 
                                     y = value, 
                                     fill = intraguild.type), 
             shape = 21, alpha = .3, size = .75,width = .2) +
  
  geom_errorbar(aes(x = intraguild.type, y = mean.value,
                    ymin = mean.value - sd.value,
                    ymax = mean.value + sd.value,
                    # ymin = mean.value - se.value,
                    # ymax = mean.value + se.value,
                    color = intraguild.type),
                width = 0.15,
                position = position_dodge(pd))+
  geom_point(aes(x = intraguild.type, 
                 y = mean.value, 
                 fill = intraguild.type),
             shape = 21, stroke = .5, size = 2,
             position = position_dodge(pd))+
  
  
  scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  facet_grid(metric~guilds, scales = "free_y")+
  labs(x = "",y = "Metric value") +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  scale_x_discrete(breaks=NULL) +
  NULL
# p3

# m.all.ii <- subset(metric.means, guilds == "All" & metric == "Intra/Intersp. ratio")

# p3.2 <- ggplot(m.all.ii) + 
#   
#   geom_errorbar(aes(x = type, y = mean.value,
#                     ymin = mean.value - sd.value,
#                     ymax = mean.value + sd.value,
#                     # ymin = mean.value - se.value,
#                     # ymax = mean.value + se.value,
#                     color = type),
#                 width = 0.15,
#                 position = position_dodge(pd))+
#   geom_point(aes(x = type, 
#                  y = mean.value, 
#                  fill = type),
#              shape = 21, stroke = .5, size = 2.5,
#              position = position_dodge(pd))+
#   scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
#   scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
#   facet_grid(metric~guilds, scales = "free_y")+
#   labs(x = "",y = "Metric value") +
#   theme_bw()+
#   # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
#   theme(strip.background = element_blank())+
#   theme(legend.position="none")+
#   # scale_y_continuous(position = "right") +
#   scale_x_discrete(breaks=NULL) +
#   NULL
# p3.2
# 
# library(patchwork)
# p3.3 <- p3/(p3.2 + plot_spacer() + plot_spacer()) + plot_layout(heights = c(4,1))

# store plot --------------------------------------------------------------

  # ggsave(filename = paste("results/images/Fig_S1",vers,".pdf",sep=""),plot = p3.3,
  #        device = cairo_pdf,
  #        width = 9,height = 10,dpi = 300)

# -------------------------------------------------------------------------

# distribution of null values against observed

metrics.null.long <- metrics.null %>% select(-intra_inter_ratio) %>% 
  pivot_longer(richness:complexity,
               names_to = "metric",
               values_to = "value")



# metrics.observed <- metrics.observed %>% filter(metric != "intra_inter_ratio") %>% 
#   group_by(guilds,intraguild.type,metric) %>%
#   mutate(scaled.value = scales::rescale(value))
#   
# metrics.null.long <- metrics.null.long %>% filter(metric != "intra_inter_ratio") %>% 
#   group_by(guilds,intraguild.type,metric,null.model) %>%
#   mutate(scaled.value = scales::rescale(value))

# p1 <- ggplot(metrics.null.all, aes(x = value)) + 
#   geom_density(aes(color = intraguild.type, linetype = null.model)) + 
#   geom_vline(data = metrics.obs.all, aes(xintercept = value, color = intraguild.type)) +
#   facet_wrap(metric~guilds, scales = "free") +
#   theme_bw() +
#   NULL
# p1

# -------------------------------------------------------------------------
# subset for the final figures

metrics.observed$intraguild.type <- recode(metrics.observed$intraguild.type, 
                                           "mean_field" = "Mean field", 
                                           "nesting_larvae" = "Resources",
                                           "nesting_larvae_phenology" = "Resources and \nphenology")

metrics.obs.all <- metrics.observed %>% filter(guilds == "all" & 
                                                 metric %in% c("avg_diagonal_dominance",
                                                               "diagonally_dominant_sp",
                                                               "qual_modularity",
                                                               "quant_modularity"))
metrics.null.long$intraguild.type <- recode(metrics.null.long$intraguild.type, 
                                           "mean_field" = "Mean field", 
                                           "nesting_larvae" = "Resources",
                                           "nesting_larvae_phenology" = "Resources and \nphenology")

metrics.null.all <- metrics.null.long %>% filter(guilds == "all" & 
                                                   metric %in% c("avg_diagonal_dominance",
                                                                 "diagonally_dominant_sp",
                                                                 "qual_modularity",
                                                                 "quant_modularity") &
                                                   null.model == "topology")
# -------------------------------------------------------------------------
# average diagonal dominance
avg.dd.obs <- metrics.obs.all %>% filter(metric == "avg_diagonal_dominance")
avg.dd.null <- metrics.null.all %>% filter(metric == "avg_diagonal_dominance")

pdd <- ggplot(avg.dd.null, aes(x = value)) + 
  geom_density(aes(color = intraguild.type, fill = intraguild.type), alpha = .4) + 
  geom_vline(data = avg.dd.obs, aes(xintercept = value, color = intraguild.type)) +
  facet_grid(intraguild.type~.) +
  scale_color_OkabeIto(darken = 0.2) +
  scale_fill_OkabeIto(darken = 0.2) +
  theme_bw() +
  theme(strip.background = element_blank())+
  theme(strip.text.y = element_blank()) +
  theme(legend.position="none")+
  labs(x = "average diagonal dominance",y = "") +
  NULL
# pdd

# -------------------------------------------------------------------------
# diagonally dominant species
# avg.ds.obs <- metrics.obs.all %>% filter(metric == "diagonally_dominant_sp")
# avg.ds.null <- metrics.null.all %>% filter(metric == "diagonally_dominant_sp")
# 
# pds <- ggplot(avg.ds.null, aes(x = value)) +
#   geom_density(aes(color = intraguild.type, fill = intraguild.type), alpha = .4) +
#   geom_vline(data = avg.ds.obs, aes(xintercept = value, color = intraguild.type)) +
#   facet_grid(intraguild.type~.) +
#   theme_bw() +
#   NULL
# pds

# -------------------------------------------------------------------------
# modularity (qualitative - i.e. accounting only for topology)
avg.m.obs <- metrics.obs.all %>% filter(metric == "qual_modularity")
avg.m.null <- metrics.null.all %>% filter(metric == "qual_modularity")

pm <- ggplot(avg.m.null, aes(x = value)) + 
  geom_density(aes(color = intraguild.type, fill = intraguild.type), alpha = .4) + 
  geom_vline(data = avg.m.obs, aes(xintercept = value, color = intraguild.type)) +
  facet_grid(intraguild.type~.) +
  scale_color_OkabeIto(darken = 0.2) +
  scale_fill_OkabeIto(darken = 0.2) +
  theme_bw() +
  theme(strip.background = element_blank())+
  theme(legend.position="none")+
  labs(x = "modularity",y = "") +
  NULL
# pm

# -------------------------------------------------------------------------
fdist <- pdd+pm

ggsave(filename = paste("results/images/Fig_metrics",vers,".pdf",sep=""),
       plot = fdist,
       device = cairo_pdf,
       width = 9,height = 7,dpi = 300)

# -------------------------------------------------------------------------
# modularity (quantitative - i.e. accounting for topology+interaction strength)
# avg.mq.obs <- metrics.obs.all %>% filter(metric == "quant_modularity")
# avg.mq.null <- metrics.null.all %>% filter(metric == "quant_modularity")
# 
# pmq <- ggplot(avg.mq.null, aes(x = value)) + 
#   geom_density(aes(color = intraguild.type, fill = intraguild.type), alpha = .4) + 
#   geom_vline(data = avg.mq.obs, aes(xintercept = value, color = intraguild.type)) +
#   facet_grid(intraguild.type~null.model, scales = "free_y") +
#   theme_bw() +
#   NULL
# pmq

