
# plot network metrics

library(tidyverse)

# read data ---------------------------------------------------------------

# vers <- "_v2"
vers <- ""

metrics.observed <- read.csv2(file = paste("results/community_metrics_observed",vers,".csv",sep=""),
                         header = TRUE,
                         stringsAsFactors = FALSE)

metrics.null <- read.csv2(file = paste("results/community_metrics_null",vers,".csv",sep=""),
                     header = TRUE,
                     stringsAsFactors = FALSE)

metrics.mf <- read.csv2(file = paste("results/community_metrics_mean_field",vers,".csv",sep=""),
                              header = TRUE,
                              stringsAsFactors = FALSE)

# join data ---------------------------------------------------------------

metrics.observed$sd <- 0
metrics.observed$type <- "Observed"

metrics.mf$sd <- 0
metrics.mf$type <- "Mean field"

metrics.null.2 <- metrics.null %>% group_by(year,plot,guilds,metric) %>%
  summarise(metrics.average = mean(value),
            sd = sd(value))
metrics.null.2$type <- "null"

names(metrics.null.2)[which(names(metrics.null.2) == "metrics.average")] <- "value"

metrics.all <- bind_rows(metrics.observed,metrics.null.2,metrics.mf)
# metrics.all <- metrics.observed

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

metrics.all$type <- recode(metrics.all$type, "observed" = "Observed","null" = "Randomized",
                           "mean field" = "Mean field")
metrics.all$type <- factor(metrics.all$type,levels = c("Observed","Randomized","Mean field"))

# subset metrics ----------------------------------------------------------

metrics.all <- subset(metrics.all,metric %in% c("richness","connectance",
                                                "intra_inter_ratio",
                                                "degree_distribution",
                                                # "kurtosis",
                                                "skewness",
                                                "complexity",
                                                "modularity"))

metrics.all$metric <- recode(metrics.all$metric, "richness" = "Richness", 
                             "connectance" = "Connectance",
                             "intra_inter_ratio" = "Intra/Intersp. ratio",
                             "degree_distribution" = "Degree dist. het.",
                             "skewness" = "Skewness",
                             "complexity" = "Complexity",
                             "modularity" = "Modularity")

# plot --------------------------------------------------------------------

library(cowplot)
library(colorblindr)


p2 <- ggplot(metrics.all, aes(x = plot, y = value, group = type)) + 
  geom_errorbar(aes(ymin = value-sd, ymax = value+sd,
                    color = type),
                position = position_dodge(.2),
                width = 0.2)+
  # geom_col(aes(fill = guilds)) + 
  geom_line(aes(color = type),position = position_dodge(.2)) +
  geom_point(aes(fill = type),position = position_dodge(.2),shape = 21) +
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
  group_by(guilds, type, metric) %>%
  summarise(mean.value = mean(value),
            sd.value = sd(value),
            se.value = sd(value)/sqrt(n()))

metric.means.v2 <- metric.means
metric.means.v2$mean.value[metric.means.v2$metric == "Intra/Intersp. ratio" &
                             metric.means.v2$guilds == "All"] <- NA
metric.means.v2$sd.value[metric.means.v2$metric == "Intra/Intersp. ratio" &
                             metric.means.v2$guilds == "All"] <- NA

pd <- 0.5
p3 <- ggplot(metric.means.v2) + 

  geom_errorbar(aes(x = type, y = mean.value,
                    ymin = mean.value - sd.value,
                    ymax = mean.value + sd.value,
                    # ymin = mean.value - se.value,
                    # ymax = mean.value + se.value,
                    color = type),
                width = 0.15,
                position = position_dodge(pd))+
  geom_point(aes(x = type, 
                 y = mean.value, 
                 fill = type),
             shape = 21, stroke = .5, size = 2.5,
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
p3

m.all.ii <- subset(metric.means, guilds == "All" & metric == "Intra/Intersp. ratio")

p3.2 <- ggplot(m.all.ii) + 
  
  geom_errorbar(aes(x = type, y = mean.value,
                    ymin = mean.value - sd.value,
                    ymax = mean.value + sd.value,
                    # ymin = mean.value - se.value,
                    # ymax = mean.value + se.value,
                    color = type),
                width = 0.15,
                position = position_dodge(pd))+
  geom_point(aes(x = type, 
                 y = mean.value, 
                 fill = type),
             shape = 21, stroke = .5, size = 2.5,
             position = position_dodge(pd))+
  scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  facet_grid(metric~guilds, scales = "free_y")+
  labs(x = "",y = "Metric value") +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  theme(legend.position="none")+
  # scale_y_continuous(position = "right") +
  scale_x_discrete(breaks=NULL) +
  NULL
p3.2

library(patchwork)
p3.3 <- p3/(p3.2 + plot_spacer() + plot_spacer()) + plot_layout(heights = c(4,1))

# store plot --------------------------------------------------------------

ggsave(filename = paste("results/images/Fig_S1",vers,".pdf",sep=""),plot = p3.3,
       device = cairo_pdf,
       width = 9,height = 10,dpi = 300)
