
## community feasibility as a function of different community metrics

library(tidyverse)
library(broom)
library(broom.mixed)
library(cowplot)
library(corrplot)
library(ggfortify)
library(DHARMa)
library(colorblindr)
library(effects)
library(patchwork)
library(lmerTest)

# -------------------------------------------------------------------------
# read data

vers <- ""
# vers <- "_v2"

# feasibility domain of the observed data
fd.obs <- read.csv2("results/feasibility_domain_observed.csv")
# forgot an "s"
names(fd.obs)[which(names(fd.obs) == "guild")] <- "guilds"

# observed network metrics
metrics.obs <- read.csv2("results/community_metrics.csv")
# null network metrics - redo this again with the 1e3 replicates
# these are in wide format to have fewer rows
metrics.null <- read_csv2("results/community_metrics_null.csv")
# -------------------------------------------------------------------------
# 1 - obtain z-scores for the observed metric values against the null distribution

# 1.1 metric means and sd

metrics.null.long <- metrics.null %>% 
  # number of diag. dom. species is not normally distributed
  select(-diagonally_dominant_sp,-complexity,-skewness,-intra_inter_ratio) %>%
  pivot_longer(richness:quant_modularity,names_to = "metric",values_to = "value") %>%
  group_by(intraguild.type,guilds,null.model,metric) %>%
  summarise(mean.value = mean(value, na.rm = T),
            sd.value = sd(value, na.rm = T)) 
# fuck it
metrics.null.long$null.model[metrics.null.long$null.model == "diagonal dominance"] <- "diag.dom"

metrics.null.wide <- metrics.null.long %>% 
  select(metric,intraguild.type,guilds,null.model,mean.value,sd.value) %>%
  pivot_wider(names_from = null.model,values_from =c(mean.value,sd.value))

metrics.obs.z <- left_join(metrics.obs,metrics.null.wide) %>%
  na.omit() %>% mutate(dd = (value-mean.value_diag.dom)/sd.value_diag.dom,
                       topo = (value-mean.value_topology)/sd.value_topology) %>%
  select(-(mean.value_diag.dom:sd.value_topology)) %>%
  pivot_longer(cols = c(dd,topo),names_to = "null_model",values_to = "z_score") %>%
  mutate(z_score = replace(z_score, z_score > 1e3 | z_score < -1e3, NA))

# -------------------------------------------------------------------------
# plot z-scores
pd <- .32
standard.z <- ggplot(metrics.obs.z, aes(x = intraguild.type, y = z_score)) + 
  geom_point(aes(color = intraguild.type, shape = null_model),
             # shape = 21, 
             position = position_dodge(pd), 
             alpha = .25,
             size = .85) +
  # geom_errorbar(data = fd.means, aes(x = intraguild.type, y = fd.avg,
  #                                    ymin = fd.avg - fd.sd,
  #                                    ymax = fd.avg + fd.sd,
  #                                    color = intraguild.type,
  #                                    linetype = type),
  #               width = 0.15,
  #               position = position_dodge(pd))+
  # geom_point(data = fd.means, aes(x = intraguild.type,
  #                                 y = fd.avg,
  #                                 fill = intraguild.type,
  #                                 shape = type),
  #            # shape = 21, 
  #            stroke = .5, 
  #            size = 2.5,
  #            position = position_dodge(pd))+
  # geom_boxplot(aes(fill = type),alpha = .6) + 
  scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  scale_shape_manual(values = c(21,23)) +
  facet_grid(metric~guilds, scales = "free_y")+
  labs(x = "",y = "z-score") +
  # xlim(0,0.26) +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  # scale_y_discrete(breaks=NULL, limits=rev) +
  # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
  guides(color=guide_legend(override.aes = list(shape=21))) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  NULL
# standard.z

# -------------------------------------------------------------------------

# final subset
metrics.z <- subset(metrics.obs.z,metric %in% c("avg_diagonal_dominance",
                                                "qual_modularity") &
                      null_model == "topo" &
                      guilds == "all" &
                      intraguild.type == "nesting_larvae_phenology")

metrics.fd <- left_join(metrics.z,fd.obs) %>%
  select(year,plot,omega_mean,omega_isotropic,isotropy_index_mean,metric,z_score) %>%
  pivot_wider(names_from = metric,values_from = z_score) %>%
  group_by(year,plot) %>%
  mutate(net = cur_group_id())

# -------------------------------------------------------------------------
# data visualization

# hist(metrics.fd$avg_diagonal_dominance)
# the modularity is skewed
# hist(metrics.fd$qual_modularity)

# no relationship among covariates
# ggplot(metrics.fd, aes(x = avg_diagonal_dominance,y = qual_modularity)) + 
#   geom_point() + 
#   NULL

# response variable looks ok
# hist(metrics.fd$omega_mean)


# -------------------------------------------------------------------------
# models

# a simple linear model and one with a random factor look similar
m.simple <- lm(omega_mean ~ avg_diagonal_dominance + qual_modularity,data = metrics.fd)
m1 <- lmerTest::lmer(omega_mean ~ avg_diagonal_dominance + qual_modularity + (1|plot),data = metrics.fd)
broom.mixed::tidy(m1)
# the residuals look good!!
DHARMa::testResiduals(m1)

# -------------------------------------------------------------------------




