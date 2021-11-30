
## community feasibility as a function of different community metrics

library(tidyverse)
library(broom)
library(cowplot)
library(corrplot)
library(ggfortify)
library(DHARMa)
library(colorblindr)
library(effects)
library(patchwork)

# -------------------------------------------------------------------------
# read data

vers <- ""
# vers <- "_v2"

obs.feas <- read.csv2(here::here("results",paste("feasibility_domain_observed",vers,".csv",sep="")))
obs.metrics <- read.csv2(here::here("results",paste("community_metrics_observed",vers,".csv",sep="")))

names(obs.feas)[3] <- "community"
names(obs.metrics)[3] <- "community"

# join data ---------------------------------------------------------------

obs.d1 <- left_join(obs.feas,obs.metrics)

# expand
obs.data <- pivot_wider(obs.d1,names_from = metric,values_from = value) %>%
  filter(!is.na(fd.average))

obs.data <- obs.data[,c("fd.average","community",
                        "richness","connectance",
                        "intra_inter_ratio","complexity",
                        "degree_distribution","skewness","modularity")]

obs.data <- subset(obs.data,!is.infinite(intra_inter_ratio))

# PCA on network metrics --------------------------------------------------
# also standard correlations

pca.data <- subset(obs.data, community == "all")

pca.data.names <- pca.data
names(pca.data.names)[5] <- "intra/inter ratio"
names(pca.data.names)[7] <- "degree distribution"

# SKEWNESS is maximal (-1) for all multitrophic observed communities
# so, drop it
pca.net <- pca.data.names %>% 
  # filter(community == "all") %>%
  dplyr::select(richness,connectance,
         complexity,
         'degree distribution',
         'intra/inter ratio',
         # skewness,
         modularity) %>%
  prcomp(scale = TRUE)


# nicer plot with ggfortify
pca.plot <- autoplot(pca.net, data = pca.data.names, fill = 'fd.average', 
                     shape = 21,
                     size = 3,
                     loadings = TRUE, 
                     loadings.colour = 'darkgrey',
                     loadings.label = TRUE, 
                     loadings.label.size = 3.5,
                     loadings.label.colour = "darkred") +
  theme_bw() + 
  scale_fill_viridis_b(name = "Feasibility\ndomain") +
  NULL
# pca.plot

# ggsave(filename = "results/images/Fig_S2.pdf",plot = pca.plot,
#        device = cairo_pdf,
#        width = 5,height = 4,dpi = 300)

# plot variance explained by each component
# pca.net %>%
#   tidy(matrix = "d") %>%
#   ggplot(aes(PC, percent)) +
#   geom_col(fill = "#56B4E9", alpha = 0.8) +
#   scale_x_continuous(breaks = 1:9) +
#   scale_y_continuous(
#     labels = scales::percent_format(),
#     expand = expansion(mult = c(0, 0.01))
#   ) +
#   theme_minimal_hgrid(12)

# model -------------------------------------------------------------------

my.data <- pca.data
my.data$community <- as.factor(my.data$community)

my.d3 <- my.data
my.d3$fd.average <- my.d3$fd.average + 0.01

# first visualization

# plot.data <- my.d3[,c("fd.average","connectance","degree_distribution")]
# plot.data <- pivot_longer(plot.data,connectance:degree_distribution,
#                           names_to = "metric",
#                           values_to = "value")
# # names(plot.data)[1] <- "feasibility_domain"
# 
# p1 <- ggplot(plot.data, aes(x = value, y = fd.average)) +
#   geom_point(show.legend = FALSE, size = 2) +
#   facet_grid(.~metric, scales = "free_x")+
#   theme_bw() +
#   theme(strip.background = element_blank())+
#   ylab("feasibility domain") + xlab("") +
#   NULL
# p1

my.d4 <- my.d3
my.d4$connectance <- c(scale(my.d4$connectance))
my.d4$degree_distribution <- c(scale(my.d4$degree_distribution))

m1 <- glm(fd.average ~ connectance + degree_distribution,
            data = my.d4,
            family = Gamma(link = "log"))

# ok
# simulateResiduals(m1,plot = TRUE)

# Table 1 from here
# summary(m1)

# -------------------------------------------------------------------------
ec <- effects::effect("connectance",m1,xlevels = 20)
ed <- effects::effect("degree_distribution",m1,xlevels = 20)
# plot(ec)
# plot(ed)

ec.df <- data.frame(ec)
ed.df <- data.frame(ed)

ec.df$connectance.unscaled <- ec.df$connectance * sd(my.d3$connectance) + mean(my.d3$connectance)
ed.df$dd.unscaled <- ed.df$degree_distribution * sd(my.d3$degree_distribution) + mean(my.d3$degree_distribution)

cp <- ggplot(ec.df,aes(y = fit, x = connectance.unscaled)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = .3) +
  geom_line()+
  geom_point(data = my.d3,aes(x = connectance, y = fd.average))+
  theme_bw() +
  labs(x="Connectance",y = "Feasibility domain") +
  ggtitle(label = "A)") +
  NULL
# cp

ddp <- ggplot(ed.df,aes(y = fit, x = dd.unscaled)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = .3) +
  geom_line()+
  geom_point(data = my.d3,aes(x = degree_distribution, y = fd.average))+
  theme_bw() +
  labs(x="Gini index of the degree distribution",y = "Feasibility domain") +
  ggtitle(label = "B)") +
  NULL
# ddp

metrics.plot <- cp + ddp

# ggsave(filename = "results/images/Fig_3.pdf",plot = metrics.plot,
#        device = cairo_pdf,
#        width = 8,height = 4,dpi = 300)

# -------------------------------------------------------------------------
# 3d scatterplot

library(plotly)

#Setup Axis
axis_x <- seq(min(my.d4$connectance), max(my.d4$connectance), length.out = 5)
axis_y <- seq(min(my.d4$degree_distribution), max(my.d4$degree_distribution), length.out = 5)
axis_x_uns <- axis_x * sd(my.d3$connectance) + mean(my.d3$connectance)
axis_y_uns <- axis_y * sd(my.d3$degree_distribution) + mean(my.d3$degree_distribution)

#Sample points
lm_surface <- expand.grid(connectance = axis_x,degree_distribution = axis_y,KEEP.OUT.ATTRS = F)
lm_surface$fd.average <- predict(m1, newdata = lm_surface, type = "response")
names(lm_surface)[1:2] <- c("c_scaled","dd_scaled")
lm_surface$connectance <- lm_surface$c_scaled * sd(my.d3$connectance) + mean(my.d3$connectance)
lm_surface$degree_distribution <- lm_surface$dd_scaled * sd(my.d3$degree_distribution) + mean(my.d3$degree_distribution)

lm_surface2 <- matrix(xtabs(fd.average ~ connectance + degree_distribution, data = lm_surface),nrow = 5, ncol = 5)

axx <- list(
  title = 'Connectance',
  nticks = 5,
  range = c(0.35,0.55)
)

axy <- list(
  title = 'Het. Deg. Dist.',
  nticks = 4,
  range = c(0.2,0.35)
)

axz <- list(
  title = 'Feasibility Domain',
  tickvals = c(0,0.05,0.1,0.15,0.2,0.25),  
  range = list(0, 0.27)
)

fig <- plot_ly(x = ~axis_x_uns, 
               y = ~axis_y_uns, 
               z = ~lm_surface2,
               type="surface", opacity = .7, showscale = FALSE) %>%
  add_trace(data=my.d3, x=~connectance, y=~degree_distribution, z=~fd.average, 
            mode="markers", type="scatter3d", opacity = 1, size = 1, color = "darkorange") %>%
            # marker = list(color=cols, opacity=0.7, symbol=105)
            # ) %>%
  add_markers() %>% 
  layout(scene = list(xaxis = axx,
                      yaxis = axy,
                      zaxis = axz),
         showlegend = FALSE)
# fig

# -------------------------------------------------------------------------
# are connectance and degree distribution correlated?
cd.test <- cor.test(my.d3$connectance,my.d3$degree_distribution)
