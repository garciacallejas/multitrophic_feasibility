
# plot feasibility results

library(tidyverse)
library(colorblindr)

# read data ---------------------------------------------------------------
# vers <- "_v2"
vers <- ""

# temp
# fd.files <- list.files(path = "results",pattern = "fd_*",full.names = T)
# fd.observed <- fd.files %>% map_dfr(read.csv2)

file.name <- "results/feasibility_domain_observed"
null.name <- "results/feasibility_domain_null"
# mean.field.name <- "results/feasibility_domain_mean_field"

file.name <- paste(file.name,vers,".csv",sep="")
null.name <- paste(null.name,vers,".csv",sep="")
# mean.field.name <- paste(mean.field.name,vers,".csv",sep="")

fd.observed <- read.csv2(file = file.name,
                         header = TRUE,
                         stringsAsFactors = FALSE)

# fd.null <- read.csv2(file = null.name,
#                      header = TRUE,
#                      stringsAsFactors = FALSE)
# 
# fd.mf <- read.csv2(file = mean.field.name,
#                      header = TRUE,
#                      stringsAsFactors = FALSE)

# join data ---------------------------------------------------------------

fd.observed$type <- "observed"
# fd.mf$type <- "mean field"
# 
# fd.null.2 <- fd.null %>% group_by(year,plot,guilds) %>%
#   summarise(fd.average.rep = mean(fd.average,na.rm = TRUE),
#             fd.sd.rep = sd(fd.average, na.rm = TRUE))
# fd.null.2$type <- "null"

# to avoid errors in repeating names
# names(fd.null.2)[c(4,5)] <- c("fd.average","fd.sd")

fd.all <- fd.observed
fd.all$fd.average <- fd.all$omega_mean
# fd.all <- bind_rows(fd.observed,fd.null.2,fd.mf)

fd.all$guild <- factor(fd.all$guild,levels = c("plants",
                                                 "floral visitors",
                                                 "herbivores",
                                                 "plants-floral visitors",
                                                 "plants-herbivores","all")) 
fd.all$guild <- recode(fd.all$guild, "plants" = "Plants", 
                        "floral visitors" = "Pollinators",
                        "herbivores" = "Herbivores",
                        "plants-floral visitors" = "Plants-Pollinators",
                        "plants-herbivores" = "Plants-Herbivores",
                        "all" = "All")

fd.all$fd.average[which(is.na(fd.all$fd.average))] <- 0
fd.all$plot <- as.factor(fd.all$plot)
fd.all$type <- as.factor(fd.all$type)

fd.all$type <- recode(fd.all$type, "null" = "Randomized", "observed" = "Observed")
fd.all$type <- factor(fd.all$type,levels = c("Observed","Randomized"))

# plot --------------------------------------------------------------------

library(cowplot)
library(colorblindr)

pd <- .4

fd.all.nona <- subset(fd.all,!is.na(fd.average))

fd.means <- fd.all.nona %>% 
  # group_by(guilds, type) %>%
  group_by(intraguild.type, guild) %>%
  summarise(fd.avg = mean(fd.average,na.rm = TRUE),
            fd.sd = sd(fd.average,na.rm = TRUE),
            fd.se = sd(fd.average,na.rm = TRUE)/sqrt(n()))
  
p3 <- ggplot(fd.all.nona, aes(y = intraguild.type, x = fd.average)) + 
  geom_point(aes(color = intraguild.type),
             # shape = 21, 
             position = position_dodge(.5), 
             alpha = .25,
             size = .85) +
  geom_errorbar(data = fd.means, aes(y = intraguild.type, x = fd.avg,
                                     xmin = fd.avg - fd.se,
                                     xmax = fd.avg + fd.se,
                                     color = intraguild.type),
                width = 0.15,
                position = position_dodge(pd))+
  geom_point(data = fd.means, aes(y = intraguild.type, 
                                  x = fd.avg, 
                                  fill = intraguild.type),
             shape = 21, stroke = .5, size = 2.5,
             position = position_dodge(pd))+
  # geom_boxplot(aes(fill = type),alpha = .6) + 
  scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  facet_grid(guild~.)+
  labs(y = "",x = "Feasibility domain") +
  # xlim(0,0.26) +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  scale_y_discrete(breaks=NULL, limits=rev) +
  # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
  NULL
# p3

# -------------------------------------------------------------------------

ip <- ggplot(fd.all) + 
  geom_point(aes(x = omega_mean, y = isotropy_index_mean, color = intraguild.type)) +
  # geom_errorbar(aes(y = isotropy_index_mean, x = guild,
  #                                    ymin = isotropy_index_lowerCI,
  #                                    ymax = isotropy_index_upperCI,
  #                                    color = intraguild.type)) +
  facet_grid(intraguild.type~.)+
  NULL


# store plot --------------------------------------------------------------

# ggsave(filename = paste("results/images/Fig_2",vers,".pdf",sep=""),plot = p3,
#        device = cairo_pdf,
#        width = 4,height = 7,dpi = 300)
