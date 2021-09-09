
# plot feasibility results

library(tidyverse)

# read data ---------------------------------------------------------------
# vers <- "_v2"
vers <- ""

file.name <- "results/feasibility_domain_observed"
null.name <- "results/feasibility_domain_null"
mean.field.name <- "results/feasibility_domain_mean_field"

file.name <- paste(file.name,vers,".csv",sep="")
null.name <- paste(null.name,vers,".csv",sep="")
mean.field.name <- paste(mean.field.name,vers,".csv",sep="")

fd.observed <- read.csv2(file = file.name,
                         header = TRUE,
                         stringsAsFactors = FALSE)

fd.null <- read.csv2(file = null.name,
                     header = TRUE,
                     stringsAsFactors = FALSE)

fd.mf <- read.csv2(file = mean.field.name,
                     header = TRUE,
                     stringsAsFactors = FALSE)

# join data ---------------------------------------------------------------

fd.observed$type <- "observed"
fd.mf$type <- "mean field"

fd.null.2 <- fd.null %>% group_by(year,plot,guilds) %>%
  summarise(fd.average.rep = mean(fd.average,na.rm = TRUE),
            fd.sd.rep = sd(fd.average, na.rm = TRUE))
fd.null.2$type <- "null"

# to avoid errors in repeating names
names(fd.null.2)[c(4,5)] <- c("fd.average","fd.sd")

fd.all <- bind_rows(fd.observed,fd.null.2,fd.mf)

fd.all$guilds <- factor(fd.all$guilds,levels = c("plants",
                                                 "floral visitors",
                                                 "herbivores",
                                                 "plants-floral visitors",
                                                 "plants-herbivores","all")) 
fd.all$guilds <- recode(fd.all$guilds, "floral visitors" = "pollinators",
                        "plants-floral visitors" = "plants-pollinators")
fd.all$fd.average[which(is.na(fd.all$fd.average))] <- 0
fd.all$plot <- as.factor(fd.all$plot)
fd.all$type <- as.factor(fd.all$type)

fd.all$type <- recode(fd.all$type, "null" = "randomized")
fd.all$type <- factor(fd.all$type,levels = c("observed","randomized","mean field"))

# plot --------------------------------------------------------------------

library(cowplot)
library(colorblindr)

pd <- .4

fd.means <- fd.all %>% 
  group_by(guilds, type) %>%
  summarise(fd.avg = mean(fd.average),
            fd.sd = sd(fd.average),
            fd.se = sd(fd.average)/sqrt(n()))
  
p3 <- ggplot(fd.all, aes(y = type, x = fd.average)) + 
  geom_point(aes(color = type),
             # shape = 21, 
             position = position_dodge(.5), 
             alpha = .25,
             size = .85) +
  geom_errorbar(data = fd.means, aes(y = type, x = fd.avg,
                                     xmin = fd.avg - fd.se,
                                     xmax = fd.avg + fd.se,
                                     color = type),
                width = 0.15,
                position = position_dodge(pd))+
  geom_point(data = fd.means, aes(y = type, 
                                  x = fd.avg, 
                                  fill = type),
             shape = 21, stroke = .5, size = 2.5,
             position = position_dodge(pd))+
  # geom_boxplot(aes(fill = type),alpha = .6) + 
  scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
  facet_grid(guilds~.)+
  labs(y = "",x = "feasibility domain") +
  # xlim(0,0.26) +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  scale_y_discrete(breaks=NULL, limits=rev) +
  scale_x_continuous(breaks=c(0,0.05,0.1,0.15,0.2,0.25), limits = c(0,0.25)) +
  NULL
# p3


# store plot --------------------------------------------------------------

# ggsave(filename = paste("results/images/Fig_2",vers,".pdf",sep=""),plot = p3,
#        device = cairo_pdf,
#        width = 4,height = 7,dpi = 300)

# temporal variability ----------------------------------------------------
# fd.tmp <- subset(fd.all,guilds == "all" & type == "observed")
fd.tmp <- subset(fd.all,guilds == "all")

fd.tmp$fd.sd <- NULL
fd.wide <- pivot_wider(fd.tmp,names_from = year,names_prefix = "fd.",
                       values_from = fd.average)

p4 <- ggplot(fd.wide,aes(x = fd.2019,y = fd.2020)) + 
  geom_abline(slope = 1,linetype = "dashed",color = "lightgrey") +
  geom_point(aes(fill = type), shape = 21, size = 3) +
  geom_text(aes(label = plot),nudge_y = 0.015) +
  labs(x = "feasibility domain - 2019", 
       y = "feasibility domain - 2020") +
  scale_fill_OkabeIto(darken = 0.2, name = "") +
  # facet_wrap(vars(plot)) +
  facet_grid(.~type) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank()) +
  NULL
# p4

# ggsave(filename = paste("results/images/fd_year_variation",vers,".pdf",sep=""),plot = p4,
#        device = cairo_pdf,
#        width = 8,height = 3,dpi = 300)
