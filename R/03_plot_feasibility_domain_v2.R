
# plot feasibility results

library(tidyverse)
library(colorblindr)
library(patchwork)

# read data ---------------------------------------------------------------
# vers <- "_v2"
vers <- ""

# temp
# fd.files <- list.files(path = "results/null_fd",pattern = "null_fd_*",full.names = T)
# fd.null <- fd.files %>% map_dfr(read_csv2)
# fd.null[,1] <- NULL
# 
# write.csv2(fd.null,"results/feasibility_domain_null_TEMP.csv",row.names = F)

file.name <- "results/feasibility_domain_observed"
# null.name <- "results/feasibility_domain_null_TEMP"
# mean.field.name <- "results/feasibility_domain_mean_field"

file.name <- paste(file.name,vers,".csv",sep="")
# null.name <- paste(null.name,vers,".csv",sep="")
# mean.field.name <- paste(mean.field.name,vers,".csv",sep="")

fd <- read.csv2(file = file.name,
                         header = TRUE,
                         stringsAsFactors = FALSE)
names(fd)[names(fd) == "guild"] <- "guilds"

# fd.null <- read.csv2(file = null.name,
#                      header = TRUE,
#                      stringsAsFactors = FALSE)

# join data ---------------------------------------------------------------

# fd.observed$type <- "observed"
# fd.mf$type <- "mean field"
# fd.observed$omega_isotropic_upperCI <- fd.observed$omega_isotropic
# fd.observed$omega_isotropic_lowerCI <- fd.observed$omega_isotropic

# quantile(small_omega_rep, prob=c(.025,.975)) %>% as.numeric()

# fd.null.2 <- fd.null %>% group_by(year,plot,guild,intraguild.type) %>%
#   summarise(
#             omega_lowerCI = as.numeric(quantile(omega_mean,prob=c(.025))),
#             omega_upperCI = as.numeric(quantile(omega_mean,prob=c(.975))),
#             omega_mean = mean(omega_mean,na.rm = TRUE),
#             omega_isotropic_lowerCI = as.numeric(quantile(omega_isotropic,prob=c(.025))),
#             omega_isotropic_upperCI = as.numeric(quantile(omega_isotropic,prob=c(.975))),
#             omega_isotropic = mean(omega_isotropic,na.rm = TRUE),
#             isotropy_index_lowerCI = as.numeric(quantile(isotropy_index_mean,prob=c(.025))),
#             isotropy_index_upperCI = as.numeric(quantile(isotropy_index_mean,prob=c(.975))),
#             isotropy_index_mean = mean(isotropy_index_mean,na.rm = TRUE))
# fd.null.2$type <- "null"

# to avoid errors in repeating names
# names(fd.null.2)[c(4,5)] <- c("fd.average","fd.sd")

# fd.all <- fd.observed
# fd.all <- bind_rows(fd.observed,fd.null.2)

fd$fd.average <- fd$omega_mean
fd$guilds <- factor(fd$guilds,levels = c("plants",
                                         "floral visitors",
                                         "herbivores",
                                         "plants-floral visitors",
                                         "plants-herbivores","all")) 
fd$guilds <- recode(fd$guilds, "plants" = "Plants", 
                    "floral visitors" = "Pollinators",
                    "herbivores" = "Herbivores",
                    "plants-floral visitors" = "Plants-Pollinators",
                    "plants-herbivores" = "Plants-Herbivores",
                    "all" = "All")

fd$fd.average[which(is.na(fd$fd.average))] <- 0
fd$plot <- as.factor(fd$plot)
# fd$type <- as.factor(fd$type)

# fd$type <- recode(fd$type, "null" = "Randomized", "observed" = "Observed")
# fd$type <- factor(fd$type,levels = c("Observed","Randomized"))
fd$intraguild.type[fd$intraguild.type == "mean_field"] <- "Mean field"
fd$intraguild.type[fd$intraguild.type == "nesting_larvae"] <- "Resources"
fd$intraguild.type[fd$intraguild.type == "nesting_larvae_phenology"] <- "Resources and \nphenology"
# plot --------------------------------------------------------------------

# fd.all.nona <- subset(fd,!is.na(fd.average))

fd.means <- fd %>% 
  # group_by(guilds, type) %>%
  group_by(intraguild.type, guilds) %>%
  summarise(fd.avg = mean(fd.average,na.rm = TRUE),
            fd.sd = sd(fd.average,na.rm = TRUE),
            fd.se = sd(fd.average,na.rm = TRUE)/sqrt(n()),
            iso.fd.avg = mean(omega_isotropic,na.rm = TRUE),
            iso.fd.sd = sd(omega_isotropic,na.rm = TRUE),
            iso.fd.se = sd(omega_isotropic,na.rm = TRUE)/sqrt(n()),
            iso.index.avg = mean(isotropy_index_mean,na.rm = TRUE),
            iso.index.sd = sd(isotropy_index_mean,na.rm = TRUE),
            iso.index.se = sd(isotropy_index_mean,na.rm = TRUE)/sqrt(n()))

pj <- .21

standard.fd <- ggplot(fd, aes(x = intraguild.type, y = fd.average)) + 
  geom_point(aes(color = intraguild.type),
             # shape = 21, 
             position=position_jitter(width = pj),
             # size = .85, 
             alpha = .7) +
  geom_boxplot(aes(fill = intraguild.type),
               outlier.shape = NA,
               alpha = .4) +
  # geom_pointrange(data = fd.means,aes(x = intraguild.type, y = fd.avg,
  #                                     ymin = fd.avg - fd.sd,
  #                                     ymax = fd.avg + fd.sd,
  #                                     fill = intraguild.type,
  #                                     color = intraguild.type,
  #                                     group = intraguild.type),
  #                 # position = position_jitterdodge(jitter.width = pj),
  #                 shape = 21)+
                    
  # geom_errorbar(data = fd.means, aes(x = guilds, y = fd.avg,
  #                                    ymin = fd.avg - fd.sd,
  #                                    ymax = fd.avg + fd.sd,
  #                                    color = intraguild.type, group = intraguild.type),
  #               position = position_dodge(width = .5),
  #               
  #               width = 0.15)+
  #               # position = position_jitterdodge())+
  # 
  # geom_point(data = fd.means, aes(x = guilds,
  #                                 y = fd.avg,
  #                                 fill = intraguild.type),
  #            shape = 21,
  #            stroke = .5, 
  #            size = 2.5,
  #            position = position_jitterdodge(jitter.width = .11))+
  # geom_boxplot(aes(fill = type),alpha = .6) + 
  scale_fill_OkabeIto(darken = 0.2, name = "Intraguild \ninteractions") +
  scale_color_OkabeIto(darken = 0.2, name = "Intraguild \ninteractions") +
  # scale_shape_manual(values = c(21,23)) +
  facet_grid(.~guilds, drop = T)+
  labs(x = "",y = "Feasibility domain") +
  # xlim(0,0.26) +
  theme_bw()+
  # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
  theme(strip.background = element_blank())+
  # scale_y_discrete(breaks=NULL, limits=rev) +
  # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
  # guides(color=guide_legend(override.aes = list(shape=21))) + 
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  NULL
standard.fd

ggsave(filename = paste("results/images/Fig_fd.pdf",sep=""),plot = standard.fd,
       device = cairo_pdf,
       width = 9,height = 4,dpi = 300)

# -------------------------------------------------------------------------
# 
# iso.fd <- ggplot(fd.all.nona, aes(x = intraguild.type, y = omega_isotropic)) + 
# 
#   geom_errorbar(data = fd.means, aes(x = intraguild.type, y = iso.fd.avg,
#                                      ymin = iso.fd.avg - iso.fd.sd,
#                                      ymax = iso.fd.avg + iso.fd.sd,
#                                      color = intraguild.type,
#                                      linetype = type),
#                 width = 0.15,
#                 position = position_dodge(pd))+
#   geom_point(data = fd.means, aes(x = intraguild.type,
#                                   y = iso.fd.avg,
#                                   fill = intraguild.type),
#              shape = 21,
#              stroke = .5, 
#              size = 2.5,
#              position = position_dodge(pd))+
#   geom_point(aes(color = intraguild.type, shape = type),
#              # shape = 21,
#              position = position_dodge(pd),
#              alpha = .25,
#              size = .85) +
#   # geom_boxplot(aes(fill = type),alpha = .6) + 
#   scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
#   scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
#   scale_shape_manual(values = c(21,23)) +
#   facet_grid(.~guild)+
#   labs(x = "",y = "Isotropic feasibility domain") +
#   # xlim(0,0.26) +
#   theme_bw()+
#   # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
#   theme(strip.background = element_blank())+
#   # scale_y_discrete(breaks=NULL, limits=rev) +
#   # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
#   guides(color=guide_legend(override.aes = list(shape=21))) + 
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
#   NULL
# iso.fd
# 
# # -------------------------------------------------------------------------
# 
# iso.index <- ggplot(fd.all.nona, aes(x = intraguild.type, y = isotropy_index_mean)) + 
#   geom_point(aes(color = intraguild.type, shape = type),
#              # shape = 21, 
#              position = position_dodge(pd), 
#              alpha = .25,
#              size = .85) +
#   geom_errorbar(data = fd.means, aes(x = intraguild.type, y = iso.index.avg,
#                                      ymin = iso.index.avg - iso.index.sd,
#                                      ymax = iso.index.avg + iso.index.sd,
#                                      color = intraguild.type,
#                                      linetype = type),
#                 width = 0.15,
#                 position = position_dodge(pd))+
#   geom_point(data = fd.means, aes(x = intraguild.type,
#                                   y = iso.index.avg,
#                                   fill = intraguild.type,
#                                   shape = type),
#              # shape = 21, 
#              stroke = .5, 
#              size = 2.5,
#              position = position_dodge(pd))+
#   # geom_boxplot(aes(fill = type),alpha = .6) + 
#   scale_fill_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
#   scale_color_OkabeIto(darken = 0.2, name = "Interaction\nstructure") +
#   scale_shape_manual(values = c(21,23)) +
#   facet_grid(.~guild)+
#   labs(x = "",y = "Isotropy index of feasibility domain") +
#   # xlim(0,0.26) +
#   theme_bw()+
#   # theme(axis.text.x  = element_text(angle = 60, hjust = 1))+
#   theme(strip.background = element_blank())+
#   # scale_y_discrete(breaks=NULL, limits=rev) +
#   # scale_x_continuous(breaks=seq(0,0.5,by = 0.05), limits = c(0,0.5)) +
#   guides(color=guide_legend(override.aes = list(shape=21))) + 
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
#   NULL
# iso.index
# 
# # store plot --------------------------------------------------------------
# 
# ggsave(filename = paste("results/images/standard_fd",vers,".pdf",sep=""),plot = standard.fd,
#        device = cairo_pdf,
#        width = 10,height = 3,dpi = 300)
# ggsave(filename = paste("results/images/isotropic_fd",vers,".pdf",sep=""),plot = iso.fd,
#        device = cairo_pdf,
#        width = 10,height = 3,dpi = 300)
# ggsave(filename = paste("results/images/isotropy_index",vers,".pdf",sep=""),plot = iso.index,
#        device = cairo_pdf,
#        width = 10,height = 3,dpi = 300)
# 
# full.plot <- standard.fd/iso.fd/iso.index
# ggsave(filename = paste("results/images/all_fd",vers,".pdf",sep=""),plot = full.plot,
#        device = cairo_pdf,
#        width = 10,height = 9,dpi = 300)
