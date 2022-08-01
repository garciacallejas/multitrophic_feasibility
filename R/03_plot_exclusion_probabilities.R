
# plot feasibility results

library(tidyverse)
library(colorblindr)

# read data ---------------------------------------------------------------
# vers <- "_v2"
vers <- ""

load("results/community_names.RData")

# temp
# exc.files <- list.files(path = "results/null_fd",pattern = "null_sp_exclusions_*",full.names = T)
# exc.null <- exc.files %>% map_dfr(read_csv2)
# exc.null[,1] <- NULL
# #
# write.csv2(exc.null,"results/sp_exclusions_null.csv",row.names = F)

file.name <- "results/exclusion_probabilities_observed"
null.name <- "results/sp_exclusions_null_TEMP"
# mean.field.name <- "results/feasibility_domain_mean_field"

file.name <- paste(file.name,vers,".csv",sep="")
null.name <- paste(null.name,vers,".csv",sep="")
# mean.field.name <- paste(mean.field.name,vers,".csv",sep="")

exc.observed <- read.csv2(file = file.name,
                         header = TRUE,
                         stringsAsFactors = FALSE)

exc.null <- read.csv2(file = null.name,
                      header = TRUE,
                      stringsAsFactors = FALSE)
# join data ---------------------------------------------------------------

exc.rank.obs <- exc.observed %>%
  filter(guild == "all") %>%
  group_by(species, guild, intraguild.type) %>%
  summarise(prob_excl_mean = mean(prob_excl_mean)) %>%
  group_by(intraguild.type) %>%
  mutate(type = "observed", 
         prob_excl_rank = rank(-prob_excl_mean,ties.method = "first"))

exc.rank.null <- exc.null %>%
  filter(guild == "all") %>%
  group_by(species, guild, intraguild.type) %>%
  summarise(prob_excl_mean = mean(prob_excl_mean)) %>%
  group_by(intraguild.type) %>%
  mutate(type = "null", 
         prob_excl_rank = rank(-prob_excl_mean,ties.method = "first")) 

# sp - guild - intraguild.type - type - average.prob - rank

exc.rank <- bind_rows(exc.rank.obs, exc.rank.null)

# add guild of each sp
all.plants <- unique(c(sp.names[[1]]$plants,sp.names[[2]]$plants))
all.herb <- unique(c(sp.names[[1]]$herbivores,sp.names[[2]]$herbivores))
all.fv <- unique(c(sp.names[[1]]$floral.visitors,sp.names[[2]]$floral.visitors))

exc.rank$sp.guild <- "unknown"
for(i.obs in 1:nrow(exc.rank)){
  my.sp <- exc.rank$species[i.obs]
  if(my.sp %in% all.herb){
    exc.rank$sp.guild[i.obs] <- "herbivores"
  }else if(my.sp %in% all.fv){
    exc.rank$sp.guild[i.obs] <- "floral visitors"
  }else if(my.sp %in% all.plants){
    exc.rank$sp.guild[i.obs] <- "plants"
  }
}# for i.obs

# -------------------------------------------------------------------------
# plot

exc.rank.plot <- ggplot(exc.rank,aes(x = prob_excl_rank, 
                                     y = prob_excl_mean, 
                                     group = intraguild.type)) + 
  geom_line() + 
  geom_point(aes(color = sp.guild), size = 1.5) +
  facet_wrap(type~intraguild.type, scales = "free")+
  theme_bw() + 
  NULL
exc.rank.plot

exc.box.plot <- ggplot(exc.rank,aes(x = sp.guild, 
                                    y = prob_excl_mean)) + 
  geom_boxplot(aes(fill = sp.guild)) + 
  # geom_point(aes(fill = sp.guild), size = 1.5) +
  facet_wrap(type~intraguild.type, scales = "free")+
  theme_bw() + 
  NULL
exc.box.plot

# store plot --------------------------------------------------------------

ggsave(filename = paste("results/images/boxplot_exclusion_probabilities",vers,".pdf",sep=""),plot = exc.box.plot,
       device = cairo_pdf,
       width = 9,height = 5,dpi = 300)
