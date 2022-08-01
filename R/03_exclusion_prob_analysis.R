
library(tidyverse)
library(DHARMa)
library(broom.mixed)

# -------------------------------------------------------------------------

exc.obs <- read.csv2("results/exclusion_probabilities_observed.csv")
names(exc.obs)[names(exc.obs) == "guild"] <- "guilds"
sp.metrics <- read.csv2("results/species_level_metrics.csv")

# -------------------------------------------------------------------------
# subset to be consistent with the community-level analysis

sp.obs <- left_join(exc.obs,sp.metrics) %>%
  filter(intraguild.type == "nesting_larvae_phenology" &
           guilds == "all") %>%
  select(species,sp.guild,year,plot,prob_excl_mean,diagonal_dominance,in_degree,out_degree)

# -------------------------------------------------------------------------
# visualization
# hist(sp.obs$diagonal_dominance)
# hist(sp.obs$in_degree)
# hist(sp.obs$out_degree)
# plot(sp.obs$diagonal_dominance,sp.obs$prob_excl_mean)
# plot(sp.obs$diagonal_dominance,log(sp.obs$prob_excl_mean))

# very skewed response
hist(sp.obs$prob_excl_mean)
sum(sp.obs$prob_excl_mean == 0)

ms <- lmerTest::lmer(prob_excl_mean ~ sp.guild + diagonal_dominance + in_degree + (1|plot),data = sp.obs)
# summary(ms)
# # residuals do not look good
# DHARMa::testResiduals(ms)

# try transforming the data - somewhat better
ms2 <- lmerTest::lmer(log(prob_excl_mean) ~ sp.guild + diagonal_dominance + in_degree + (1|plot),data = sp.obs)
# ms2 <- lmerTest::lmer(log(prob_excl_mean) ~ sp.guild + diagonal_dominance + in_degree + (1|plot),data = subset(sp.obs,diagonal_dominance < 0.2))
summary(ms2)
DHARMa::testResiduals(ms2)
broom.mixed::tidy(ms2)
