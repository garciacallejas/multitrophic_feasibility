
library(tidyverse)
library(ineq) # for the gini index
library(robustbase)
library(igraph)
# devtools::install_github("RadicalCommEcol/MultitrophicFun",
# auth_token = "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")
library(MultitrophicFun)

source("R/aux_linkRankModularity.R")
source("R/aux_allowed_interactions.R")

# see https://ecological-complexity-lab.github.io/infomap_ecology_package/
library(infomapecology)
source("R/aux_run_infomap_monolayer2.R")
check_infomap()

# include null matrices ---------------------------------------------------

include.null <- TRUE
include.mean.field <- TRUE

# read community matrices -------------------------------------------------

# vers <- "_v2"
vers <- ""

# -------------------------------------------------------------------------

load(file = paste("results/community_matrices_observed",vers,".Rdata",sep=""))

# load sp names as well
load(file = paste("results/community_names",vers,".RData",sep=""))
plants <- sort(unique(c(sp.names[["2019"]][["plants"]],sp.names[["2020"]][["plants"]])))
fv <- sort(unique(c(sp.names[["2019"]][["floral.visitors"]],sp.names[["2020"]][["floral.visitors"]])))
herb <- sort(unique(c(sp.names[["2019"]][["herbivores"]],sp.names[["2020"]][["herbivores"]])))

# read structural zeros ---------------------------------------------------

sz <- read.csv2(paste("data/potential_interactions",vers,".csv",sep=""))

# sz are divided by year, so join them together in a single set
# such that pairs that are structural zero all years are kept, and
# those that have at least one observation in one year are not structural zeros

sz.combined <- sz %>% pivot_wider(names_from = year,values_from = c(structural.zero,observed.freq))

sz.combined$structural.zero_2019[which(is.na(sz.combined$structural.zero_2019))] <- TRUE
sz.combined$structural.zero_2020[which(is.na(sz.combined$structural.zero_2020))] <- TRUE
sz.combined$structural.zero <- sz.combined$structural.zero_2019 * sz.combined$structural.zero_2020

sz.combined$observed.freq_2019[which(is.na(sz.combined$observed.freq_2019))] <- 0
sz.combined$observed.freq_2020[which(is.na(sz.combined$observed.freq_2020))] <- 0
sz.combined$observed.freq <- sz.combined$observed.freq_2019 + sz.combined$observed.freq_2020
sz.combined$structural.zero[sz.combined$observed.freq > 0] <- 0

sz.combined <- sz.combined[,c("plant","animal","type","overlap","observed.freq","structural.zero")]
sz.combined$structural.zero <- as.logical(sz.combined$structural.zero)

# define range ------------------------------------------------------------

years <- c(2019,2020)
plots <- 1:9

guild.combinations <- c("plants",
                        "floral visitors",
                        "herbivores",
                        "plants-floral visitors",
                        "plants-herbivores",
                        "all")
# define metrics ----------------------------------------------------------

# constants for the modularity algorithm
lr.damping <- .85
lr.alg <- "prpack"

# minimum number of nodes to calculate modularity
min.node.modularity <- 2

metric.names <- c("richness",
                  "connectance",
                  "intra_inter_ratio",
                  "degree_distribution",
                  # "kurtosis",
                  "skewness",
                  "modularity",
                  "complexity")

community.metrics <- expand.grid(year = years, 
                                 plot = plots, 
                                 guilds = guild.combinations,
                                 metric = metric.names,value = NA)

module_members <- NULL

# calculate metrics -------------------------------------------------------

for(i.year in 1:length(years)){
  for(i.plot in plots){
    
    year.plot.matrix <- 
      community_matrices[[as.character(years[i.year])]][[i.plot]]
    
    plant.positions <- which(rownames(year.plot.matrix) %in% 
                               sp.names[[i.year]][["plants"]])
    herb.positions <- which(rownames(year.plot.matrix) %in% 
                              sp.names[[i.year]][["herbivores"]])
    fv.positions <- which(rownames(year.plot.matrix) %in% 
                            sp.names[[i.year]][["floral.visitors"]])
    
    for(i.guild in 1:length(guild.combinations)){
      
      cat(years[i.year],"-",i.plot,"-",guild.combinations[i.guild],"\n")
      
      if(guild.combinations[i.guild] == "plants"){
        
        my.matrix <- year.plot.matrix[plant.positions,plant.positions]
        
      }else if(guild.combinations[i.guild] == "plants-floral visitors"){
        
        my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
                                      c(plant.positions,fv.positions)]
        
      }else if(guild.combinations[i.guild] == "floral visitors"){
        
        my.matrix <- year.plot.matrix[fv.positions,fv.positions]
        
      }else if(guild.combinations[i.guild] == "herbivores"){
        
        my.matrix <- year.plot.matrix[herb.positions,herb.positions]
        
      }else if(guild.combinations[i.guild] == "plants-herbivores"){
        
        my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
                                      c(plant.positions,herb.positions)]
        
      }else if(guild.combinations[i.guild] == "all"){
        
        my.matrix <- year.plot.matrix
      }
      
      # convert to igraph object ------------------------------------------------
      
      my.network <- graph.adjacency(my.matrix,
                                    mode="directed",
                                    weighted=TRUE)

      # richness ----------------------------------------------------------------

      my.richness <- length(unique(c(rownames(my.matrix),colnames(my.matrix))))

      # intra/inter ratio -------------------------------------------------------

      my.ratio <- sum(diag(my.matrix))/sum(my.matrix[col(my.matrix) != row(my.matrix)])

      # connectance -------------------------------------------------------------
      # for now, binary
      
      sz.all <- get_structural_zeros(A = my.matrix,
                                 plants = plants, 
                                 fv = fv, 
                                 herb = herb,
                                 dfz = sz.combined)

      my.connectance <- connectance(interaction.matrix = my.matrix,
                                    quant = FALSE,
                                    structural_zeros = sz.all)

      # degree distribution -----------------------------------------------------
      # first take: gini index
      my.degrees <- igraph::degree(my.network)
      my.deg.dist <- ineq::ineq(my.degrees,type = "Gini")

      # kurtosis ----------------------------------------------------------------

      # my.kurtosis <- moments::kurtosis(as.numeric(my.matrix))
      my.skew <- robustbase::mc(as.numeric(my.matrix))

      # May's complexity --------------------------------------------------------
      
      # connectance explicitly deals with structural zeros
      # so sigma must, as well.
      
      all.interactions <- expand.grid(row = 1:nrow(my.matrix),
                                      col = 1:ncol(my.matrix))
      all.interactions$value <- c(my.matrix)
      
      # return all interactions wihout a match in sz
      valid.interactions <- anti_join(all.interactions,sz.all)
      
      # sd of valid interaction strengths
      my.sigma <- sd(valid.interactions$value)
      
      my.complexity <- MultitrophicFun::May_complexity(S = nrow(my.matrix),
                                                       C = my.connectance,
                                                       sigma = my.sigma)
      
      # modularity --------------------------------------------------------------
      
      if(nrow(my.matrix) >= min.node.modularity){
      
      g  <- graph.adjacency(my.matrix,weighted=TRUE)
      
      # To run infomap it seems that weights should be non-negative
      my.edge.list <- get.data.frame(g) %>% mutate(weight=abs(weight)) 
      
      nodes <- my.edge.list$from %>% unique()
      
      nodes.ID <- tibble(node_id=as.numeric(1:length(nodes)),species=nodes)
      
      # Preparing edge.list and nodes.ID to run infomap
      
      my.edge.list.ID <- my.edge.list %>% rename(species=from) %>%
        left_join(nodes.ID,by="species") %>% select(-species) %>% 
        rename(from=node_id,species=to) %>%
        left_join(nodes.ID,by="species") %>% 
        select(-species) %>% rename(to=node_id) %>% select(from,to,weight)
      
      nodes.ID2 <- nodes.ID %>% rename(node_name=species)
      
      infomap_mono <- create_monolayer_object(x = my.edge.list.ID, 
                                              directed = T, 
                                              bipartite = F, 
                                              node_metadata = nodes.ID2)
      
      infomap_mono2 <- infomap_mono
      infomap_mono2$nodes$node_name <- as.numeric(infomap_mono2$nodes$node_name)
      
      # Run Infomap
      
      modules_relax_rate <- run_infomap_monolayer2(x = infomap_mono2, 
                                                   flow_model = 'directed', 
                                                   silent=T,
                                                   trials=1000, 
                                                   two_level=T, 
                                                   seed=200952)
      
      my.modularity <- modules_relax_rate[["L"]] # modularity in bits
      
      # Extract module information
      modules <- modules_relax_rate$modules %>% 
        dplyr::select(node_id,module_level1) %>% 
        rename(module=module_level1) %>%
        left_join(nodes.ID,by="node_id")
      
      
      modules.aux <- tibble(year = years[i.year],
                            plot = plots[i.plot],
                            guild = guild.combinations[i.guild],
                            species = modules$species,
                            module = modules$module) 
      
      module_members <- bind_rows(module_members,modules.aux)
      
      # If we dont want to use bits, I guess that we could translate the previous partition
      # to other units by uning the alternative definition of modularity
      # E. A. Leicht and M. E. J. Newman, Phys. Rev. Lett. 100, 118703 (2008).
      
      # Note 2: infomap do not optimise this generalized modularity function. Thus,
      # I dont expect high values.
      
      # According to results for plot 1, year 2019 and guild == "floral visitors",
      # it seems that it seems that isolated nodes module is NA (see species "Apoidea")
      # For each one of those nodes, we will add a new module. We denote such partition
      # as "corrected partition"
      
      # CORRECTED partition from infomap
      module_max <- max(modules$module,na.rm = T)
      
      for(i in 1:nrow(modules)){
        
        if(is.na(modules$module[i])){
          module_max <- module_max + 1
          modules$module[i] <- module_max
        }
      }
      
      # linkrank modularity -----------------------------------------------------
      # Note 1: I will use non-negative weights to be consistent with 
      # the inputs that we used to feed the infomap algorithm.
      
      g  <- graph.adjacency(abs(my.matrix),weighted=TRUE)
      
      # this returns a "mask", a matrix of the same dimensions as my.matrix
      # with zeros if an interaction is forbidden, 1 if it is allowed.
      allowed.interactions.mat <- allowed.interactions(my.matrix, year.plot.matrix,
                                                       guild.combinations[i.guild],
                                                       sp.names, i.year)
      
      linkrank_modularity <- linkRankModularity(g,
                                                partition=modules$module, 
                                                allowed.interactions.mat = allowed.interactions.mat,
                                                damping = lr.damping, 
                                                pr.algo = lr.alg)
      }else{ # less than 5 sp
        linkrank_modularity <- NA_real_
      }
      # store -------------------------------------------------------------------
      
      pos <- which(community.metrics$year == years[i.year] &
                     community.metrics$plot == i.plot &
                     community.metrics$guilds == guild.combinations[i.guild])
      
      community.metrics$value[pos[which(community.metrics$metric[pos] == "richness")]] <- my.richness
      community.metrics$value[pos[which(community.metrics$metric[pos] == "connectance")]] <- my.connectance
      community.metrics$value[pos[which(community.metrics$metric[pos] == "intra_inter_ratio")]] <- my.ratio
      community.metrics$value[pos[which(community.metrics$metric[pos] == "degree_distribution")]] <- my.deg.dist
      community.metrics$value[pos[which(community.metrics$metric[pos] == "complexity")]] <- my.complexity
      community.metrics$value[pos[which(community.metrics$metric[pos] == "modularity")]] <- linkrank_modularity
      # community.metrics$value[pos[which(community.metrics$metric[pos] == "kurtosis")]] <- my.kurtosis
      community.metrics$value[pos[which(community.metrics$metric[pos] == "skewness")]] <- my.skew
      
    }# for i.guild
  }# for i.plot
}# for i.year

# store metrics -----------------------------------------------------------

# year - plot - guilds - metric - value
# store a different file for observed and null networks

write.csv2(community.metrics,file = paste("results/community_metrics_observed",vers,".csv",sep=""),
           row.names = FALSE)
# write.csv2(module_members,file = "results/infomap_modular_partitions.csv",
#            row.names = FALSE)

# repeat for null matrices ------------------------------------------------

if(include.null){
  
  # load data ---------------------------------------------------------------
  
  load(paste("results/community_matrices_null",vers,".Rdata",sep=""))
  
  replicates <- 100#length(null_matrices)
  
  # results dataframe -------------------------------------------------------
  
  community.null.metrics <- expand.grid(year = years, 
                                        plot = plots, 
                                        replicate = 1:replicates,
                                        guilds = guild.combinations,
                                        metric = metric.names,value = NA)
  
  for(i.rep in 1:replicates){
    
    # calculate metrics -------------------------------------------------------
    
    for(i.year in 1:length(years)){
      for(i.plot in plots){
        
        year.plot.matrix <- 
          null_matrices[[i.rep]][[as.character(years[i.year])]][[i.plot]]
        
        plant.positions <- which(rownames(year.plot.matrix) %in% 
                                   sp.names[[i.year]][["plants"]])
        herb.positions <- which(rownames(year.plot.matrix) %in% 
                                  sp.names[[i.year]][["herbivores"]])
        fv.positions <- which(rownames(year.plot.matrix) %in% 
                                sp.names[[i.year]][["floral.visitors"]])
        
        for(i.guild in 1:length(guild.combinations)){
          
          cat(i.rep,"-",years[i.year],"-",i.plot,"-",
              guild.combinations[i.guild],"\n")
          
          if(guild.combinations[i.guild] == "plants"){
            
            my.matrix <- year.plot.matrix[plant.positions,plant.positions]
            
          }else if(guild.combinations[i.guild] == "plants-floral visitors"){
            
            my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
                                          c(plant.positions,fv.positions)]
            
          }else if(guild.combinations[i.guild] == "floral visitors"){
            
            my.matrix <- year.plot.matrix[fv.positions,fv.positions]
            
          }else if(guild.combinations[i.guild] == "herbivores"){
            
            my.matrix <- year.plot.matrix[herb.positions,herb.positions]
            
          }else if(guild.combinations[i.guild] == "plants-herbivores"){
            
            my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
                                          c(plant.positions,herb.positions)]
            
          }else if(guild.combinations[i.guild] == "all"){
            
            my.matrix <- year.plot.matrix
          }
          
          # convert to igraph object ------------------------------------------------
          
          my.network <- graph.adjacency(my.matrix, 
                                        mode="directed", 
                                        weighted=TRUE)
          # richness ----------------------------------------------------------------
          
          my.richness <- length(unique(c(rownames(my.matrix),colnames(my.matrix))))
          
          # intra/inter ratio -------------------------------------------------------
          
          my.ratio <- sum(diag(my.matrix))/sum(my.matrix[col(my.matrix) != row(my.matrix)])
          
          # connectance -------------------------------------------------------------
          # for now, binary
          
          sz <- get_structural_zeros(A = my.matrix,
                                     plants = plants, 
                                     fv = fv, 
                                     herb = herb)
          
          my.connectance <- connectance(interaction.matrix = my.matrix,
                                        quant = FALSE,
                                        structural_zeros = sz)
          
          # degree distribution -----------------------------------------------------
          # first take: gini index
          my.degrees <- igraph::degree(my.network)
          my.deg.dist <- ineq::ineq(my.degrees,type = "Gini")
          
          # kurtosis ----------------------------------------------------------------
          
          # my.kurtosis <- moments::kurtosis(as.numeric(my.matrix))      
          my.skew <- robustbase::mc(as.numeric(my.matrix))
          
          # May's complexity --------------------------------------------------------
          
          # connectance explicitly deals with structural zeros
          # so sigma must, as well.
          
          all.interactions <- expand.grid(row = 1:nrow(my.matrix),
                                          col = 1:ncol(my.matrix))
          all.interactions$value <- c(my.matrix)
          
          # return all interactions wihout a match in sz
          valid.interactions <- anti_join(all.interactions,sz)
          
          # sd of valid interaction strengths
          my.sigma <- sd(valid.interactions$value)
          
          my.complexity <- MultitrophicFun::May_complexity(S = nrow(my.matrix),
                                                           C = my.connectance,
                                                           sigma = my.sigma)
          
          # modularity --------------------------------------------------------------
          
          if(nrow(my.matrix) >= min.node.modularity){
          
          # we compute the infomap partition and extract the associated Q
          # including any potential isolated components
          # the methodology comes from
          # https://doi.org/10.1103/PhysRevLett.100.118703
          
          g  <- graph.adjacency(my.matrix,weighted=TRUE)
          
          # To run infomap it seems that weights should be non-negative
          my.edge.list <- get.data.frame(g) %>% mutate(weight=abs(weight)) 
          
          nodes <- my.edge.list$from %>% unique()
          
          nodes.ID <- tibble(node_id=as.numeric(1:length(nodes)),species=nodes)
          
          # Preparing edge.list and nodes.ID to run infomap
          
          my.edge.list.ID <- my.edge.list %>% rename(species=from) %>%
            left_join(nodes.ID,by="species") %>% select(-species) %>% 
            rename(from=node_id,species=to) %>%
            left_join(nodes.ID,by="species") %>% 
            select(-species) %>% rename(to=node_id) %>% select(from,to,weight)
          
          nodes.ID2 <- nodes.ID %>% rename(node_name=species)
          
          infomap_mono <- create_monolayer_object(x = my.edge.list.ID, 
                                                  directed = T, 
                                                  bipartite = F, 
                                                  node_metadata = nodes.ID2)
          
          infomap_mono2 <- infomap_mono
          infomap_mono2$nodes$node_name <- as.numeric(infomap_mono2$nodes$node_name)
          
          # Run Infomap
          
          modules_relax_rate <- run_infomap_monolayer2(x = infomap_mono2, 
                                                       flow_model = 'directed', 
                                                       silent=T,
                                                       trials=1000, 
                                                       two_level=T, 
                                                       seed=200952)
          
          my.modularity <- modules_relax_rate[["L"]] # modularity in bits
          
          # Extract module information
          modules <- modules_relax_rate$modules %>% 
            dplyr::select(node_id,module_level1) %>% 
            rename(module=module_level1) %>%
            left_join(nodes.ID,by="node_id")
          
          
          modules.aux <- tibble(year = years[i.year],
                                plot = plots[i.plot],
                                guild = guild.combinations[i.guild],
                                species = modules$species,
                                module = modules$module) 
          
          module_members <- bind_rows(module_members,modules.aux)
          
          # If we dont want to use bits, I guess that we could translate the previous partition
          # to other units by uning the alternative definition of modularity
          # E. A. Leicht and M. E. J. Newman, Phys. Rev. Lett. 100, 118703 (2008).
          
          # Note 2: infomap do not optimise this generalized modularity function. Thus,
          # I dont expect high values.
          
          # According to results for plot 1, year 2019 and guild == "floral visitors",
          # it seems that it seems that isolated nodes module is NA (see species "Apoidea")
          # For each one of those nodes, we will add a new module. We denote such partition
          # as "corrected partition"
          
          # CORRECTED partition from infomap
          module_max <- max(modules$module,na.rm = T)
          
          for(i in 1:nrow(modules)){
            
            if(is.na(modules$module[i])){
              module_max <- module_max + 1
              modules$module[i] <- module_max
            }
          }
          
          # linkrank modularity -----------------------------------------------------
          # Note 1: I will use non-negative weights to be consistent with 
          # the inputs that we used to feed the infomap algorithm.
          
          g  <- graph.adjacency(abs(my.matrix),weighted=TRUE)
          
          allowed.interactions.mat <- allowed.interactions(my.matrix, year.plot.matrix,
                                                           guild.combinations[i.guild],
                                                           sp.names, i.year)
          
          linkrank_modularity <- linkRankModularity(g,
                                                    partition=modules$module, 
                                                    allowed.interactions.mat = allowed.interactions.mat,
                                                    damping = lr.damping, 
                                                    pr.algo = lr.alg)
          }else{ # less than 5 sp
            linkrank_modularity <- NA_real_
          }
          # store -------------------------------------------------------------------
          
          pos <- which(community.null.metrics$year == years[i.year] &
                         community.null.metrics$plot == i.plot &
                         community.null.metrics$replicate == i.rep &
                         community.null.metrics$guilds == guild.combinations[i.guild])
          
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "richness")]] <- my.richness
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "connectance")]] <- my.connectance
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "intra_inter_ratio")]] <- my.ratio
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "degree_distribution")]] <- my.deg.dist
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "modularity")]] <- linkrank_modularity
          # community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "kurtosis")]] <- my.kurtosis
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "skewness")]] <- my.skew
          community.null.metrics$value[pos[which(community.null.metrics$metric[pos] == "complexity")]] <- my.complexity
          
        }# for i.guild
      }# for i.plot
    }# for i.year
  }# for i.rep
  
  # store null --------------------------------------------------------------
  
  write.csv2(community.null.metrics,file = paste("results/community_metrics_null",vers,".csv",sep=""))
  
}# if include.null

# mean field matrices -----------------------------------------------------

if(include.mean.field){
  
  # load data ---------------------------------------------------------------
  
  load(paste("results/community_matrices_mean_field",vers,".Rdata",sep=""))
  
  # results dataframe -------------------------------------------------------
  
  community.mean.field.metrics <- expand.grid(year = years, 
                                        plot = plots, 
                                        guilds = guild.combinations,
                                        metric = metric.names,value = NA)
  
  for(i.year in 1:length(years)){
    for(i.plot in plots){
      
      year.plot.matrix <- 
        mean_field_matrices[[as.character(years[i.year])]][[i.plot]]
      
      plant.positions <- which(rownames(year.plot.matrix) %in% 
                                 sp.names[[i.year]][["plants"]])
      herb.positions <- which(rownames(year.plot.matrix) %in% 
                                sp.names[[i.year]][["herbivores"]])
      fv.positions <- which(rownames(year.plot.matrix) %in% 
                              sp.names[[i.year]][["floral.visitors"]])
      
      for(i.guild in 1:length(guild.combinations)){
        
        cat(years[i.year],"-",i.plot,"-",guild.combinations[i.guild],"\n")
        
        if(guild.combinations[i.guild] == "plants"){
          
          my.matrix <- year.plot.matrix[plant.positions,plant.positions]
          
        }else if(guild.combinations[i.guild] == "plants-floral visitors"){
          
          my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
                                        c(plant.positions,fv.positions)]
          
        }else if(guild.combinations[i.guild] == "floral visitors"){
          
          my.matrix <- year.plot.matrix[fv.positions,fv.positions]
          
        }else if(guild.combinations[i.guild] == "herbivores"){
          
          my.matrix <- year.plot.matrix[herb.positions,herb.positions]
          
        }else if(guild.combinations[i.guild] == "plants-herbivores"){
          
          my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
                                        c(plant.positions,herb.positions)]
          
        }else if(guild.combinations[i.guild] == "all"){
          
          my.matrix <- year.plot.matrix
        }
        
        
        # convert to igraph object ------------------------------------------------
        
        my.network <- graph.adjacency(my.matrix, 
                                      mode="directed", 
                                      weighted=TRUE)
        
        # richness ----------------------------------------------------------------
        
        my.richness <- length(unique(c(rownames(my.matrix),colnames(my.matrix))))
        
        # intra/inter ratio -------------------------------------------------------
        
        # TODO: check again this. Use intra/sum(inter)
        my.ratio <- sum(diag(my.matrix))/sum(my.matrix[col(my.matrix) != row(my.matrix)])
        
        # connectance -------------------------------------------------------------
        # for now, binary
        
        sz <- get_structural_zeros(A = my.matrix,
                                   plants = plants, 
                                   fv = fv, 
                                   herb = herb)
        
        my.connectance <- connectance(interaction.matrix = my.matrix,
                                      quant = FALSE,
                                      structural_zeros = sz)
        
        # degree distribution -----------------------------------------------------
        # first take: gini index
        my.degrees <- igraph::degree(my.network)
        my.deg.dist <- ineq::ineq(my.degrees,type = "Gini")
        
        # kurtosis ----------------------------------------------------------------
        
        # my.kurtosis <- moments::kurtosis(as.numeric(my.matrix))      
        my.skew <- robustbase::mc(as.numeric(my.matrix))
        
        # May's complexity --------------------------------------------------------
        
        # connectance explicitly deals with structural zeros
        # so sigma must, as well.
        
        all.interactions <- expand.grid(row = 1:nrow(my.matrix),
                                        col = 1:ncol(my.matrix))
        all.interactions$value <- c(my.matrix)
        
        # return all interactions wihout a match in sz
        valid.interactions <- anti_join(all.interactions,sz)
        
        # sd of valid interaction strengths
        my.sigma <- sd(valid.interactions$value)
        
        my.complexity <- MultitrophicFun::May_complexity(S = nrow(my.matrix),
                                                         C = my.connectance,
                                                         sigma = my.sigma)
        
        # modularity --------------------------------------------------------------
        
        if(nrow(my.matrix) >= min.node.modularity){
          
        # we compute the infomap partition and extract the associated Q
        # including any potential isolated components
        # the methodology comes from
        # https://doi.org/10.1103/PhysRevLett.100.118703
        
        g  <- graph.adjacency(my.matrix,weighted=TRUE)
        
        # To run infomap it seems that weights should be non-negative
        my.edge.list <- get.data.frame(g) %>% mutate(weight=abs(weight)) 
        
        nodes <- my.edge.list$from %>% unique()
        
        nodes.ID <- tibble(node_id=as.numeric(1:length(nodes)),species=nodes)
        
        # Preparing edge.list and nodes.ID to run infomap
        
        my.edge.list.ID <- my.edge.list %>% rename(species=from) %>%
          left_join(nodes.ID,by="species") %>% select(-species) %>% 
          rename(from=node_id,species=to) %>%
          left_join(nodes.ID,by="species") %>% 
          select(-species) %>% rename(to=node_id) %>% select(from,to,weight)
        
        nodes.ID2 <- nodes.ID %>% rename(node_name=species)
        
        infomap_mono <- create_monolayer_object(x = my.edge.list.ID, 
                                                directed = T, 
                                                bipartite = F, 
                                                node_metadata = nodes.ID2)
        
        infomap_mono2 <- infomap_mono
        infomap_mono2$nodes$node_name <- as.numeric(infomap_mono2$nodes$node_name)
        
        # Run Infomap
        
        modules_relax_rate <- run_infomap_monolayer2(x = infomap_mono2, 
                                                     flow_model = 'directed', 
                                                     silent=T,
                                                     trials=1000, 
                                                     two_level=T, 
                                                     seed=200952)
        
        my.modularity <- modules_relax_rate[["L"]] # modularity in bits
        
        # Extract module information
        modules <- modules_relax_rate$modules %>% 
          dplyr::select(node_id,module_level1) %>% 
          rename(module=module_level1) %>%
          left_join(nodes.ID,by="node_id")
        
        
        modules.aux <- tibble(year = years[i.year],
                              plot = plots[i.plot],
                              guild = guild.combinations[i.guild],
                              species = modules$species,
                              module = modules$module) 
        
        module_members <- bind_rows(module_members,modules.aux)
        
        # If we dont want to use bits, I guess that we could translate the previous partition
        # to other units by uning the alternative definition of modularity
        # E. A. Leicht and M. E. J. Newman, Phys. Rev. Lett. 100, 118703 (2008).

        # Note 2: infomap do not optimise this generalized modularity function. Thus,
        # I dont expect high values.
        
        # According to results for plot 1, year 2019 and guild == "floral visitors",
        # it seems that it seems that isolated nodes module is NA (see species "Apoidea")
        # For each one of those nodes, we will add a new module. We denote such partition
        # as "corrected partition"
        
        # CORRECTED partition from infomap
        module_max <- max(modules$module,na.rm = T)
        
        for(i in 1:nrow(modules)){
          
          if(is.na(modules$module[i])){
            module_max <- module_max + 1
            modules$module[i] <- module_max
          }
        }
        
        # linkrank modularity -----------------------------------------------------
        # Note 1: I will use non-negative weights to be consistent with 
        # the inputs that we used to feed the infomap algorithm.
        
        g  <- graph.adjacency(abs(my.matrix),weighted=TRUE)
        
        allowed.interactions.mat <- allowed.interactions(my.matrix, year.plot.matrix,
                                                         guild.combinations[i.guild],
                                                         sp.names, i.year)
        
        linkrank_modularity <- linkRankModularity(g,
                                                  partition=modules$module, 
                                                  allowed.interactions.mat = allowed.interactions.mat,
                                                  damping = lr.damping, 
                                                  pr.algo = lr.alg)
        }else{ # less than 5 sp
          linkrank_modularity <- NA_real_
        }
        # store -------------------------------------------------------------------
        
        pos <- which(community.mean.field.metrics$year == years[i.year] &
                       community.mean.field.metrics$plot == i.plot &
                       community.mean.field.metrics$guilds == guild.combinations[i.guild])
        
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "richness")]] <- my.richness
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "connectance")]] <- my.connectance
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "intra_inter_ratio")]] <- my.ratio
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "degree_distribution")]] <- my.deg.dist
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "modularity")]] <- linkrank_modularity
        # community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "kurtosis")]] <- my.kurtosis
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "skewness")]] <- my.skew
        community.mean.field.metrics$value[pos[which(community.mean.field.metrics$metric[pos] == "complexity")]] <- my.complexity
      }# for i.guild
    }# for i.plot
  }# for i.year
  
  # store mean field --------------------------------------------------------------
  
  write.csv2(community.mean.field.metrics,
             file = paste("results/community_metrics_mean_field",vers,".csv",sep=""))
  
}# if include.mean.field
