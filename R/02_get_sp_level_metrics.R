# obtain a couple of species-level metrics for
# the exclusion probabilities regression

load(file = paste("results/community_matrices.RData",sep=""))

# load sp names as well
load(file = paste("results/community_names.RData",sep=""))
plants <- sort(unique(c(sp.names[["2019"]][["plants"]],sp.names[["2020"]][["plants"]])))
fv <- sort(unique(c(sp.names[["2019"]][["floral.visitors"]],sp.names[["2020"]][["floral.visitors"]])))
herb <- sort(unique(c(sp.names[["2019"]][["herbivores"]],sp.names[["2020"]][["herbivores"]])))

years <- c(2019,2020)
plots <- 1:9

guild.combinations <- c("plants",
                        "floral visitors",
                        "herbivores",
                        "plants-floral visitors",
                        "plants-herbivores",
                        "all")

intraguild.types <- names(community_matrices)

# -------------------------------------------------------------------------
sp.metrics.list <- list()

# test
# i.type <- i.year <- i.plot <- i.guild <- 2

for(i.type in 1:length(intraguild.types)){
  for(i.year in 1:length(years)){
    for(i.plot in plots){
      
      year.plot.matrix <- 
        community_matrices[[i.type]][[as.character(years[i.year])]][[i.plot]]
      
      plant.positions <- which(rownames(year.plot.matrix) %in% 
                                 sp.names[[i.year]][["plants"]])
      herb.positions <- which(rownames(year.plot.matrix) %in% 
                                sp.names[[i.year]][["herbivores"]])
      fv.positions <- which(rownames(year.plot.matrix) %in% 
                              sp.names[[i.year]][["floral.visitors"]])
      
      for(i.guild in 1:length(guild.combinations)){
        
        cat(intraguild.types[i.type],"-",years[i.year],"-",i.plot,"-",guild.combinations[i.guild],"\n")
        
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
        
        present.sp <- rownames(my.matrix)
        
        df.sp <- expand.grid(year = years[i.year],
                            plot = i.plot,
                            sp.guild = NA,
                            guilds = guild.combinations[i.guild],
                            intraguild.type = intraguild.types[i.type],
                            species = present.sp,
                            diagonal_dominance = NA,
                            in_degree = NA,
                            out_degree = NA,stringsAsFactors = FALSE)
        
        # diagonal dominance does not make sense for bipartite matrices
        if(guild.combinations[i.guild] %in% c("plants","floral visitors","herbivores","all")){
          # diagonal dominance
          # diagonal elements
          d <- diag(my.matrix)
          # rowsums
          nd <- rowSums(my.matrix)
          # but without the diagonal
          nd <- nd - d
          
          df.sp$diagonal_dominance <- abs(d) - abs(nd)
        }
        
        # to harmonize uni and bipartite nets, better to go one sp at a time
        for(i.sp in 1:nrow(df.sp)){
          
          # first, assign guild
          my.sp <- df.sp$species[i.sp]
          
          if(my.sp %in% plants){
            df.sp$sp.guild[i.sp] <- "plants"
          }else if(my.sp %in% fv){
            df.sp$sp.guild[i.sp] <- "floral visitors"
          }else if(my.sp %in% herb){
            df.sp$sp.guild[i.sp] <- "herbivores"
          }
          
          # exclude the diagonal
          df.sp$in_degree[i.sp] <- sum(my.matrix[my.sp,] != 0) - 1
          df.sp$out_degree[i.sp] <- sum(my.matrix[,my.sp] != 0) - 1
        }
        
        sp.metrics.list[[length(sp.metrics.list)+1]] <- df.sp
        
      }# i.guild
    }# i.plot
  }# i.year
}# i.type

sp.metrics <- bind_rows(sp.metrics.list)

# -------------------------------------------------------------------------
write.csv2(sp.metrics,"results/species_level_metrics.csv",row.names = F)


