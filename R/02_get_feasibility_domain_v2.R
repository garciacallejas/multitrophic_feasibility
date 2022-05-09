
# read community matrices -------------------------------------------------
library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)

# -------------------------------------------------------------------------
# which version do we read?
# vers <- "_v2"
vers <- ""

# calculate the feasibility domain of the null/randomised matrices?
calculate.null <- TRUE

# -------------------------------------------------------------------------

file.name <- "community_matrices"
file.name <- paste(file.name,vers,".RData",sep="")

# community_matrices[[intraguild.type]][[year]][[plot]]
load(paste("results/",file.name,sep=""))

null.name <- "community_matrices_null"
null.name <- paste(null.name,vers,".RData",sep="")

# for null: community_matrices[[intraguild.type]][[replicate]][[year]][[plot]]
load(paste("results/",null.name,sep=""))

# positions for the blocks are dynamic, species pool is different for 
# each plot/year
load(paste("results/community_names",vers,".RData",sep=""))

intraguild.types <- names(community_matrices)

# set important constants -------------------------------------------------
# number of null replicates
null.replicates <- length(community_matrices_null[[1]])

# number of bootstrap replicates for the feasibility calculation
bootstrap.replicates <- 10#100
noise.type <- "random"
# how much noise to add relative to the minimum observed values
# e.g. an order of magnitude of 1% 
# relative to the minimum of the matrix elements
noise.relative.magnitude <- .1

years <- as.numeric(names(community_matrices[[1]]))
plots <- 1:9
guild.combinations <- c("plants","floral visitors","herbivores",
                        "plants-floral visitors","plants-herbivores",
                        "all")

# if(calculate.obs){
  
  # define results ----------------------------------------------------------

fd <- list()
extinction.probs <- list()

# calculate feasibility domain  -------------------------------------------
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
        
        cat(intraguild.types[i.type],"-",years[i.year],"-",i.plot,"-",guild.combinations[i.guild])
        
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
        
        A <- -my.matrix
        # diag(A) <- 1
        
        my.noise.threshold <- c(0,min(abs(my.matrix[which(my.matrix != 0)])) * 
          noise.relative.magnitude)
        
        my.matrix.feas <- get_matrix_feasibility(A = A,
                                                 bootstrap.replicates = bootstrap.replicates,
                                                 noise.threshold = my.noise.threshold,
                                                 noise.type = noise.type)
        
        # TODO Update with Alfonso's
        my.iso.feas <- 0
        
        fd[[length(fd)+1]] <- data.frame(year = years[i.year],
                                         plot = i.plot,
                                         guild = guild.combinations[i.guild],
                                         intraguild.type = intraguild.types[i.type],
                                         fd.average = mean(my.matrix.feas$feasibility.domain),
                                         fd.sd = sd(my.matrix.feas$feasibility.domain),
                                         iso.fd = my.iso.feas)
        
        cat(" ...completed\n")
      }# for i.guilds
    }# for i.plot
  }# for i.year
}# for i.type - intraguild types
  
feasibility.df <- bind_rows(fd)
  
# store results -----------------------------------------------------------

file.name <- "feasibility_domain_observed"
file.name <- paste("results/",file.name,vers,".csv",sep="")

write.csv2(x = feasibility.df,file = file.name,
           row.names = FALSE)

# repeat for null matrices ------------------------------------------------
if(calculate.null){
  
  file.name <- "community_matrices_null"
  file.name <- paste(file.name,vers,".RData",sep="")
  
  load(paste("results/",file.name,sep=""))
  
  # define results ----------------------------------------------------------
  
  null.fd <- list()
  
  # calculate feasibility domain  -------------------------------------------
  for(i.type in 1:intraguild.types){
    for(i.rep in 1:null.replicates){
      for(i.year in 1:length(years)){
        for(i.plot in plots){
          
          year.plot.matrix <-
            null_matrices[[i.type]][[i.rep]][[as.character(years[i.year])]][[i.plot]]
          
          if(!is.null(year.plot.matrix)){
            
            plant.positions <- which(rownames(year.plot.matrix) %in%
                                       sp.names[[i.year]][["plants"]])
            herb.positions <- which(rownames(year.plot.matrix) %in%
                                      sp.names[[i.year]][["herbivores"]])
            fv.positions <- which(rownames(year.plot.matrix) %in%
                                    sp.names[[i.year]][["floral.visitors"]])
            
            for(i.guild in 1:length(guild.combinations)){
              
              cat("null rep",i.rep,"-",years[i.year],"-",i.plot,"-",
                  guild.combinations[i.guild])
              
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
              
              A <- -my.matrix
              # diag(A) <- 1
              
              my.matrix.feas <- get_matrix_feasibility(A = A,
                                                       bootstrap.replicates = bootstrap.replicates,
                                                       noise.threshold = noise.threshold[[i.year]][[i.plot]],
                                                       noise.type = noise.type)
              
              # cat(" fd:",round(mean(my.matrix.feas$feasibility.domain),3))
              
              null.fd[[length(fd)+1]] <- data.frame(year = years[i.year],
                                                    plot = i.plot,
                                                    guild = guild.combinations[i.guild],
                                                    replicate = i.rep,
                                                    intraguild.type = intraguild.types[i.type],
                                                    fd.average = mean(my.matrix.feas$feasibility.domain),
                                                    fd.sd = sd(my.matrix.feas$feasibility.domain),
                                                    iso.fd = my.iso.feas)
              
              cat(" ...completed\n")
            }# if !is.null
            
          }# for i.guilds
        }# for i.plot
      }# for i.year
    }# for i.rep
  }# for i.type
  feasibility.df.null <- bind_rows(fd)
  
  # store results -----------------------------------------------------------
  
  file.name <- "feasibility_domain_null"
  file.name <- paste("results/",file.name,vers,".csv",sep="")
  
  write.csv2(x = feasibility.df.null,file = file.name,
             row.names = FALSE)
  
}# if calculate.null
