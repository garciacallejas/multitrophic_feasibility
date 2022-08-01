
# read community matrices -------------------------------------------------
library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)
library(matlib) # to multiply matrices
library(nleqslv) # to solve non-linear equations
library(zipfR) # beta incomplete function to estimate the area of d-dimensional spherical caps 
library(pracma) # to solve n-dimensional cross products
library(boot) # to bootstrap
library(CValternatives) # to estimate PV index

library(foreach)
library(doParallel)

# source auxiliary functions
list.files("R/feasibility_functions/", full.names = TRUE) %>% map(source)

# set number of cores -----------------------------------------------------

workers <- 10
cl <- makeCluster(workers)
# register the cluster for using foreach
registerDoParallel(cl)

# -------------------------------------------------------------------------
# which version do we read?
# vers <- "_v2"
vers <- ""

# # calculate the observed matrices?
# calculate.obs <- FALSE
# 
# # calculate the feasibility domain of the null/randomised matrices?
# calculate.null <- FALSE

# -------------------------------------------------------------------------

file.name <- "community_matrices"
file.name <- paste(file.name,vers,".RData",sep="")

# community_matrices[[intraguild.type]][[year]][[plot]]
load(paste("results/",file.name,sep=""))

# null.name <- "community_matrices_null"
# null.name <- paste(null.name,vers,".RData",sep="")
# 
# # for null: community_matrices[[intraguild.type]][[replicate]][[year]][[plot]]
# load(paste("results/",null.name,sep=""))

# positions for the blocks are dynamic, species pool is different for 
# each plot/year
load(paste("results/community_names",vers,".RData",sep=""))

intraguild.types <- names(community_matrices)

# set important constants -------------------------------------------------
# number of null replicates
# null.replicates <- length(community_matrices_null[[1]])

# replicates for the feasibility calculations
omega.replicates <- 100
bootstrap.replicates <- 100

# number of noise replicates for the feasibility calculation
# this is for the case in which the original matrix is not invertible
# so that I add a small amount of noise to make it so.
noise.replicates <- 10
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

# this is the combined id to loop over in parallel
id <- expand.grid(intraguild.types,years,plots,guild.combinations)
id.char <- paste(id[,1],"_",id[,2],"_",id[,3],"_",id[,4],sep="")

# calculate feasibility metrics  -------------------------------------------

comb.fun <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

# if(calculate.obs){
  
  feasibility.metrics <- foreach(i.id = 1:length(id.char),
                                 .combine=comb.fun, 
                                 .packages = c("tidyverse","foreach","matlib",
                                               "nleqslv","zipfR","pracma",
                                               "boot","CValternatives",
                                               "MultitrophicFun")) %dopar% {
                                                 
                                                 # cat(id.char[i.id],"- started\n")
                                                 
                                                 list.files("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/feasibility_functions/", 
                                                            full.names = TRUE) %>% map(source)
                                                 
                                                 # first, recover the parameters of each matrix
                                                 
                                                 my.type <- NA
                                                 my.year <- NA
                                                 my.plot <- NA
                                                 my.guild <- NA
                                                 
                                                 # intraguild matrix type
                                                 if(grepl("mean_field",id.char[i.id])){
                                                   my.type <- "mean_field"
                                                 }else if(grepl("nesting_larvae_phenology",id.char[i.id])){
                                                   my.type <- "nesting_larvae_phenology"
                                                 }else if(grepl("nesting_larvae",id.char[i.id])){
                                                   my.type <- "nesting_larvae"
                                                 }else if(grepl("phenology",id.char[i.id])){
                                                   my.type <- "phenology"
                                                 }else if(grepl("nesting",id.char[i.id])){
                                                   my.type <- "nesting"
                                                 }else if(grepl("larvae",id.char[i.id])){
                                                   my.type <- "larvae"
                                                 }
                                                 
                                                 # guild
                                                 if(grepl("plants-floral visitors",id.char[i.id])){
                                                   my.guild <- "plants-floral visitors"
                                                 }else if(grepl("plants-herbivores",id.char[i.id])){
                                                   my.guild <- "plants-herbivores"
                                                 }else if(grepl("floral visitors",id.char[i.id])){
                                                   my.guild <- "floral visitors"
                                                 }else if(grepl("herbivores",id.char[i.id])){
                                                   my.guild <- "herbivores"
                                                 }else if(grepl("plants",id.char[i.id])){
                                                   my.guild <- "plants"
                                                 }else{
                                                   my.guild <- "all"
                                                 }
                                                 
                                                 # year
                                                 if(grepl("2019",id.char[i.id])){
                                                   my.year <- "2019"
                                                 }else if(grepl("2020",id.char[i.id])){
                                                   my.year <- "2020"
                                                 }
                                                 
                                                 # plot
                                                 my.plot1 <- stringr::str_remove(id.char[i.id],my.type)
                                                 my.plot2 <- stringr::str_remove(my.plot1,my.guild)
                                                 my.plot3 <- stringr::str_remove(my.plot2,my.year)
                                                 my.plot <- as.numeric(stringr::str_remove_all(my.plot3,"_"))
                                                 
                                                 # recover the matrix
                                                 year.plot.matrix <- 
                                                   community_matrices[[my.type]][[my.year]][[my.plot]]
                                                 
                                                 plant.positions <- which(rownames(year.plot.matrix) %in% 
                                                                            sp.names[[my.year]][["plants"]])
                                                 herb.positions <- which(rownames(year.plot.matrix) %in% 
                                                                           sp.names[[my.year]][["herbivores"]])
                                                 fv.positions <- which(rownames(year.plot.matrix) %in% 
                                                                         sp.names[[my.year]][["floral.visitors"]])
                                                 
                                                 if(my.guild == "plants"){
                                                   
                                                   my.matrix <- year.plot.matrix[plant.positions,plant.positions]
                                                   
                                                 }else if(my.guild == "plants-floral visitors"){
                                                   
                                                   my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
                                                                                 c(plant.positions,fv.positions)]
                                                   
                                                 }else if(my.guild == "floral visitors"){
                                                   
                                                   my.matrix <- year.plot.matrix[fv.positions,fv.positions]
                                                   
                                                 }else if(my.guild == "herbivores"){
                                                   
                                                   my.matrix <- year.plot.matrix[herb.positions,herb.positions]
                                                   
                                                 }else if(my.guild == "plants-herbivores"){
                                                   
                                                   my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
                                                                                 c(plant.positions,herb.positions)]
                                                   
                                                 }else if(my.guild == "all"){
                                                   
                                                   my.matrix <- year.plot.matrix
                                                 }
                                                 
                                                 # -------------------------------------------------------------------------
                                                 # obtain feasibility domains and species exclusion probabilities
                                                 
                                                 # TODO check with Alfonso
                                                 A <- -my.matrix
                                                 # diag(A) <- 1
                                                 
                                                 my.noise.threshold <- c(0,min(abs(A[which(A != 0)])) * 
                                                                           noise.relative.magnitude)
                                                 if(nrow(A)>2){
                                                   omega.df <- feasibility_metrics(A = A,
                                                                                   omega.replicates = omega.replicates,
                                                                                   bootstrap.replicates = bootstrap.replicates,
                                                                                   noise.replicates = noise.replicates,
                                                                                   noise.threshold = my.noise.threshold,
                                                                                   noise.type = noise.type)
                                                 }else{
                                                   omega.df <- data.frame(omega_mean = NA,
                                                                          omega_lowerCI = NA,
                                                                          omega_upperCI = NA,
                                                                          omega_isotropic = NA,
                                                                          isotropy_index_mean = NA,
                                                                          isotropy_index_lowerCI = NA,
                                                                          isotropy_index_upperCI = NA)
                                                 }
                                                 # this is a test for checking some features of the matrices
                                                 # omega.df <- data.frame(diag.mean = mean(diag(A)),
                                                 #                        diag.pos = sum(diag(A) > 0),
                                                 #                        diag.dom = mean(diag(A)-rowSums(A)),
                                                 #                        n.zeros = sum(A == 0),
                                                 #                        total.strength = sum(A))
                                                 
                                                 omega.df$year <- my.year
                                                 omega.df$plot <- my.plot
                                                 omega.df$guild <- my.guild
                                                 omega.df$intraguild.type <- my.type
                                                 
                                                 if(nrow(A)>2){
                                                   sp.exclusions <- exclusion_probabilities(A = A,
                                                                                            omega.replicates = omega.replicates,
                                                                                            bootstrap.replicates = bootstrap.replicates,
                                                                                            noise.replicates = noise.replicates,
                                                                                            noise.threshold = noise.threshold,
                                                                                            noise.type = noise.type)
                                                 }else{
                                                   sp.exclusions <- data.frame(species = NA,
                                                                               prob_excl_mean = NA,
                                                                               prob_excl_lowerCI = NA,
                                                                               prob_excl_upperCI = NA)
                                                 }
                                                 # sp.exclusions <- data.frame(hey = 1)
                                                 
                                                 sp.exclusions$year <- my.year
                                                 sp.exclusions$plot <- my.plot
                                                 sp.exclusions$guild <- my.guild
                                                 sp.exclusions$intraguild.type <- my.type
                                                 
                                                 # cat(id.char[i.id],"- completed\n")
                                                 # write.csv2(omega.df, file = paste("results/fd_",id.char[i.id],".csv",sep=""))
                                                 # write.csv2(sp.exclusions, file = paste("results/sp_exclusions_",id.char[i.id],".csv",sep=""))
                                                 
                                                 # return
                                                 list(omega.df,sp.exclusions)
                                                 
                                               }
  
  # feasibility.df <- feasibility.metrics[[1]]
  # exclusions.df <- feasibility.metrics[[2]]
  
  # test for displaying some matrix features, 
  # mainly diagonal dominances, diagonal values, 
  # which seem to strongly drive fd
  
  # tt <- feasibility.metrics[[1]] %>%
  #   group_by(guild,intraguild.type) %>%
  #   summarise(mean.zeros = mean(n.zeros),mean.diag = mean(diag.mean),
  #             mean.diag.dom = mean(diag.dom))
  # ggplot(tt, aes(x = intraguild.type, y = mean.diag)) + 
  #   geom_bar(aes(fill = guild),stat = "identity",position=position_dodge()) + 
  #   theme_bw() + 
  #   NULL
  
  # store results -----------------------------------------------------------
  # 
  fd.name <- "feasibility_domain_observed"
  fd.name <- paste("results/",fd.name,vers,".csv",sep="")
  
  exc.name <- "exclusion_probabilities_observed"
  exc.name <- paste("results/",exc.name,vers,".csv",sep="")
  
  write.csv2(x = feasibility.metrics[[1]],file = fd.name,
             row.names = FALSE)
  write.csv2(x = feasibility.metrics[[2]],file = exc.name,
             row.names = FALSE)
  
# }# if calculate.obs

# repeat for null matrices ------------------------------------------------
# if(calculate.null){
#   
#   file.name <- "community_matrices_null"
#   file.name <- paste(file.name,vers,".RData",sep="")
#   
#   load(paste("results/",file.name,sep=""))
#   
#   null.replicates <- length(community_matrices_null[[1]])
#   
#   null.id <- expand.grid(intraguild.types,years,plots,guild.combinations,paste("r",1:null.replicates,sep=""))
#   null.id <- subset(null.id,Var1 != "mean_field")
# 
#   null.id.char <- paste(null.id[,1],"_",null.id[,2],"_",null.id[,3],"_",null.id[,4],"_",null.id[,5],sep="")
#   # null.id.char <- null.id.char[21371:length(null.id.char)]
#   # calculate feasibility metrics  -------------------------------------------
#   
#   # null.id.char <- null.id.char[1:20]
#   
#   null.feasibility.metrics <- foreach(i.id = 1:length(null.id.char),
#                                       # .combine=comb.fun, 
#                                       .packages = c("tidyverse","foreach","matlib",
#                                                     "nleqslv","zipfR","pracma",
#                                                     "boot","CValternatives",
#                                                     "MultitrophicFun")) %dopar% {
#                                                       
#                                                       # cat(id.char[i.id],"- started\n")
#                                                       
#                                                       list.files("/home/david/Work/Projects/EBD/multitrophic_feasibility/R/feasibility_functions/", 
#                                                                  full.names = TRUE) %>% map(source)
#                                                       
#                                                       my.type <- NA
#                                                       my.year <- NA
#                                                       my.plot <- NA
#                                                       my.guild <- NA
#                                                       my.rep <- NA
#                                                       
#                                                       # intraguild matrix type
#                                                       if(grepl("mean_field",null.id.char[i.id])){
#                                                         my.type <- "mean_field"
#                                                       }else if(grepl("nesting_larvae_phenology",null.id.char[i.id])){
#                                                         my.type <- "nesting_larvae_phenology"
#                                                       }else if(grepl("nesting_larvae",null.id.char[i.id])){
#                                                         my.type <- "nesting_larvae"
#                                                       }else if(grepl("phenology",null.id.char[i.id])){
#                                                         my.type <- "phenology"
#                                                       }else if(grepl("nesting",null.id.char[i.id])){
#                                                         my.type <- "nesting"
#                                                       }else if(grepl("larvae",null.id.char[i.id])){
#                                                         my.type <- "larvae"
#                                                       }
#                                                       
#                                                       # guild
#                                                       if(grepl("plants-floral visitors",null.id.char[i.id])){
#                                                         my.guild <- "plants-floral visitors"
#                                                       }else if(grepl("plants-herbivores",null.id.char[i.id])){
#                                                         my.guild <- "plants-herbivores"
#                                                       }else if(grepl("floral visitors",null.id.char[i.id])){
#                                                         my.guild <- "floral visitors"
#                                                       }else if(grepl("herbivores",null.id.char[i.id])){
#                                                         my.guild <- "herbivores"
#                                                       }else if(grepl("plants",null.id.char[i.id])){
#                                                         my.guild <- "plants"
#                                                       }else{
#                                                         my.guild <- "all"
#                                                       }
#                                                       
#                                                       # year
#                                                       if(grepl("2019",null.id.char[i.id])){
#                                                         my.year <- "2019"
#                                                       }else if(grepl("2020",null.id.char[i.id])){
#                                                         my.year <- "2020"
#                                                       }
#                                                       
#                                                       # replicate
#                                                       rep.pos <- gregexpr("_r",null.id.char[i.id])[[1]][1]
#                                                       my.rep <- as.numeric(substr(null.id.char[i.id],rep.pos+2,nchar(null.id.char[i.id])))
#                                                       
#                                                       # plot
#                                                       my.plot1 <- stringr::str_remove(null.id.char[i.id],my.type)
#                                                       my.plot2 <- stringr::str_remove(my.plot1,my.guild)
#                                                       my.plot3 <- stringr::str_remove(my.plot2,my.year)
#                                                       my.plot4 <- stringr::str_remove(my.plot3,paste("_r",my.rep,sep=""))
#                                                       my.plot <- as.numeric(stringr::str_remove_all(my.plot4,"_"))
#                                                       
#                                                       
#                                                       year.plot.matrix <-
#                                                         community_matrices_null[[my.type]][[my.rep]][[my.year]][[my.plot]]
#                                                       
#                                                       if(!any(is.na(year.plot.matrix))){
#                                                         
#                                                         plant.positions <- which(rownames(year.plot.matrix) %in%
#                                                                                    sp.names[[my.year]][["plants"]])
#                                                         herb.positions <- which(rownames(year.plot.matrix) %in%
#                                                                                   sp.names[[my.year]][["herbivores"]])
#                                                         fv.positions <- which(rownames(year.plot.matrix) %in%
#                                                                                 sp.names[[my.year]][["floral.visitors"]])
#                                                         
#                                                         # cat("null rep",i.rep,"-",years[i.year],"-",i.plot,"-",
#                                                         #     guild.combinations[i.guild])
#                                                         
#                                                         if(my.guild == "plants"){
#                                                           
#                                                           my.matrix <- year.plot.matrix[plant.positions,plant.positions]
#                                                           
#                                                         }else if(my.guild == "plants-floral visitors"){
#                                                           
#                                                           my.matrix <- year.plot.matrix[c(plant.positions,fv.positions),
#                                                                                         c(plant.positions,fv.positions)]
#                                                           
#                                                         }else if(my.guild == "floral visitors"){
#                                                           
#                                                           my.matrix <- year.plot.matrix[fv.positions,fv.positions]
#                                                           
#                                                         }else if(my.guild == "herbivores"){
#                                                           
#                                                           my.matrix <- year.plot.matrix[herb.positions,herb.positions]
#                                                           
#                                                         }else if(my.guild == "plants-herbivores"){
#                                                           
#                                                           my.matrix <- year.plot.matrix[c(plant.positions,herb.positions),
#                                                                                         c(plant.positions,herb.positions)]
#                                                           
#                                                         }else if(my.guild == "all"){
#                                                           
#                                                           my.matrix <- year.plot.matrix
#                                                         }
#                                                         
#                                                         # TODO check with Alfonso
#                                                         A <- -my.matrix
#                                                         # diag(A) <- 1
#                                                         
#                                                         my.noise.threshold <- c(0,min(abs(my.matrix[which(my.matrix != 0)])) * 
#                                                                                   noise.relative.magnitude)
#                                                         
#                                                         omega.df <- feasibility_metrics(A = A,
#                                                                                         omega.replicates = omega.replicates,
#                                                                                         bootstrap.replicates = bootstrap.replicates,
#                                                                                         noise.replicates = noise.replicates,
#                                                                                         noise.threshold = my.noise.threshold,
#                                                                                         noise.type = noise.type)
#                                                         # omega.df <- data.frame(id.char = null.id.char[i.id])
#                                                         
#                                                         omega.df$year <- my.year
#                                                         omega.df$plot <- my.plot
#                                                         omega.df$guild <- my.guild
#                                                         omega.df$intraguild.type <- my.type
#                                                         omega.df$replicate <- my.rep
#                                                         
#                                                         sp.exclusions <- exclusion_probabilities(A = A,
#                                                                                                  omega.replicates = omega.replicates,
#                                                                                                  bootstrap.replicates = bootstrap.replicates,
#                                                                                                  noise.replicates = noise.replicates,
#                                                                                                  noise.threshold = noise.threshold,
#                                                                                                  noise.type = noise.type)
#                                                         # sp.exclusions <- data.frame(id.char = null.id.char[i.id])
#                                                         
#                                                         sp.exclusions$year <- my.year
#                                                         sp.exclusions$plot <- my.plot
#                                                         sp.exclusions$guild <- my.guild
#                                                         sp.exclusions$intraguild.type <- my.type
#                                                         sp.exclusions$replicate <- my.rep
#                                                         
#                                                         write.csv2(omega.df, file = paste("results/null_fd/null_fd_",null.id.char[i.id],".csv",sep=""))
#                                                         write.csv2(sp.exclusions, file = paste("results/null_fd/null_sp_exclusions_",null.id.char[i.id],".csv",sep=""))
#                                                         
#                                                         list(omega.df,sp.exclusions)
#                                                       }# if !is.na the matrix
#                                                     }
#   
#   feas.metrics <- sapply(null.feasibility.metrics,function(x) x[1])
#   sp.metrics <- sapply(null.feasibility.metrics,function(x) x[2])
#   
#   feas.df <- bind_rows(feas.metrics)
#   sp.df <- bind_rows(sp.metrics)
#   
#   # store results -----------------------------------------------------------
#   
#   null.fd.name <- "feasibility_domain_null"
#   null.fd.name <- paste("results/",null.fd.name,vers,".csv",sep="")
#   
#   null.exc.name <- "exclusion_probabilities_null"
#   null.exc.name <- paste("results/",null.exc.name,vers,".csv",sep="")
#   
#   write.csv2(x = feas.df,file = null.fd.name,
#              row.names = FALSE)
#   write.csv2(x = sp.df,file = null.exc.name,
#              row.names = FALSE)
#   
# }# if calculate.null

stopCluster(cl)



