
# combine plant-plant, plant-floral visitor, plant-herbivore matrices
# in a block matrix. Also, calculate intraguild matrices for
# herbivores and floral visitors
# 
# In parallel, store all names of the different species of each guild

library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun")
library(MultitrophicFun)

source("R/aux_combine_matrices.R")

# -------------------------------------------------------------------------
# add a version suffix?
# vers <- "_v2"
vers <- ""

# -------------------------------------------------------------------------
# which intraguild matrix types?
intraguild.types <- c("mean.field","phenology","phenology_nesting",
                      "phenology_larvae")


# -------------------------------------------------------------------------
# mean field coefficients
mean.field.offdiag <- .2 # Saavedra et al. 2013
mean.field.diag <- 1

# compute null matrices? --------------------------------------------------
# null matrices will be based on reshuffling visits while keeping
# col and row totals fixed
# null matrices are obtained for each intraguild matrix type

include.null <- TRUE
replicates <- 100

# read data ---------------------------------------------------------------

years <- c(2019,2020)
plots <- 1:9

plant.phenology <- read.csv2("data/plant_phenology_categories.csv")
animal.phenology <- read.csv2(file = paste("data/species_phenology_taxonomy",vers,".csv",sep=""),
                     stringsAsFactors = FALSE)

animal.info <- read.csv2("data/species_nest_larval_info.csv")
animal.nesting.info <- animal.info[,c("ID","nesting")]
animal.larval.info <- animal.info[,c("ID","larval.food.requirements")]

pp.all.years <- list()
ph.all.years <- list()
pfv.all.years <- list()

for(i.year in 1:length(years)){
  
  load(paste("./data/plant_plant_matrices_",
             years[i.year],vers,".RData",sep=""))
  load(paste("./data/plant_floral_visitor_matrices_",
             years[i.year],vers,".RData",sep=""))
  load(paste("./data/plant_herbivore_matrices_",
             years[i.year],vers,".RData",sep=""))
  
  pp.all.years[[i.year]] <- p_p
  ph.all.years[[i.year]] <- p_h
  pfv.all.years[[i.year]] <- p_fv
}

names(pp.all.years) <- years
names(ph.all.years) <- years
names(pfv.all.years) <- years

# remove empty rows and columns -------------------------------------------
# keep only species that appear in a given year and plot

for(i.year in 1:length(years)){
  for(i.plot in 1:length(plots)){
    for(i.guild in c("pp","ph","pfv")){
      
      if(i.guild == "pp"){
        my.matrix <- pp.all.years[[i.year]][[i.plot]]
        
        my.valid.rows <- apply(my.matrix,1,sum)
        my.valid.rows <- names(my.valid.rows)[which(my.valid.rows != 0)]
        
        my.valid.cols <- apply(my.matrix,2,sum)
        my.valid.cols <- names(my.valid.cols)[which(my.valid.cols != 0)]
        
        # slightly different from ph,pfv, because this is a plant-plant
        # matrix, so one sp cannot be only on rows/cols.
        my.valid.sp <- intersect(my.valid.rows,my.valid.cols)
        my.matrix <- my.matrix[my.valid.sp,my.valid.sp]
        
        pp.all.years[[i.year]][[i.plot]] <- my.matrix
        
      }else if(i.guild == "ph"){
        my.matrix <- ph.all.years[[i.year]][[i.plot]]
        
        my.valid.rows <- apply(my.matrix,1,sum)
        my.valid.rows <- names(my.valid.rows)[which(my.valid.rows != 0)]
        
        my.valid.cols <- apply(my.matrix,2,sum)
        my.valid.cols <- names(my.valid.cols)[which(my.valid.cols != 0)]
        
        ph.all.years[[i.year]][[i.plot]] <- my.matrix[my.valid.rows,
                                                      my.valid.cols]
        
      }else if(i.guild == "pfv"){
        my.matrix <- pfv.all.years[[i.year]][[i.plot]]
        
        my.valid.rows <- apply(my.matrix,1,sum)
        my.valid.rows <- names(my.valid.rows)[which(my.valid.rows != 0)]
        
        my.valid.cols <- apply(my.matrix,2,sum)
        my.valid.cols <- names(my.valid.cols)[which(my.valid.cols != 0)]
        
        pfv.all.years[[i.year]][[i.plot]] <- my.matrix[my.valid.rows,
                                                      my.valid.cols]
      }# if i.guild
    }# for i.guild
  }# for i.plot
}# for i.year

# -------------------------------------------------------------------------
# obtain the block matrix associated with each intraguild matrix type

community_matrices <- list()

for(i.type in 1:length(intraguild.types)){
  
  community_matrices[[i.type]] <- list()
  
  # -------------------------------------------------------------------------
  # obtain "observed" matrix of i.type
  
  my.observed.matrix <- aux_combine_matrices(pp.all.years = pp.all.years,
                                             ph.all.years = ph.all.years,
                                             pfv.all.years = pfv.all.years,
                                             plant.phenology = plant.phenology,
                                             animal.phenology = animal.phenology,
                                             animal.nesting.info = animal.nesting.info,
                                             animal.larval.info = animal.larval.info,
                                             randomize = FALSE,
                                             intraguild.type = intraguild.types[i.type],
                                             mean.field.offdiag = mean.field.offdiag,
                                             mean.field.diag = mean.field.diag)

  # -------------------------------------------------------------------------
  community_matrices[[i.type]][[1]] <- my.observed.matrix[[1]]
  
  # retrieve the names as well
  # this only needs to be done once, as species composition does not change
  if(i.type == 1){
    sp.names <- my.observed.matrix[[2]]
  }
  
  if(include.null){
    
    community_matrices[[i.type]][[2]] <- list()
    
    for(i.rep in 1:replicates){
      
      # obtain null block matrix of i.type
      my.null.matrix <- aux_combine_matrices(pp.all.years = pp.all.years,
                                                 ph.all.years = ph.all.years,
                                                 pfv.all.years = pfv.all.years,
                                                 plant.phenology = plant.phenology,
                                                 animal.phenology = animal.phenology,
                                                 animal.nesting.info = animal.nesting.info,
                                                 animal.larval.info = animal.larval.info,
                                                 randomize = TRUE,
                                                 intraguild.type = intraguild.types[i.type],
                                                 mean.field.offdiag = mean.field.offdiag,
                                                 mean.field.diag = mean.field.diag)
      
      # add to the list
      community_matrices[[i.type]][[2]][[i.rep]] <- my.null.matrix[[1]]
      
    }# for i.rep
    
    names(community_matrices[[i.type]]) <- c("observed","null.replicates")
    
  }# if include.null
}# for i.type
names(community_matrices) <- intraguild.types

# -------------------------------------------------------------------------

# save(community_matrices,
#      file = paste("results/community_matrices",vers,"RData",sep=""))

# save(sp.names,
#      file = paste("results/community_names",vers,".RData",sep=""))


# -------------------------------------------------------------------------
# old code, delete when ready

# cm <- aux_combine_matrices(pp.all.years = pp.all.years,
#                            ph.all.years = ph.all.years,
#                            pfv.all.years = pfv.all.years,
#                            plant.phenology = plant.phenology,
#                            animal.phenology = animal.phenology,
#                            animal.nesting.info = animal.nesting.info,
#                            animal.larval.info = animal.larval.info,
#                            randomize = FALSE,
#                            intraguild.type = "phenology")
# 
# # store block matrix ------------------------------------------------------
# 
# file.name <- "community_matrices_observed"
# file.name <- paste("results/",file.name,vers,".Rdata",sep="")
# 
# community_matrices <- cm[[1]]
# sp.names <- cm[[2]]
# 
# save(community_matrices,
#      file = file.name)
# 
# save(sp.names,
#      file = paste("results/community_names",vers,".RData",sep=""))
# 
# # repeat for null matrices ------------------------------------------------
# 
# if(include.null){
#   
#   null_matrices <- list()
#   
#   for(i.rep in 1:replicates){
#     
#     cm.null <- aux_combine_matrices(pp.all.years = pp.all.years,
#                                     ph.all.years = ph.all.years,
#                                     pfv.all.years = pfv.all.years,
#                                     plant.phenology = plant.phenology,
#                                     sp.data = sp.data,
#                                     taxo.in = taxo.in,
#                                     include.overlap = include.overlap,
#                                     randomize = TRUE)
#     
#     null_matrices[[i.rep]] <- cm.null[[1]]
#     
#   }# for i.replicate
#   
#   # store block matrix ------------------------------------------------------
#   
#   file.name <- "community_matrices_null"
#   file.name <- paste("results/",file.name,vers,".Rdata",sep="")
#   
#   save(null_matrices,
#        file = file.name)
#   
# }# if include.null
# 
# # repeat for mean field matrices ------------------------------------------
# 
# if(include.mean.field){
#   
#     cmmf <- aux_combine_matrices(pp.all.years = pp.all.years,
#                                  ph.all.years = ph.all.years,
#                                  pfv.all.years = pfv.all.years,
#                                  plant.phenology = plant.phenology,
#                                  sp.data = sp.data,
#                                  taxo.in = taxo.in,
#                                  # overlap is not compatible with mean-field
#                                  # so it does not matter what value it takes here
#                                  include.overlap = FALSE,
#                                  randomize = FALSE,
#                                  mean.field.intraguild = TRUE,
#                                  mean.field.offdiag = mean.field.offdiag,
#                                  mean.field.diag = mean.field.diag)
#     
#   # store block matrix ------------------------------------------------------
#   
#   file.name <- "community_matrices_mean_field"
#   file.name <- paste("results/",file.name,vers,".Rdata",sep="")
#   
#   mean_field_matrices <- cmmf[[1]]
#   
#   save(mean_field_matrices,
#        file = file.name)
#   
# }# if include.mean.field
