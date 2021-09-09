
# combine plant-plant, plant-floral visitor, plant-herbivore matrices
# in a block matrix. Also, calculate resource overlap matrices for
# herbivores and floral visitors
# 
# In parallel, store all names of the different species of each guild

library(tidyverse)
# devtools::install_github("RadicalCommEcol/MultitrophicFun",
# auth_token = "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")
library(MultitrophicFun)

source("R/aux_combine_matrices.R")

# -------------------------------------------------------------------------
# add a version suffix?
# vers <- "_v2"
vers <- ""

# which type of overlap - if any?
include.overlap <- TRUE

# include phylogenetic overlap?
# this means in practice that, if included, only taxa from the same order
# are assumed to compete among them
taxo.in <- FALSE

# compute null matrices? --------------------------------------------------
# null matrices will be based on reshuffling visits while keeping
# col and row totals fixed

# also, compute mean field matrices?
# aij will be the same for every i,j, mean(aij). I do not expect
# other aij values would make a difference.
# note that mean field will be produced for within-guild matrices only.

include.null <- TRUE
replicates <- 100

include.mean.field <- TRUE
mean.field.value <- .2 # Saavedra et al. 2013
mean.field.diag <- 1

# read data ---------------------------------------------------------------

years <- c(2019,2020)
plots <- 1:9

plant.phenology <- read.csv2("data/plant_phenology_categories.csv")
sp.data <- read.csv2(file = paste("data/species_phenology_taxonomy",vers,".csv",sep=""),
                     stringsAsFactors = FALSE)

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

cm <- aux_combine_matrices(pp.all.years = pp.all.years,
                           ph.all.years = ph.all.years,
                           pfv.all.years = pfv.all.years,
                           plant.phenology = plant.phenology,
                           sp.data = sp.data,
                           taxo.in = taxo.in,
                           include.overlap = TRUE,
                           randomize = FALSE) 

# store block matrix ------------------------------------------------------

file.name <- "community_matrices_observed"
file.name <- paste("results/",file.name,vers,".Rdata",sep="")

community_matrices <- cm[[1]]
sp.names <- cm[[2]]

save(community_matrices,
     file = file.name)

save(sp.names,
     file = paste("results/community_names",vers,".RData",sep=""))

# repeat for null matrices ------------------------------------------------

if(include.null){
  
  null_matrices <- list()
  
  for(i.rep in 1:replicates){
    
    cm.null <- aux_combine_matrices(pp.all.years = pp.all.years,
                                    ph.all.years = ph.all.years,
                                    pfv.all.years = pfv.all.years,
                                    plant.phenology = plant.phenology,
                                    sp.data = sp.data,
                                    taxo.in = taxo.in,
                                    include.overlap = include.overlap,
                                    randomize = TRUE)
    
    null_matrices[[i.rep]] <- cm.null[[1]]
    
  }# for i.replicate
  
  # store block matrix ------------------------------------------------------
  
  file.name <- "community_matrices_null"
  file.name <- paste("results/",file.name,vers,".Rdata",sep="")
  
  save(null_matrices,
       file = file.name)
  
}# if include.null

# repeat for mean field matrices ------------------------------------------

if(include.mean.field){
  
    cmmf <- aux_combine_matrices(pp.all.years = pp.all.years,
                                 ph.all.years = ph.all.years,
                                 pfv.all.years = pfv.all.years,
                                 plant.phenology = plant.phenology,
                                 sp.data = sp.data,
                                 taxo.in = taxo.in,
                                 # overlap is not compatible with mean-field
                                 # so it does not matter what value it takes here
                                 include.overlap = FALSE,
                                 randomize = FALSE,
                                 mean.field.intraguild = TRUE,
                                 mean.field.offdiag = mean.field.value,
                                 mean.field.diag = mean.field.diag)
    
  # store block matrix ------------------------------------------------------
  
  file.name <- "community_matrices_mean_field"
  file.name <- paste("results/",file.name,vers,".Rdata",sep="")
  
  mean_field_matrices <- cmmf[[1]]
  
  save(mean_field_matrices,
       file = file.name)
  
}# if include.mean.field
