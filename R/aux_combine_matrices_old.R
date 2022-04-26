
aux_combine_matrices <- function(pp.all.years,
                             ph.all.years,
                             pfv.all.years,
                             plant.phenology = NULL,
                             sp.data,
                             taxo.in,
                             include.overlap = TRUE,
                             randomize = FALSE,
                             mean.field.intraguild = FALSE,
                             mean.field.offdiag = NULL,
                             mean.field.diag = NULL){
  
  years <- names(pp.all.years)
  plots <- 1:length(pp.all.years[[1]])
  
  # if matrices are randomized, apply the randomization to the raw data
  if(randomize){
    pp.all.years.null <- list()
    ph.all.years.null <- list()
    pfv.all.years.null <- list()
    
    # r2dtable allows generating a number of replicates in a list, 
    # but doing so one by one I can exactly replicate the code for the
    # observed matrices. Lazy and prone to errors, but i don't expect
    # to further modify this code.
    for(i.year in 1:length(years)){
      
      pp.all.years.null[[i.year]] <- list()
      ph.all.years.null[[i.year]] <- list()
      pfv.all.years.null[[i.year]] <- list()
      
      for(i.plot in plots){
        pp.all.years.null[[i.year]][[i.plot]] <- 
          r2dtable(1,rowSums(pp.all.years[[i.year]][[i.plot]]),
                   colSums(pp.all.years[[i.year]][[i.plot]]))[[1]]
        
        colnames(pp.all.years.null[[i.year]][[i.plot]]) <- 
          colnames(pp.all.years[[i.year]][[i.plot]])
        rownames(pp.all.years.null[[i.year]][[i.plot]]) <- 
          rownames(pp.all.years[[i.year]][[i.plot]])
        
        ph.all.years.null[[i.year]][[i.plot]] <- 
          r2dtable(1,rowSums(ph.all.years[[i.year]][[i.plot]]),
                   colSums(ph.all.years[[i.year]][[i.plot]]))[[1]]
        
        colnames(ph.all.years.null[[i.year]][[i.plot]]) <- 
          colnames(ph.all.years[[i.year]][[i.plot]])
        rownames(ph.all.years.null[[i.year]][[i.plot]]) <- 
          rownames(ph.all.years[[i.year]][[i.plot]])
        
        pfv.all.years.null[[i.year]][[i.plot]] <- 
          r2dtable(1,rowSums(pfv.all.years[[i.year]][[i.plot]]),
                   colSums(pfv.all.years[[i.year]][[i.plot]]))[[1]]
        
        colnames(pfv.all.years.null[[i.year]][[i.plot]]) <- 
          colnames(pfv.all.years[[i.year]][[i.plot]])
        rownames(pfv.all.years.null[[i.year]][[i.plot]]) <- 
          rownames(pfv.all.years[[i.year]][[i.plot]])
      }# for i.plot
    }# for i.year
    
    names(pp.all.years.null) <- years
    names(ph.all.years.null) <- years
    names(pfv.all.years.null) <- years
    
    # bad practice, but i don't want to change all the rest
    pp.all.years <- pp.all.years.null
    ph.all.years <- ph.all.years.null
    pfv.all.years <- pfv.all.years.null
  }
  
  # -------------------------------------------------------------------------
  # convert every intraguild matrix to a square matrix
  # 1 - floral visitors
  fv.all.years <- pfv.all.years
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      A <- pfv.all.years[[i.year]][[i.plot]]
      
      visits_matrix <- matrix(0,nrow = ncol(A),ncol = ncol(A),
                              dimnames = list(colnames(A),colnames(A)))
      # diag: colsum
      # non-diag: first, select non-zero rows (rows with shared visits)
      # second, sum of the minimums (min number of shared visits)
      for(i in 1:nrow(visits_matrix)){
        for(j in 1:ncol(visits_matrix)){
          if(i == j){
            visits_matrix[i,j] <- colSums(A)[i]
          }else{
            # only selected columns
            subm <- A[,c(i,j)]
            # # only non-zero rows
            # subm <- subm[apply(subm,1,FUN = function(x){all(x != 0)}),]
            # note: the above is not necessary, actually
            # sum of the minimum values
            visits_matrix[i,j] <- sum(apply(subm,1,FUN = function(x){min(x)}))
          }
        }# for j
      }# for i
      
      fv.all.years[[i.year]][[i.plot]] <- visits_matrix
      
    }# for i.plot
  }# for i.year
  
  # 2 - herbivores
  h.all.years <- ph.all.years
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      A <- ph.all.years[[i.year]][[i.plot]]
      
      visits_matrix <- matrix(0,nrow = ncol(A),ncol = ncol(A),
                              dimnames = list(colnames(A),colnames(A)))
      # diag: colsum
      # non-diag: first, select non-zero rows (rows with shared visits)
      # second, sum of the minimums (min number of shared visits)
      for(i in 1:nrow(visits_matrix)){
        for(j in 1:ncol(visits_matrix)){
          if(i == j){
            visits_matrix[i,j] <- colSums(A)[i]
          }else{
            # only selected columns
            subm <- A[,c(i,j)]
            # # only non-zero rows
            # subm <- subm[apply(subm,1,FUN = function(x){all(x != 0)}),]
            # note: the above is not necessary, actually
            # sum of the minimum values
            visits_matrix[i,j] <- sum(apply(subm,1,FUN = function(x){min(x)}))
          }
        }# for j
      }# for i
      
      h.all.years[[i.year]][[i.plot]] <- visits_matrix
      
    }# for i.plot
  }# for i.year
  
  # -------------------------------------------------------------------------
  # apply overlap masks before standardizing
  # so that, when standardizing, mathematical properties hold (e.g. all elements sum 1)
  
  if(include.overlap){
    
    # the overlap is slightly different for plants and for animal communities
    # plants have three phenological modes, see "plant.phenology"
    # so that sp in same category have overlap 1, adjacent 0.5, and other 0
    p_overlap <- list()
    
    # the phenological overlap is constant across years for plants
    # i.e. early, middle, and late species are conserved in different years
    
    for(i.year in 1:length(years)){
      
      p_overlap[[i.year]] <- list()
      
      for(i.plot in 1:length(plots)){
        my.plants <- sort(unique(row.names(pp.all.years[[i.year]][[i.plot]])))
        p.overlap.matrix <- matrix(0,nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                                   ncol = ncol(pp.all.years[[i.year]][[i.plot]]),
                                   dimnames = list(my.plants,my.plants))
        
        # fill up plant overlap matrix
        for(i.plant in 1:nrow(p.overlap.matrix)){
          for(j.plant in 1:ncol(p.overlap.matrix)){
            icat <- plant.phenology$pheno.cat[plant.phenology$sp == rownames(p.overlap.matrix)[i.plant]]
            jcat <- plant.phenology$pheno.cat[plant.phenology$sp == rownames(p.overlap.matrix)[j.plant]]
            
            if(icat == jcat){
              my.overlap <- 1
            }else if(icat == "middle" | jcat == "middle"){
              my.overlap <- 0.5
            }else{
              my.overlap <- 0
            }
            
            p.overlap.matrix[i.plant,j.plant] <- my.overlap
            
          }# for j
        }# for i
        
        # weight the plant observations by phenological overlap
        p_overlap[[i.year]][[i.plot]] <- pp.all.years[[i.year]][[i.plot]] * p.overlap.matrix
        
      }# for i.plot
    }# for i.year
    
    names(p_overlap) <- years
    
    # now, animals
    fv_overlap <- list()
    h_overlap <- list()
    
    # phenological and taxonomic potential overlap - for each year
    taxo.matrix <- get_taxonomic_overlap(sp.data = sp.data[sp.data$year == years[1],])
    
    pheno.matrix <- list()
    mask <- list() 
    
    for(i.year in 1:length(years)){
      my.sp.data <- subset(sp.data,!is.na(min.month) & year == years[i.year])
      pheno.matrix[[i.year]] <- get_phenologic_overlap(sp.data = my.sp.data)
      
      # their combination gives the overall potential overlap,
      # that will be used to constrain resource overlap
      year.names <- my.sp.data$sp.ID
      taxo.matrix.year <- taxo.matrix[year.names,year.names]
      
      if(taxo.in){
        mask[[i.year]] <- pheno.matrix[[i.year]] * taxo.matrix.year
      }else{
        mask[[i.year]] <- pheno.matrix[[i.year]]
      }
    }
    
    # here, apply the overlap mask to each pollinator and each herbivore community
    for(i.year in 1:length(years)){
      
      fv_overlap[[i.year]] <- list()
      h_overlap[[i.year]] <- list()
      
      for(i.plot in 1:length(plots)){
        
        plot.fv.names <- colnames(fv.all.years[[i.year]][[i.plot]])
        plot.h.names <- colnames(h.all.years[[i.year]][[i.plot]])
        
        fv.plot.mask <- mask[[i.year]][plot.fv.names,plot.fv.names]
        h.plot.mask <- mask[[i.year]][plot.h.names,plot.h.names]
        
        fv_overlap[[i.year]][[i.plot]] <- fv.all.years[[i.year]][[i.plot]] * fv.plot.mask
        
        h_overlap[[i.year]][[i.plot]] <- h.all.years[[i.year]][[i.plot]] * h.plot.mask
      }# for i.plot
    }# for i.year
    
    names(fv_overlap) <- years
    names(h_overlap) <- years
    
  }else{ # no overlap
    p_overlap <- pp.all.years
    fv_overlap <- fv.all.years
    h_overlap <- h.all.years
  }  
  
  # -------------------------------------------------------------------------
  # now, standardize the intraguild AND interguild matrices
  # there are different options
  
  # here I implement a standardization independent for each community/interaction type
  # so that each matrix is assumed to have the same amount of "interaction effort" 
  # basically meaning each matrix sums to 1
  
  # natural or log values? there is strong variability among them
  # in principle, natural.
  
  pp.all.years.norm <- p_overlap
  fv.all.years.norm <- fv_overlap
  h.all.years.norm <- h_overlap
  pfv.all.years.norm <- pfv.all.years
  ph.all.years.norm <- ph.all.years
  
  for(i.year in 1:length(years)){
    for(i.plot in 1:length(plots)){
      
      # pp.all.years
      pp.all.years.norm[[i.year]][[i.plot]] <- p_overlap[[i.year]][[i.plot]]/sum(p_overlap[[i.year]][[i.plot]])
      
      # fv.all.years
      fv.all.years.norm[[i.year]][[i.plot]] <- fv_overlap[[i.year]][[i.plot]]/sum(fv_overlap[[i.year]][[i.plot]])
      
      # h.all.years
      h.all.years.norm[[i.year]][[i.plot]] <- h_overlap[[i.year]][[i.plot]]/sum(h_overlap[[i.year]][[i.plot]])
      
      # pfv.all.years
      pfv.all.years.norm[[i.year]][[i.plot]] <- pfv.all.years[[i.year]][[i.plot]]/sum(pfv.all.years[[i.year]][[i.plot]])
      
      # ph.all.years
      ph.all.years.norm[[i.year]][[i.plot]] <- ph.all.years[[i.year]][[i.plot]]/sum(ph.all.years[[i.year]][[i.plot]])
      
    }# for i.plot
  }# for i.year
  

# -------------------------------------------------------------------------
# impose a mean-field structure to intra-guild matrices?
  
  if(mean.field.intraguild){
    for(i.year in 1:length(years)){
      for(i.plot in 1:length(plots)){
        
        # plants
        pp.num <- nrow(pp.all.years[[i.year]][[i.plot]]) * 
          ncol(pp.all.years[[i.year]][[i.plot]])
        
        pp.all.years.norm[[i.year]][[i.plot]] <- 
          # matrix(data = rep(mean(pp.all.years[[i.year]][[i.plot]]),pp.num),
          matrix(data = rep(mean.field.offdiag,pp.num),
                 nrow = nrow(pp.all.years[[i.year]][[i.plot]]),
                 dimnames = list(rownames(pp.all.years[[i.year]][[i.plot]]),
                                 colnames(pp.all.years[[i.year]][[i.plot]])))
        
        diag(pp.all.years.norm[[i.year]][[i.plot]]) <- mean.field.diag
        
        # floral visitors and herbivores
        plot.fv.names <- colnames(fv.all.years[[i.year]][[i.plot]])
        plot.h.names <- colnames(h.all.years[[i.year]][[i.plot]])
        
        fv.all.years.norm[[i.year]][[i.plot]] <- 
              matrix(data = rep(mean.field.offdiag,length(plot.fv.names)^2),
                     nrow = length(plot.fv.names),
                     dimnames = list(plot.fv.names,plot.fv.names))

            diag(fv.all.years.norm[[i.year]][[i.plot]]) <- mean.field.diag

            h.all.years.norm[[i.year]][[i.plot]] <-
              matrix(data = rep(mean.field.offdiag,length(plot.h.names)^2),
                     nrow = length(plot.h.names),
                     dimnames = list(plot.h.names,plot.h.names))

            diag(h.all.years.norm[[i.year]][[i.plot]]) <- mean.field.diag
        
      }
    }
    
  }
  
  # impose signs ------------------------------------------------------------
  
  # only positive should be plant-floral visitor
  # and herb-plant!!
  pp_sign <- -1
  pfv_sign <- 1
  ph_sign <- -1
  hp_sign <- 1
  hh_sign <- -1
  fvfv_sign <- -1
  
  # combine all matrices ----------------------------------------------------
  
  community_matrices <- list()
  sp.names <- list()
  
  for(i.year in 1:length(years)){
    
    # year-plot list, each element a block matrix
    community_matrices[[i.year]] <- list()
    
    # also store the set of plants, herbivores, floral visitors
    # per year
    plants.year <- NULL
    herbivores.year <- NULL
    floral.visitors.year <- NULL
    
    for(i.plot in 1:length(plots)){
      
      # components of the block matrix
      pp.matrix <- pp_sign * pp.all.years.norm[[i.year]][[i.plot]]
      ph.matrix <- ph_sign * ph.all.years.norm[[i.year]][[i.plot]]
      pfv.matrix <- pfv_sign * pfv.all.years.norm[[i.year]][[i.plot]]
      fv.overlap.matrix <- fvfv_sign * fv.all.years.norm[[i.year]][[i.plot]]
      h.overlap.matrix <- hh_sign * h.all.years.norm[[i.year]][[i.plot]]
      
      # get the set of species for this community
      valid.names <- get_valid_sp(pp.matrix = pp.matrix,
                                  ph.matrix = ph.matrix,
                                  pfv.matrix = pfv.matrix)
      
      # construct the empty matrix
      block.matrix <- build_block_matrix(plants = valid.names[[1]],
                                         herbivores = valid.names[[2]],
                                         floral.visitors = valid.names[[3]])
      
      # fill the matrix
      # Note that in the function "fill_block_matrix" there is an argument
      # switch_herb_sign that is set to TRUE by default (so not specified here), 
      # that sets plant effects on herbivores as positive, but not the 
      # other way around
      block.matrix.filled <- fill_block_matrix(block.matrix = block.matrix,
                                               pp.matrix = pp.matrix,
                                               ph.matrix = ph.matrix,
                                               pfv.matrix = pfv.matrix,
                                               fv.overlap.matrix = fv.overlap.matrix,
                                               h.overlap.matrix = h.overlap.matrix)
      
      # store the matrix
      community_matrices[[i.year]][[i.plot]] <- block.matrix.filled
      
      # store the valid names of this community
      plants.year <- unique(c(plants.year,valid.names[[1]]))
      herbivores.year <- unique(c(herbivores.year,valid.names[[2]]))
      floral.visitors.year <- unique(c(floral.visitors.year,valid.names[[3]]))
    }
    
    sp.names[[i.year]] <- list(plants = sort(plants.year),
                               herbivores = sort(herbivores.year),
                               floral.visitors = sort(floral.visitors.year))
  }
  
  names(community_matrices) <- years
  names(sp.names) <- years
  
  return(list(community_matrices = community_matrices,
       sp.names = sp.names))
  
}
