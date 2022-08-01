
library(MultitrophicFun)

# -------------------------------------------------------------------------

horizontal_community_matrix <- function(S = 5,
                                        c = 0.5,
                                        tau = 1.5,
                                        min.diag.dom = 0,
                                        restricted.positive = TRUE,
                                        int.mean = 0,
                                        int.sd = 1){
  
  a.rows <- S
  a.cols <- S
  l <- round(c * (a.rows*a.cols))
  
  A <- matrix(0,nrow = a.rows,ncol = a.cols)
  if(restricted.positive){
    ints <- abs(gamlss.dist::rSHASHo(l, mu = int.mean,
                                     sigma = int.sd, nu = 0, tau = tau))
  }else{
    ints <- gamlss.dist::rSHASHo(l, mu = int.mean,
                                 sigma = int.sd, nu = 0, tau = tau)
  }
  
  # randomly assign interaction strengths outside the diagonal
  for(i in 1:l){
    my.sample.row <- sample(1:a.rows,1,replace = T)
    my.sample.col <- sample(1:a.cols,1,replace = T)
    
    while(A[my.sample.row,my.sample.col] != 0 &
          my.sample.row == my.sample.col){
      my.sample.row <- sample(1:a.rows,1,replace = T)
      my.sample.col <- sample(1:a.cols,1,replace = T)
    }
    A[my.sample.row,my.sample.col] <- ints[i]
  }# for i
  
  # diag values
  for(i.row in 1:a.rows){
    non.diag <- abs(sum(A[i.row,]))
    if(min.diag.dom > 0){
      # values around that needed to achieve dominance in this row
      A[i.row,i.row] <- abs(rnorm(1,mean = (non.diag + min.diag.dom),sd = .1))
    }else{
      # values from the same distribution as the rest
      A[i.row,i.row] <- abs(gamlss.dist::rSHASHo(1, mu = int.mean,
                                                 sigma = int.sd, nu = 0, tau = tau))
    }
    # cat(i.row,"diag:",A[i.row,i.row],"-non diag:",non.diag,"-dominance:",A[i.row,i.row] - non.diag,"\n")
  }
  return(A)
}

# -------------------------------------------------------------------------
S <- 5
c <- 0.5

A <- round(horizontal_community_matrix(S,c)*100)
bipartite::grouplevel(A,index = "niche overlap",level = "both")
tt <- bipartite::NOS(A,TRUE,FALSE)
library(nos)
t2 <- nos::NOSM_undir(net = nos::freqMat_2_edge(A),sl = 0)
summary(t2)
# I need: 
# species-level

avg_niche_overlap <- function(A,sp = NULL,overlap.type = c("rows",
                                                           "columns",
                                                           "all")){
  
}

