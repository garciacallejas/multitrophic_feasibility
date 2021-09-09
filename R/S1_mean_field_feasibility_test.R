
# test the feasibility domain of mean-field communities of increasing richness
library(ggplot2)
library(MultitrophicFun)
rows <- 2:50

# vers <- "_v2"
vers <- ""

df <- data.frame(size = rows, fd = NA_real_)

for(i.row in rows){
  B <- matrix(.2,nrow = i.row,ncol = i.row)
  diag(B) <- 1
  
  my.fd <- get_matrix_feasibility(B,bootstrap.replicates = 10)
    
  df$fd[which(df$size == i.row)] <- mean(my.fd$feasibility.domain)
}

fd.plot <- ggplot(df,aes(x = size, y = fd)) + 
  geom_point() +
  theme_bw() + 
  labs(x = "community size", y = "feasibility domain") +
  NULL

ggsave(filename = paste("results/images/Fig_S3",vers,".pdf",sep=""),plot = fd.plot,
       device = cairo_pdf,
       width = 5,height = 5,dpi = 300)

