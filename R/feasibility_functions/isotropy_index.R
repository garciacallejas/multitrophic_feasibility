isotropy_index <- function(probability_vector){
  
  number_species <- length(probability_vector)
  
  probability_vector_positive <- probability_vector[probability_vector > 0]
  
  log_prob <- log(probability_vector_positive)
  
  return(-sum(probability_vector_positive*log_prob)/log(number_species))
  
}
