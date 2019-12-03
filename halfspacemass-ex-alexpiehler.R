get_normalized_vector <- function(dimensions){

  randn_unif <- runif(dimensions, -1, 1)
  normalized_vector <- matrix(randn_unif, ncol = 1) * (1/norm(randn_unif, type = "2"))
  
}

sample_without_replacement <- function(data, fraction){
  rows_to_sample <- sample(nrow(data), size = floor(nrow(data) * fraction))
  
  data[rows_to_sample, , drop = FALSE]
}

