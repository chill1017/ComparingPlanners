library(tidyverse)

wd = "-- your wd here --"

# probability of m1 < m2
# takes: (equal length)
#   vector of m1 values
#   vector of m2 values
# returns:
#   approximation to P( m1<m2 )

ParamCompare <- function(m1, m2){
  N <- length(m1)^2
  n <- length(m1)
  c <- 0
  
  for(i in 1:n){
    for(j in 1:n){
      if(m1[i] < m2[j]){ c <- c+1}
    }
  }
  return(c/N)
}