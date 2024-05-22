library(tidyverse)

wd = "-- your wd here --"

# get planner dataframe and sort
# takes:
#   path to planner data
#returns:
#   sorted dataframe 
get.dat <- function(path){
  
  dat <- read.csv(path)
  dat.sorted <- dat[order(dat$instance),]
  return(dat.sorted)
}