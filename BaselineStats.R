# Mann-Whitney (Wilcoxon Rank-Sum) test
# gripper.csv
# first runs
p.ffpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/ffpo-gripper.csv"
p.cgpo <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cgpo-gripper.csv"
p.cg <- "/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cg-gripper.csv"
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
dat.ffpo$instance <- 15 + 5*dat.ffpo$instance
dat.cgpo$instance <- 15 + 5*dat.cgpo$instance
dat.cg$instance <- 15 + 5*dat.cg$instance

wilcox.test(dat.ffpo$runtime, dat.cgpo$runtime, alternative="g")
wilcox.test(dat.ffpo$runtime, dat.cg$runtime, alternative="g")
wilcox.test(dat.cg$runtime, dat.cgpo$runtime, alternative="g")

##################################################

# Coverage data
# total solved by each planner from all domains

setwd("/Users/calebhill/Documents/UNH/MinorProject/code/")
files <- list.files(path=paste(getwd(), "/planner_data/data_1/", sep=""))
cov.cg <- 0
cov.cgpo <- 0
cov.ffpo <- 0
total.instances <- 30*length(files)/3
for(f in files){
  if(substring(f, 1,3)=="cg-"){
    dat <- read.csv( paste(getwd(),"/planner_data/data_1/",f,sep="") )
    cov.cg <- cov.cg + length(dat$instance)
  }
  else if(substring(f, 1,3)=="cgp"){
    dat <- read.csv( paste(getwd(),"/planner_data/data_1/",f,sep="") )
    cov.cgpo <- cov.cgpo + length(dat$instance)
  }
  else if(substring(f, 1,3)=="ffp"){
    dat <- read.csv( paste(getwd(),"/planner_data/data_1/",f,sep="") )
    cov.ffpo <- cov.ffpo + length(dat$instance)
  }
  rm(dat,f)
}


##################################################

# McNemar's test
# cg vs cgpo

p.ffpo.grip <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/ffpo-gripper.csv")
p.ffpo.block <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/ffpo-blocksworld.csv")
p.ffpo.child <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/ffpo-childsnack.csv")
p.ffpo.park <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/ffpo-parking.csv")
p.ffpo.scan <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/ffpo-scanalyzer.csv")
p.cgpo.grip <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cgpo-gripper.csv")
p.cgpo.block <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cgpo-blocksworld.csv")
p.cgpo.child <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cgpo-childsnack.csv")
p.cgpo.park <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cgpo-parking.csv")
p.cgpo.scan <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cgpo-scanalyzer.csv")
p.cg.grip <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_1/cg-gripper.csv")
p.cg.block <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cg-blocksworld.csv")
p.cg.child <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cg-childsnack.csv")
p.cg.park <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cg-parking.csv")
p.cg.scan <- get.dat("/Users/calebhill/Documents/UNH/MinorProject/code/planner_data/data_2/cg-scanalyzer.csv")

#        cg+ cg-
# ffpo+  a11 a12
# ffpo-  a21 a22
a11 <- 0; a12 <- 0; a21 <- 0; a22 <- 0
# gripper
for(i in 1:30){
  if( (i %in% p.cg.grip$instance) & (i %in% p.cffpo.grip$instance)){ a11 <- a11+1}
  else if( (i %in% p.cg.grip$instance) & !(i %in% p.ffpo.grip$instance)){ a21 <- a21+1}
  else if( !(i %in% p.cg.grip$instance) & (i %in% p.ffpo.grip$instance)){ a12 <- a12+1}
  else if( !(i %in% p.cg.grip$instance) & !(i %in% p.ffpo.grip$instance)){ a22 <- a22+1}
}
# parking
for(i in 1:30){
  if( (i %in% p.cg.park$instance) & (i %in% p.ffpo.park$instance)){ a11 <- a11+1}
  else if( (i %in% p.cg.park$instance) & !(i %in% p.ffpo.park$instance)){ a21 <- a21+1}
  else if( !(i %in% p.cg.park$instance) & (i %in% p.ffpo.park$instance)){ a12 <- a12+1}
  else if( !(i %in% p.cg.park$instance) & !(i %in% p.ffpo.park$instance)){ a22 <- a22+1}
}
# blocksworld
for(i in 1:30){
  if( (i %in% p.cg.block$instance) & (i %in% p.ffpo.block$instance)){ a11 <- a11+1}
  else if( (i %in% p.cg.block$instance) & !(i %in% p.ffpo.block$instance)){ a21 <- a21+1}
  else if( !(i %in% p.cg.block$instance) & (i %in% p.ffpo.block$instance)){ a12 <- a12+1}
  else if( !(i %in% p.cg.block$instance) & !(i %in% p.ffpo.block$instance)){ a22 <- a22+1}
}
# childsnack
for(i in 1:30){
  if( (i %in% p.cg.child$instance) & (i %in% p.ffpo.child$instance)){ a11 <- a11+1}
  else if( (i %in% p.cg.child$instance) & !(i %in% p.ffpo.child$instance)){ a21 <- a21+1}
  else if( !(i %in% p.cg.child$instance) & (i %in% p.ffpo.child$instance)){ a12 <- a12+1}
  else if( !(i %in% p.cg.child$instance) & !(i %in% p.ffpo.child$instance)){ a22 <- a22+1}
}

# scanalyzer
for(i in 1:30){
  if( (i %in% p.cg.scan$instance) & (i %in% p.ffpo.scan$instance)){ a11 <- a11+1}
  else if( (i %in% p.cg.scan$instance) & !(i %in% p.ffpo.scan$instance)){ a21 <- a21+1}
  else if( !(i %in% p.cg.scan$instance) & (i %in% p.ffpo.scan$instance)){ a12 <- a12+1}
  else if( !(i %in% p.cg.scan$instance) & !(i %in% p.ffpo.scan$instance)){ a22 <- a22+1}
}


# perform the test
conting <- matrix(data=c(a11,a21,a12,a22), ncol=2,nrow=2)
mcnemar.test(conting)




# size adjustment for parking
x = c(9,10,11,12,13,14,16,17,18,19,20,22,23,24,25,26,28,29,30,31,32,34,35,36,37,39,40,41,42,43)
p.ffpo.park$instance[1:8] <- x[1:8]
p.ffpo.park$instance[9] <- x[10]
for(i in 1:length(p.ffpo.park$instance)){
  rho1 <- p.ffpo.park$instance[i]
  p.ffpo.park$instance[i] <- 1+ sqrt( (9-rho1)^2 + (16-(2*rho1-2))^2 )
  rm(rho1)
}
p.cgpo.park$instance <- p.ffpo.park$instance[1:8]
p.cg.park$instance <- p.ffpo.park$instance[1:8]

# size adjustment for gripper
p.cgpo.grip$instance <- 15 + 5*p.cgpo.grip$instance
p.cg.grip$instance <- 15 + 5*p.cg.grip$instance

stats.park <- SimulStats(dat.ffpo=p.ffpo.park, dat.cgpo=p.cgpo.park, dat.cg=p.cg.park, num.iters = 1000)
stats.grip <- SimulStats(dat.ffpo=p.ffpo.park, dat.cgpo=p.cgpo.park, dat.cg=p.cg.park, num.iters = 1000)
SlopeCompare(stats.park$m2, stats.park$m3) 
SlopeCompare(stats.grip$m2, stats.grip$m3) 
