# everything needed to do stuff for all things

# arguments
###################################################################
###################################################################
###################################################################
# This script is an example of how to parse command line arguments when running an R script at
# the linux command line using the Rscript binary.
# It has the following command line usage syntax:
# $ Rscript get_args.r <DOMAIN> <FFPO_FILENAME> <CGPO_FILENAME> <CG_FILENAME>
# For example, the following is correct usage:
# $ Rscript get_args.r gripper ffpo-gripper.csv cgpo-gripper.csv cg-gripper.csv
args = commandArgs()

# One file for each of ffpo, cgpo, cg
expected_number_of_filenames = 3

# Rscript adds in extra args after the binary name (which is always the first arg), so we have to count from the back
index_of_first_filename = length(args) - expected_number_of_filenames + 1
index_of_domain_name = index_of_first_filename - 1

# Use this variable to determine which problem size ordering to use.
# (If R does not have a dictionary/hashtable kind of datastructure, 
# then you might just need several if-statements, one to test this 
# variable for the name of each domain for which we have a problem size order for so far.)
domain_name = args[index_of_domain_name]

write("\nGot this domain name:", stdout())
write(domain_name, stdout())

ffpo_data_filename = args[index_of_first_filename]
cgpo_data_filename = args[index_of_first_filename + 1]
cg_data_filename = args[index_of_first_filename + 2]

write("\nGot these filenames:", stdout())

write("\nfor ffpo:", stdout())
write(ffpo_data_filename, stdout())
stopifnot(startsWith(ffpo_data_filename, "ffpo-"))

write("\nfor cgpo:", stdout())
write(cgpo_data_filename, stdout())
stopifnot(startsWith(cgpo_data_filename, "cgpo-"))

write("\nfor cg:", stdout())
write(cg_data_filename, stdout())
stopifnot(startsWith(cg_data_filename, "cg-"))

# For the future, but not pre-WIPC (I tried but R is ennervating):
# TODO ensure all three filenames are on the same domain
# TODO stop hard-coding algorithms, read in variable number of files, 
#   parse alg name and domain name from each file, make list of alg names, 
#   set of domain names, ensure set of domain names has cardinality 1. 
#   (This would change command line argument interface for the script, 
#   and batch system would have to be updated accordingly.)

install.packages("fitdistrplus")
library(fitdistrplus)

# functions
###################################################################
###################################################################
###################################################################

# Planner analysis function that takes:
# dataframes
# number of bootstrap iterations
# initial m1, m2, m3, k values
# and returns:
# df with
# m1, m2, m3, k values
SimulStats <- function(dat.ffpo, dat.cgpo, dat.cg, num.iters=500, rand=FALSE){ #, init.vals=c(0.1, 0.1, 0.1, 0.1)){
  x1 <- dat.ffpo$instance
  x2 <- dat.cgpo$instance
  x3 <- dat.cg$instance
  
  y1 <- dat.ffpo$runtime
  y2 <- dat.cgpo$runtime
  y3 <- dat.cg$runtime
  
  # https://www.dataquest.io/blog/how-to-create-a-dataframe-in-r/
  output <- data.frame(
    matrix(rep(0,4*num.iters),
           nrow=num.iters,
           ncol=4)
  )
  names(output) <- c("m.ffpo", "m.cgpo", "m.cg", "k")
  
  if(rand==FALSE){
    # prelim nls for m1, m2, m3, k initial values
    fit1.init <- nls(y1~m1*2^(k1*x1), data=dat.ffpo, start=list(m1=0.1,k1=0.1))
    m1.init <- summary(fit1.init)$coefficients[1]
    k1.init <- summary(fit1.init)$coefficients[2]
    
    fit2.init <- nls(y2~m2*2^(k2*x2), data=dat.ffpo, start=list(m2=0.1,k2=0.1))
    m2.init <- summary(fit2.init)$coefficients[1]
    k2.init <- summary(fit2.init)$coefficients[2]
    
    fit3.init <- nls(y3~m3*2^(k3*x3), data=dat.ffpo, start=list(m3=0.1,k3=0.1))
    m3.init <- summary(fit3.init)$coefficients[1]
    k3.init <- summary(fit3.init)$coefficients[2]
    
    k.found <- mean(k1.init,k2.init,k3.init)
    init.vals <- c(m1.init, m2.init, m3.init, k.found)
    
  }else{
    init.vals <- runif(4,0.01,0.5)
    k.found <- init.vals[4]
  }
  # simultaneous cost function
  BigCost <- function(params){
    m1 <- params[1]
    m2 <- params[2]
    m3 <- params[3]
    k <- params[4]                                 
    c <- 0
    # cost for ffpo
    c1 <- 0
    len <- length(dat.ffpo)
    for(i in 1:len){
      x <- dat.ffpo$instance[i]
      y <- dat.ffpo$runtime[i]
      c1 <- c1 + (y - ( m1*2^(k*x) ) )^2
    }
    #cost for cgpo
    c2 <- 0
    len <- length(dat.cgpo)
    for(i in 1:len){
      x <- dat.cgpo$instance[i]
      y <- dat.cgpo$runtime[i]
      c2 <- c2 + (y - ( m2*2^(k*x) ) )^2
    }
    #cost for cg
    c3 <- 0
    len <- length(dat.cg)
    for(i in 1:len){
      x <- dat.cg$instance[i]
      y <- dat.cg$runtime[i]
      c3 <- c3 + (y - ( m3*2^(k*x) ) )^2           
    }
    c <- c1+c2+c3                      
    return(c)
  }
  # first parameter estimation
  first <- optim(par=init.vals,
                 fn=BigCost,
                 method=c("L-BFGS"),
                 lower=c(0.000000001, 0.000000001, 0.000000001, 0.000000001)
                )
  params.found <- first$par
  m1.found <- params.found[1]
  m2.found <- params.found[2]
  m3.found <- params.found[3]
  k.found <- params.found[4]              
  
  
  
  
  
  # switch to percentage error!
  
  # residual distributions
  y1.est <- m1.found*2^(k.found*x1); y2.est <- m2.found*2^(k.found*x2); y3.est <- m3.found*2^(k.found*x3);
  res1 <- y1/y1.est; res2 <- y2/y2.est; res3 <- y3/y3.est;
  # mu1 <- mean(res1); mu2 <- mean(res2); mu3 <- mean(res3);
  # sd1 <- sqrt(var(res1)); sd2 <- sqrt(var(res2)); sd3 <- sqrt(var(res3)); 
  fit.gamma1 <- fitdist(res1,distr="gamma",method="mle")
  fit.gamma2 <- fitdist(res2,distr="gamma",method="mle")
  fit.gamma3 <- fitdist(res3,distr="gamma",method="mle")
  
  
  # printing residuals for debugging
  #pdf("./plots/resid-plots.pdf")
  #hist(res1, col="blue", xlim=c(0,50), probability=TRUE, breaks=30)
  #hist(res2, col="red", probability=FALSE,breaks=10, add=TRUE)
  #hist(res3, col="green", probability=FALSE,breaks=10, add=TRUE)

  # n.domain <- seq(-1,1,length=50)
  # n.fit1 <- dnorm(n.domain,mean=mu1,sd=sd1)
  # lines(n.domain,n.fit1,col="blue",lwd=2)
  # n.fit2 <- dnorm(n.domain,mean=mu2,sd=sd2)
  # lines(n.domain,n.fit2,col="red",lwd=2)
  # n.fit3 <- dnorm(n.domain,mean=mu3,sd=sd3)
  # lines(n.domain,n.fit3,col="green",lwd=2)
  #dev.off()

  # bootstrap for param distribution
  for(i in 1:num.iters){
    # new simulated data each time
    y1.new <- m1.found * 2^(x1*k.found)*rgamma(length(x1),                   #+ rnorm(length(x1),mu1,sd1)
                                               shape=fit.gamma1$estimate[1],
                                               rate=fit.gamma1$estimate[2]
                                               ) 
    dat.ffpo.new <- dat.ffpo
    dat.ffpo.new$runtime <- y1.new
    
    y2.new <- m2.found * 2^(x2*k.found) *rgamma(length(x2),                   #+ rnorm(length(x2),mu2,sd2)
                                                shape=fit.gamma2$estimate[1],
                                                rate=fit.gamma2$estimate[2]
                                                )
    dat.cgpo.new <- dat.cgpo
    dat.cgpo.new$runtime <- y2.new
    
    y3.new <- m3.found * 2^(x3*k.found) *rgamma(length(x3),                   #+ rnorm(length(x3),mu3,sd3)
                                                shape=fit.gamma3$estimate[1],
                                                rate=fit.gamma3$estimate[2]
                                                )
    dat.cg.new <- dat.cg
    dat.cg.new$runtime <- y3.new
    
    # find new params
    model.iter <- optim(par=runif(4,0.01,.9),               
                        fn=BigCost,
                        method=c("L-BFGS"),
                        lower=c(0.001, 0.001, 0.001, 0.001)  
    )
    params.i <- model.iter$par
    output[i,] <- params.i
    
    # plot simulated data and its fit
    pdf(paste("./plots/simulated_data_",i,".pdf",sep=""))
    crv <- function(x){params.i[1]*2^(x*params.i[4])}
    plot(x1,y1.new,xlab="size",ylab="simulated runtime")
    curve(crv, add=TRUE)
    dev.off()
    
    # output total loss from fit for debugging
    loss <- sum( abs(y1.new-params.i[1]*2^(x1*params.i[4])))+ sum( abs(y2.new-params.i[2]*2^(x2*params.i[4])))+ sum( abs(y3.new-params.i[3]*2^(x3*params.i[4])))
    print(loss)
  }
  return(output)
}
###################################################################

# probability of m1 < m2
# takes:
#   vector of m1 values
#   vector of m2 values
# returns:
#   approximation to P( m1<m2 )

# assumes equal length
SlopeCompare <- function(m1, m2){
  N <- length(m1)^2
  n <- length(m1)
  c <- 0
  
  # loop over all n^2
  # count cases where m1<m2
  # return total/n^2
  for(i in 1:n){
    for(j in 1:n){
      if(m1[i] < m2[j]){ c <- c+1}
    }
  }
  return(c/N)
}
###################################################################

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

# parameter values
###################################################################
###################################################################
###################################################################
params.blocksworld <- c(9,13,18,23,28,33,38,43,48,53,57,62,67,72,77,82,87,92,97,102,107,111,116,121,126,131136,141,146,151)
params.childsnack <- c(5,6,7,9,10,12,13,14,16,17,19,20,21,23,24,25,27,28,30,31,32,34)
params.gripper <- c(20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165)
params.scanalyzer <- c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
x = c(9,10,11,12,13,14,16,17,18,19,20,22,23,24,25,26,28,29,30,31,32,34,35,36,37,39,40,41,42,43)
params.parking <- 1:30
for(i in 1:30){ params.parking[i] <- 1+ sqrt( (9-x[i])^2 + (16-(2*x[i]-2))^2 )}
rm(x)

# grab data
dat.ffpo <- get.dat(ffpo_data_filename); dat.cgpo <- get.dat(cgpo_data_filename); dat.cg <- get.dat(cg_data_filename)
# transform sizes
if(domain_name=="blocksworld"){
  for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.blocksworld[i]}
  for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.blocksworld[i]}
  for(i in dat.cg$instance){dat.cg$instance[i] <- params.blocksworld[i]}
}else if(domain_name=="childsnack"){
  for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.childsnack[i]}
  for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.childsnack[i]}
  for(i in dat.cg$instance){dat.cg$instance[i] <- params.childsnack[i]}
}else if(domain_name=="gripper"){
  for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.gripper[i]}
  for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.gripper[i]}
  for(i in dat.cg$instance){dat.cg$instance[i] <- params.gripper[i]}
}else if(domain_name=="parking"){
  for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.parking[i]}
  for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.parking[i]}
  for(i in dat.cg$instance){dat.cg$instance[i] <- params.parking[i]}  
}else if(domain_name=="scanalyzer"){
  for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.scanalyzer[i]}
  for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.scanalyzer[i]}
  for(i in dat.cg$instance){dat.cg$instance[i] <- params.scanalyzer[i]}  
}else{
  stop("invalid domain.")
}


# slope params csv, histograms, runtime plots
###################################################################
###################################################################
###################################################################

# stats
stats <- SimulStats(dat.ffpo=dat.ffpo, dat.cgpo=dat.cgpo, dat.cg=dat.cg, num.iters=1000)
# export slope and scale parameters
write.csv(stats, paste("slope_scale_params-", domain_name,".csv"), row.names=FALSE)

# histogram
# smallest rectangle needed
h1 <- hist(stats$m1, plot=FALSE, breaks=50); h1$counts <- h1$counts/sum(h1$counts)
h2 <- hist(stats$m2, plot=FALSE, breaks=50); h2$counts <- h2$counts/sum(h2$counts)
h3 <- hist(stats$m3, plot=FALSE, breaks=50); h3$counts <- h3$counts/sum(h3$counts)
x.upper <- max(h1$breaks, h2$breaks, h3$breaks)
y.upper <- max(h1$counts, h2$counts, h3$counts)
# plot
pdf( paste("slope_hist-",domain_name,".pdf") )
p1 <- plot(h1,xlim=c(0,x.upper), ylim=c(0,y.upper), col="blue", xlab="", ylab="",main="", axes=TRUE)
p2 <- plot(h2,xlim=c(0,x.upper), ylim=c(0,y.upper), col="red", xlab="", ylab="", add=TRUE)
p3 <- plot(h3,xlim=c(0,x.upper), ylim=c(0,y.upper), col="green", xlab="", ylab="", add=TRUE)
legend("topright",
       col=c("blue","red", "green"),
       legend=c("ffpo", "cgpo", "cg"),
       pch=c(15,15,15)
      )
dev.off()
rm(x.upper, y.upper)


# runtimes
# smallest rectangle
x.upper <- max(dat.ffpo$instance, dat.cgpo$instance, dat.cg$instance)
y.upper <- max(dat.ffpo$runtime, dat.cgpo$runtime, dat.cg$runtime)
# plot
pdf( paste("runtimes-", domain_name, ".pdf") )
plot(dat.ffpo$instance, dat.ffpo$runtime, xlim=c(0,x.upper), ylim=c(0,y.upper), col="blue", xlab="size", ylab="runtime", main="")
points(dat.cgpo$instance, dat.cgpo$runtime, col="red")
points(dat.cg$instance, dat.cg$runtime, col="green")
legend("topleft",
       col=c("blue","red", "green"),
       legend=c("ffpo", "cgpo", "cg"),
       pch=c(1,1,1)
)
dev.off()
rm(x.upper, y.upper)

# power plot
###################################################################
###################################################################
###################################################################


ffpo.vs.cgpo <- rep(NA,30)
cgpo.vs.cg <- rep(NA,30)
cg.vs.ffpo <- rep(NA,30)
t.s <- Sys.time()
for(i in 1:30){
  tryCatch(
    stats <- SimulStats(dat.ffpo=dat.ffpo[1:i],
                      dat.cg=dat.cg[1:i,],
                      dat.cgpo=dat.cgpo[1:i,],
                      num.iters=500,
                      rand=TRUE),
    error=function(e){}
  )
  ffpo.vs.cgpo[i] <- SlopeCompare(stats$m.ffpo, stats$m.cgpo)
  cgpo.vs.cg[i] <- SlopeCompare(stats$m.cgpo, stats$m.cg)
  cg.vs.ffpo[i] <- SlopeCompare(stats$m.cg, stats$m.ffpo)
}
t.e <- Sys.time()
t.e-t.s

NonNAindex.ffpo.vs.cgpo <- which(!is.na(ffpo.vs.cgpo))
lastNonNA.ffpo.vs.cgpo <- max(NonNAindex.ffpo.vs.cgpo)
NonNAindex.cgpo.vs.cg <- which(!is.na(cgpo.vs.cg))
lastNonNA.cgpo.vs.cg <- max(NonNAindex.cgpo.vs.cg)
NonNAindex.cg.vs.ffpo <- which(!is.na(cg.vs.ffpo))
lastNonNA.cg.vs.ffpo <- max(NonNAindex.cg.vs.ffpo)
max.instance <- max(lastNonNA.ffpo.vs.cgpo, lastNonNA.cgpo.vs.cg, lastNonNA.cg.vs.ffpo)

pdf( paste("quasi_power-",domain_name,".pdf") )
plot(1:max.instance,ffpo.vs.cgpo[1:max.instance], ylim=c(0,1.4), xlab="problem instance number", ylab="probability", col="black")
points(1:max.instance, cgpo.vs.cg[1:max.instance], col="seagreen1")
points(1:max.instance, cg.vs.ffpo[1:max.instance], col="darkorange2")
legend("topleft",
       col=c("black","seagreen1", "darkorange2"),
       legend=c("P(m_ffpo < m_cgpo)", "P(m_cgpo < m_cg)", "P(m_cg < m_ffpo)"),
       pch=c(1,1)
)







