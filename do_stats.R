# v 4
setwd(" ")


# functions
###################################################################
###################################################################
###################################################################

# Planner analysis function that takes:
#   dataframes
#   number of bootstrap iterations
#   rand - initial fit with random values?
# and returns:
# df with
# m1, m2, m3, k values
SimulStats <- function(dat.ffpo, dat.cgpo, dat.cg, num.iters=500, rand=FALSE, clean.thresh=0.5, plot=FALSE){
  dat.ffpo.clean <- subset(dat.ffpo, runtime>clean.thresh)
  dat.cgpo.clean <- subset(dat.cgpo, runtime>clean.thresh)
  dat.cg.clean <- subset(dat.cg, runtime>clean.thresh)
  
  x1 <- dat.ffpo.clean$instance
  x2 <- dat.cgpo.clean$instance
  x3 <- dat.cg.clean$instance
  
  y1 <- dat.ffpo.clean$runtime
  y2 <- dat.cgpo.clean$runtime
  y3 <- dat.cg.clean$runtime
  
  output <- data.frame(
    matrix(rep(0,6*num.iters),
           nrow=num.iters,
           ncol=6)
  )
  names(output) <- c("m.ffpo", "m.cgpo", "m.cg", "k.ffpo", "k.cgpo", "k.cg")
  
  if(rand==FALSE){
    # prelim nls for m1, m2, m3, k1, k2, k3 initial values
    ldat.ffpo <- dat.ffpo.clean; ldat.ffpo$runtime <- log(ldat.ffpo$runtime)
    fit1.init <- lm(ldat.ffpo$runtime ~ ldat.ffpo$instance, data=ldat.ffpo)
    m1.init <- exp( fit1.init$coeff[[1]] )
    k1.init <- fit1.init$coeff[[2]]/log(2)
    
    ldat.cgpo <- dat.cgpo.clean; ldat.cgpo$runtime <- log(ldat.cgpo$runtime)
    fit2.init <- lm(ldat.cgpo$runtime ~ ldat.cgpo$instance, data=ldat.cgpo)
    m2.init <- exp( fit2.init$coeff[[1]] )
    k2.init <- fit2.init$coeff[[2]]/log(2)
    
    ldat.cg <- dat.cg.clean; ldat.cg$runtime <- log(ldat.cg$runtime)
    fit3.init <- lm(ldat.cg$runtime ~ ldat.cg$instance, data=ldat.cg)
    m3.init <- exp( fit3.init$coeff[[1]] )
    k3.init <- fit3.init$coeff[[2]]/log(2)
    
    init.vals <- c(m1.init, m2.init, m3.init, k1.init, k2.init, k3.init)
  }else{
    init.vals <- runif(6,0.01,0.5)
  }
  
  
  # residual distributions
  y1.est <- m1.init*2^(k1.init*x1)
  y2.est <- m2.init*2^(k2.init*x2)
  y3.est <- m3.init*2^(k3.init*x3)
  res1 <- log(y1)-log(y1.est)
  res2 <- log(y2)-log(y2.est)
  res3 <- log(y3)-log(y3.est)
  mu1 <- mean(res1)
  mu2 <- mean(res2)
  mu3 <- mean(res3)
  sd1 <- sqrt(var(res1))
  sd2 <- sqrt(var(res2))
  sd3 <- sqrt(var(res3))
  
  # plot residuals for debugging
  # plot(x1,exp(res1), col="blue", xlab="problem size", ylab="multiplicative deviation", main=paste("multiplicative residuals",domain_name))
  # points(x2,exp(res2), col="red",add=TRUE )
  # points(x3,exp(res3), col="green", add=TRUE )
  # legend("topleft",
  #        col=c("blue","red", "green"),
  #        legend=c("ffpo", "cgpo", "cg"),
  #        pch=c(1,1,1)
  # )
  
  # printing residuals for debugging
  if(plot){
    d1 <- density(res1)
    d2 <- density(res2)
    d3 <- density(res3)
    x.range <- range(d1$x, d2$x, d3$x)
    y.range <- range(d1$y, d2$y, d3$y)
    pdf(paste("./plots/resid_plots-",domain_name,".pdf",sep=""))
    plot(d1, col="blue", xlim=x.range, ylim=y.range, main=paste(domain_name, "residuals: ffpo"))
    lines(d2,col="red", main=paste(domain_name, "residuals: cgpo"))
    lines(d3,col="green", main=paste(domain_name, "residuals: cg"))
    dev.off()
  }
  
  # bootstrap for param distribution
  for(i in 1:num.iters){
    # new simulated data each time
    ldat.ffpo.iter <- dat.ffpo.clean
    ldat.ffpo.iter$runtime <- log(ldat.ffpo.iter$runtime) + rnorm(length(x1),mean=mu1,sd=sd1)
    fit1.iter <- lm(ldat.ffpo.iter$runtime ~ ldat.ffpo.iter$instance, data=ldat.ffpo.iter)
    m1.iter <- exp( fit1.iter$coeff[[1]] )
    k1.iter <- fit1.iter$coeff[[2]]/log(2)
    
    ldat.cgpo.iter <- dat.cgpo.clean
    ldat.cgpo.iter$runtime <- log(ldat.cgpo.iter$runtime) + rnorm(length(x2),mean=mu2,sd=sd2)
    fit2.iter <- lm(ldat.cgpo.iter$runtime ~ ldat.cgpo.iter$instance, data=ldat.cgpo.iter)
    m2.iter <- exp( fit2.iter$coeff[[1]] )
    k2.iter <- fit2.iter$coeff[[2]]/log(2)
    
    ldat.cg.iter <- dat.cg.clean
    ldat.cg.iter$runtime <- log(ldat.cg.iter$runtime) + rnorm(length(x3),mean=mu3,sd=sd3)
    fit3.iter <- lm(ldat.cg.iter$runtime ~ ldat.cg.iter$instance, data=ldat.cg.iter)
    m3.iter <- exp( fit3.iter$coeff[[1]] )
    k3.iter <- fit3.iter$coeff[[2]]/log(2)
    
    output[i,] <- c(m1.iter, m2.iter, m3.iter, k1.iter, k2.iter, k3.iter)
    
    # plot simulated data and its fit
    if(plot & i<=100){ 
      x.upper <- max(dat.ffpo.clean$instance, dat.cgpo.clean$instance, dat.cg.clean$instance)
      y.upper <- max(dat.ffpo$runtime, dat.cgpo$runtime, dat.cg$runtime)
      crv1 <- function(x){m1.iter*2^(x*k1.iter)}
      crv2 <- function(x){m2.iter*2^(x*k2.iter)}
      crv3 <- function(x){m3.iter*2^(x*k3.iter)}
      pdf(paste("./plots/simulated_data_",i,".pdf",sep=""))
      plot(x1,exp(ldat.ffpo.iter$runtime), xlim=c(0,x.upper), ylim=c(0,y.upper), col="blue",
           xlab="size", ylab="simulated runtime", main=paste(domain_name,"bootstrap iter:",i))
      points(x2,exp(ldat.cgpo.iter$runtime), col="red")
      points(x3,exp(ldat.cg.iter$runtime), col="green")
      curve(crv1,col="blue",add=TRUE)
      curve(crv2,col="red",add=TRUE)
      curve(crv3,col="green",add=TRUE)
      legend("topleft",
             col=c("blue","red", "green"),
             legend=c("ffpo", "cgpo", "cg"),
             pch=c(1,1,1)
      )
      dev.off()
    }
    
  }
  if(plot){
    library(qpdf)
    simulation.plots <- rep("",min(num.iters,100))
    for(i in 1:min(num.iters,100)){ simulation.plots[i] <- paste("./plots/simulated_data_",i,".pdf",sep="") }
    qpdf::pdf_combine(input=simulation.plots,
                      output=paste("./combined_plots/simulated-",domain_name,".pdf",sep="")
    )
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
ParamCompare <- function(m1, m2){
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
domain_name <- "parking"
p.ffpo <- " "
p.cgpo <- " "
p.cg <- " "
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
x = c(9,10,11,12,13,14,16,17,18,19,20,22,23,24,25,26,28,29,30,31,32,34,35,36,37,39,40,41,42,43)
dat.ffpo$instance[1:8] <- x[1:8]
dat.ffpo$instance[9] <- x[10]
for(i in 1:length(dat.ffpo$instance)){
  rho1 <- dat.ffpo$instance[i]
  dat.ffpo$instance[i] <- 1+ sqrt( (9-rho1)^2 + (16-(2*rho1-2))^2 )
}
dat.cgpo$instance <- dat.ffpo$instance[1:8]
dat.cg$instance <- dat.ffpo$instance[1:8]
rm(p.ffpo, p.cgpo, p.cg,i,x,rho1)


domain_name <- "gripper"
p.ffpo <- " "
p.cgpo <- " "
p.cg <- " "
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
dat.ffpo$instance <- 15 + 5*dat.ffpo$instance
dat.cgpo$instance <- 15 + 5*dat.cgpo$instance
dat.cg$instance <- 15 + 5*dat.cg$instance
rm(p.ffpo, p.cgpo, p.cg)


domain_name <- "scanalyzer"
p.ffpo <- " "
p.cgpo <- " "
p.cg <- " "
dat.ffpo <- get.dat(p.ffpo)
dat.cgpo <- get.dat(p.cgpo)
dat.cg <- get.dat(p.cg)
params.scanalyzer <- c(4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33)
for(i in dat.ffpo$instance){dat.ffpo$instance[i] <- params.scanalyzer[i]}
for(i in dat.cgpo$instance){dat.cgpo$instance[i] <- params.scanalyzer[i]}
for(i in dat.cg$instance){dat.cg$instance[i] <- params.scanalyzer[i]}
rm(p.ffpo, p.cgpo, p.cg, i)


# params csv, distributions, runtime plots
###################################################################
###################################################################
###################################################################

# stats
stats <- SimulStats(dat.ffpo=dat.ffpo, dat.cgpo=dat.cgpo, dat.cg=dat.cg,
                    num.iters=1000, plot=TRUE, clean.thresh=0.25)

# scale parameter distributions
dk1 <- density(stats$k.ffpo)
dk2 <- density(stats$k.cgpo)
dk3 <- density(stats$k.cg)
x.range <- range(dk1$x, dk2$x, dk3$x)
y.range <- range(dk1$y, dk2$y, dk3$y)
# plot
scale.distr.path <- paste("./plots/scale_param_dist-",domain_name,".pdf",sep="")
pdf( scale.distr.path )
plot(dk1, lwd=2, col="blue", xlab="scale parameter", ylab="probability",main=paste("scale parameters for ", domain_name,sep=""),
     xlim=x.range, ylim=y.range)
lines(dk2, lwd=2, col="red", xlab="", ylab="", add=TRUE)
lines(dk3, lwd=2, col="green", xlab="", ylab="", add=TRUE)
legend("topright",
       col=c("blue","red", "green"),
       legend=c("ffpo", "cgpo", "cg"),
       pch=c(15,15,15)
)
dev.off()
rm(dk1,dk2,dk3, x.range,y.range)

# slope parameter distributions
dm1 <- density(stats$m.ffpo)
dm2 <- density(stats$m.cgpo)
dm3 <- density(stats$m.cg)
x.range <- range(dm1$x, dm2$x, dm3$x)
y.range <- range(dm1$y, dm2$y, dm3$y)
# plot
slope.distr.path <- paste("./plots/slope_param_dist-",domain_name,".pdf",sep="")
pdf( slope.distr.path )
plot(dm1, lwd=2, col="blue", xlab="slope parameter", ylab="probability",main=paste("slope parameters for ", domain_name,sep=""),
     xlim=x.range, ylim=y.range)
lines(dm2, lwd=2, col="red", xlab="", ylab="", add=TRUE)
lines(dm3, lwd=2, col="green", xlab="", ylab="", add=TRUE)
legend("topright",
       col=c("blue","red", "green"),
       legend=c("ffpo", "cgpo", "cg"),
       pch=c(15,15,15)
)
dev.off()
rm(dm1,dm2,dm3, x.range,y.range)

# combine slope and scale PDFs
qpdf::pdf_combine(input=c(slope.distr.path,scale.distr.path),
                  output=paste("./combined_plots/param_distr-",domain_name,".pdf",sep=""))





# actual runtimes
# smallest rectangle
x.upper <- max(dat.ffpo$instance, dat.cgpo$instance, dat.cg$instance)
y.upper <- max(dat.ffpo$runtime, dat.cgpo$runtime, dat.cg$runtime)
# plot
actual.runtime.path <- paste("./plots/actual_runtimes-", domain_name, ".pdf", sep="")
pdf( actual.runtime.path )
plot(dat.ffpo$instance, dat.ffpo$runtime, xlim=c(0,x.upper), ylim=c(0,y.upper), col="blue",
     xlab="size", ylab="runtime", main=paste("runtimes on ", domain_name,sep=""))
points(dat.cgpo$instance, dat.cgpo$runtime, col="red")
points(dat.cg$instance, dat.cg$runtime, col="green")
legend("topleft",
       col=c("blue","red", "green"),
       legend=c("ffpo", "cgpo", "cg"),
       pch=c(1,1,1)
)
dev.off()
rm(x.upper, y.upper)







# "power" plot
###################################################################
###################################################################
###################################################################

# our method
ffpo.vs.cgpo.k <- rep(NA,30)
cgpo.vs.cg.k <- rep(NA,30)
cg.vs.ffpo.k <- rep(NA,30)
ffpo.vs.cgpo.m <- rep(NA,30)
cgpo.vs.cg.m <- rep(NA,30)
cg.vs.ffpo.m <- rep(NA,30)
for(i in 4:30){
  stats.iter <- SimulStats(dat.ffpo=dat.ffpo[1:i,],
                           dat.cg=dat.cg[1:i,],
                           dat.cgpo=dat.cgpo[1:i,],
                           num.iters=1000, rand=FALSE, clean.thresh=0.01, plot=FALSE)
  ffpo.vs.cgpo.k[i] <- ParamCompare(stats.iter$k.ffpo, stats.iter$k.cgpo)
  cgpo.vs.cg.k[i] <- ParamCompare(stats.iter$k.cgpo, stats.iter$k.cg)
  cg.vs.ffpo.k[i] <- ParamCompare(stats.iter$k.cg, stats.iter$k.ffpo)
  
  ffpo.vs.cgpo.m[i] <- ParamCompare(stats.iter$m.ffpo, stats.iter$m.cgpo)
  cgpo.vs.cg.m[i] <- ParamCompare(stats.iter$m.cgpo, stats.iter$m.cg)
  cg.vs.ffpo.m[i] <- ParamCompare(stats.iter$m.cg, stats.iter$m.ffpo)
}
non.na.ffpo.vs.cgpo.k <- which(!is.na(ffpo.vs.cgpo.k), arr.ind=TRUE)
non.na.cgpo.vs.cg.k <- which(!is.na(cgpo.vs.cg.k), arr.ind=TRUE)
non.na.cg.vs.ffpo.k <- which(!is.na(cg.vs.ffpo.k), arr.ind=TRUE)
non.na.ffpo.vs.cgpo.m <- which(!is.na(ffpo.vs.cgpo.m), arr.ind=TRUE)
non.na.cgpo.vs.cg.m <- which(!is.na(cgpo.vs.cg.m), arr.ind=TRUE)
non.na.cg.vs.ffpo.m <- which(!is.na(cg.vs.ffpo.m), arr.ind=TRUE)

quasi.pwr.k.path <- paste("./plots/quasi_power_k-",domain_name,".pdf",sep="")
quasi.pwr.m.path <- paste("./plots/quasi_power_m-",domain_name,".pdf",sep="")

# scale parameter
pdf( quasi.pwr.k.path )
plot(non.na.ffpo.vs.cgpo.k ,ffpo.vs.cgpo.k[non.na.ffpo.vs.cgpo.k],xlim=c(0,30), ylim=c(0,1.4),col="black",type="l",
     xlab="problem instance number", ylab="probability", main=paste("bootstrap method",domain_name) )
lines(non.na.cgpo.vs.cg.k, cgpo.vs.cg.k[non.na.cgpo.vs.cg.k], col="seagreen1")
lines(non.na.cg.vs.ffpo.k, cg.vs.ffpo.k[non.na.cg.vs.ffpo.k], col="darkorange2")
legend("topleft",
       col=c("black","seagreen1", "darkorange2"),
       legend=c("P(k_ffpo < k_cgpo)", "P(k_cgpo < k_cg)", "P(k_cg < k_ffpo)"),
       pch=c(15,15,15)
)
dev.off()

# slope parameter
pdf( quasi.pwr.m.path )
plot(non.na.ffpo.vs.cgpo.m ,ffpo.vs.cgpo.m[non.na.ffpo.vs.cgpo.m], ylim=c(0,1.4),col="black",type="l",
     xlab="problem instance number",xlim=c(0,30), ylab="probability", main=paste("bootstrap method",domain_name) )
lines(non.na.cgpo.vs.cg.m, cgpo.vs.cg.m[non.na.cgpo.vs.cg.m], col="seagreen1")
lines(non.na.cg.vs.ffpo.m, cg.vs.ffpo.m[non.na.cg.vs.ffpo.m], col="darkorange2")
legend("topleft",
       col=c("black","seagreen1", "darkorange2"),
       legend=c("P(m_ffpo < m_cgpo)", "P(m_cgpo < m_cg)", "P(m_cg < m_ffpo)"),
       pch=c(15,15,15)
)
dev.off()




# sign test
library(BSDA)
ffpo.vs.cgpo.sign <- rep(NA,30)
cgpo.vs.cg.sign <- rep(NA,30)
cg.vs.ffpo.sign <- rep(NA,30)
for(i in 1:30){
  ffpo.vs.cgpo.sign[i] <- SIGN.test(dat.ffpo$runtime[1:i], dat.cgpo$runtime[1:i])$p.value[1]
  cgpo.vs.cg.sign[i] <- SIGN.test(dat.cgpo$runtime[1:i], dat.cg$runtime[1:i])$p.value[1]
  cg.vs.ffpo.sign[i] <- SIGN.test(dat.cg$runtime[1:i], dat.ffpo$runtime[1:i])$p.value[1]
}
NonNAindex.ffpo.vs.cgpo.sign <- which(!is.na(ffpo.vs.cgpo.sign))
lastNonNA.ffpo.vs.cgpo.sign <- max(NonNAindex.ffpo.vs.cgpo.sign)
NonNAindex.cgpo.vs.cg.sign <- which(!is.na(cgpo.vs.cg.sign))
lastNonNA.cgpo.vs.cg.sign <- max(NonNAindex.cgpo.vs.cg.sign)
NonNAindex.cg.vs.ffpo.sign <- which(!is.na(cg.vs.ffpo.sign))
lastNonNA.cg.vs.ffpo.sign <- max(NonNAindex.cg.vs.ffpo.sign)
max.instance.sign <- max(lastNonNA.ffpo.vs.cgpo.sign, lastNonNA.cgpo.vs.cg.sign, lastNonNA.cg.vs.ffpo.sign)
sign.pval.path <- paste("./plots/sign_p_val-",domain_name,".pdf",sep="")
pdf( sign.pval.path )
plot(1:max.instance.sign,ffpo.vs.cgpo.sign[1:max.instance.sign], ylim=c(0,1.4),col="black",type="l",
     xlab="problem instance number", ylab="p-value", main=paste("sign test",domain_name) )
lines(1:max.instance.sign, cgpo.vs.cg.sign[1:max.instance.sign], col="seagreen1")
lines(1:max.instance.sign, cg.vs.ffpo.sign[1:max.instance.sign], col="darkorange2")
legend("topleft",
       col=c("black","seagreen1", "darkorange2"),
       legend=c("k_ffpo v k_cgpo", "k_cgpo v k_cg", "k_cg v k_ffpo"),
       pch=c(15,15,15)
)
dev.off()

# wilcox test
ffpo.vs.cgpo.wilcox <- rep(NA,30)
cgpo.vs.cg.wilcox <- rep(NA,30)
cg.vs.ffpo.wilcox <- rep(NA,30)
for(i in 1:30){
  ffpo.vs.cgpo.wilcox[i] <- wilcox.test(dat.ffpo$runtime[1:i], dat.cgpo$runtime[1:i])$p.value[1]
  cgpo.vs.cg.wilcox[i] <- wilcox.test(dat.cgpo$runtime[1:i], dat.cg$runtime[1:i])$p.value[1]
  cg.vs.ffpo.wilcox[i] <- wilcox.test(dat.cg$runtime[1:i], dat.ffpo$runtime[1:i])$p.value[1]
}
NonNAindex.ffpo.vs.cgpo.wilcox <- which(!is.na(ffpo.vs.cgpo.wilcox))
lastNonNA.ffpo.vs.cgpo.wilcox <- max(NonNAindex.ffpo.vs.cgpo.wilcox)
NonNAindex.cgpo.vs.cg.wilcox <- which(!is.na(cgpo.vs.cg.wilcox))
lastNonNA.cgpo.vs.cg.wilcox <- max(NonNAindex.cgpo.vs.cg.wilcox)
NonNAindex.cg.vs.ffpo.wilcox <- which(!is.na(cg.vs.ffpo.wilcox))
lastNonNA.cg.vs.ffpo.wilcox <- max(NonNAindex.cg.vs.ffpo.wilcox)
max.instance.wilcox <- max(lastNonNA.ffpo.vs.cgpo.wilcox, lastNonNA.cgpo.vs.cg.wilcox, lastNonNA.cg.vs.ffpo.wilcox)
wilcox.pval.path <- paste("./plots/wilcox_p_val-",domain_name,".pdf",sep="")
pdf( wilcox.pval.path )
plot(1:max.instance.wilcox,ffpo.vs.cgpo.wilcox[1:max.instance.wilcox], ylim=c(0,1.4),col="black",type="l",
     xlab="problem instance number", ylab="p-value", main=paste("Wilcoxon rank-sum",domain_name) )
lines(1:max.instance.wilcox, cgpo.vs.cg.wilcox[1:max.instance.wilcox], col="seagreen1")
lines(1:max.instance.wilcox, cg.vs.ffpo.wilcox[1:max.instance.wilcox], col="darkorange2")
legend("topleft",
       col=c("black","seagreen1", "darkorange2"),
       legend=c("k_ffpo v k_cgpo", "k_cgpo v k_cg", "k_cg v k_ffpo"),
       pch=c(15,15,15)
)
dev.off()

qpdf::pdf_combine(input=c(quasi.pwr.k.path, quasi.pwr.m.path, sign.pval.path, wilcox.pval.path),
                  output=paste("./combined_plots/pwr_analysis-",domain_name,".pdf",sep=""))





# combine all into report-domain_name.pdf

res.path <- paste("./combined_plots/resid_plots-",domain_name,".pdf",sep="")
pwr <- paste("./combined_plots/pwr_analysis-",domain_name,".pdf",sep="")
sims.path <- paste("./combined_plots/simulated-",domain_name,".pdf",sep="")
param.path <- paste("./combined_plots/param_distr-",domain_name,".pdf",sep="")

qpdf::pdf_combine(input=c(param.path, res.path, pwr, actual.runtime.path, sims.path),
                  output=paste("./combined_plots/report-",domain_name,".pdf",sep=""))