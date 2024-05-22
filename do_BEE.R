library(tidyverse)

wd = "-- your wd here --"

# function that performs BEE
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