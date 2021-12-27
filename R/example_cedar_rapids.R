locationOfTHisFile <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(locationOfTHisFile)

# fits to Hydrological data used in Hamed & Rao - 1998

mainTestRealData <- function() {

  if (!exists("AR1LikeNull_alpha",mode="function")){source("AR1LikeNull_alpha.R")}
  if (!exists("AR1LikeAlt_alpha",mode="function")){source("AR1LikeAlt_alpha.R")}

  # load Data
  dataSet <- read.csv("usgs_05464500_cedar_rapids.csv",header=TRUE)
  ind <- (dataSet$year<=1992) # Hamed et al only used 1992 data
  #ind <- (dataSet$year<=2017) & (dataSet$year>1967) # Hamed et al only used 1992 data

#  ind <- (dataSet$year<=2017)
  dataSet <- dataSet[ind,]
  nT <- dim(dataSet)[1]
  data <- data.frame(x=c(1:nT),y=dataSet$riverflow)

  #png("Figure1.png",width=900,height=600,units="px")
  #par(mai=c(1,1.5,0,0),oma=c(0,0,1,1))
 #   plot(dataSet$year,dataSet$riverflow,type="l",xlab="Year",ylab=expression(Discharge~ (f^3/s)),lty=2)
  #plot(dataSet$year,dataSet$riverflow,type="l",
  #       xlab="Year",ylab=expression(Discharge~ (f^3/s)),
  #       cex.lab=2.5,cex.axis=2,lwd=2)
  #  dev.off()

  nBootSims = 1000

  # fit under the null and alternative
  null <- AR1LikeNull_alpha(data,nT)
  alt <- AR1LikeAlt_alpha(data,nT)

  LRstat <- vector(mode="numeric",length=nBootSims) # likelihood ratio statistic

  LRstat[1] <- -2*(null$likelihood-alt$likelihood)
  print(paste0("LR stat = ",LRstat[1]))
  #return()
  # pvalue using chi quare approcimation
  #pValChi2[isim] <- 1-pchisq(LRstat[1],1) # uses distributional theory

          # bootstrapping
  for (iboot in 2:nBootSims) {
    # simulate under Null
    bootdata <- simulateData(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT)
    # fit under null and alt
    nullBoot <- AR1LikeNull_alpha(bootdata,nT)
    altBoot <- AR1LikeAlt_alpha(bootdata,nT)
    # statisicic
    LRstat[iboot] <- -2*(nullBoot$likelihood-altBoot$likelihood)
  } # end bootstrap

  # now we can calculate the p-value based on the bootstrapping
  pVal_boot <- sum(LRstat >= LRstat[1])/nBootSims
  print(paste0("pval_boot = ",pVal_boot))

  # plot fits under the null and alternative
  #png("Figure7.png",width=900,height=600,units="px")
  par(mai=c(1,1.5,0,0),oma=c(0,0,1,1))
  plot(dataSet$year,dataSet$riverflow,type="l",xlab="Year",ylab=expression(Discharge~ (f^3/s)),
       cex.lab=2.5,cex.axis=2,lwd=2)
  lines(dataSet$year,rep(null$betaEst,nT),col="black",lty=2,lwd=2)
  lines(dataSet$year,alt$betaEst[1]+alt$betaEst[2]*c(1:nT),col="black",lty=3,lwd=2)
  #dev.off()
  return(list(null=null, alt=alt))
}




simulateData <- function(alpha,beta=0,sigma,rho,nT){
  # model->   y_t = b*t + z_t
  # where z_t = rho* z_t-1 + e_t   (e_t ~ N(0,sigma^2))
  xt <- c(1:nT)
  # simulate AR - error process
  zt <- vector(mode = "numeric",length=nT)
  zt[1] <- rnorm(1,mean=0,sd=sqrt((sigma^2)/(1-rho^2)))
  for (it in 2:nT) {
    zt[it] <- rho*zt[it-1]+rnorm(1,0,sigma)
  }

  # simulate y
  yt <- rep(alpha,nT) + xt*beta + zt
  data <- data.frame(x=xt,y=yt)
  return(data)
}
