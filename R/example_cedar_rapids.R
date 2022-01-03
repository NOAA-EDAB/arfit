#' Run test on example data set
#'
#' Use data from USGS 05464500 Cedar River at Cedar Rapids, IA
#'
#' @param dataSet Data frame. Two columns (year, riverflow)
#' @param nBootSims Numeric scalar. Number of bootstrap samples to perform
#'
#'
#'@export

example_cedar_rapids <- function(dataSet=arfit::cedar_rapids,nBootSims=999) {

  ind <- (dataSet$year<=1992) # Hamed et al only used 1992 data

  dataSet <- dataSet[ind,]
  nT <- dim(dataSet)[1]
  data <- data.frame(x=c(1:nT),y=dataSet$riverflow)

  # fit under the null and alternative
  null <- fit_ar1_grid(data,hypothesis ="null")
  alt <- fit_ar1_grid(data,hypothesis="alt")
  # null <- fit_ar1_opt(data,rho=0.5,hypothesis ="null")
  # alt <- fit_ar1_opt(data,null$rhoEst,hypothesis="alt")
  # preallocate likelihood ratio statistic vector
  LRstat <- vector(mode="numeric",length=nBootSims+1)
  # LR stat for data
  LRstat[1] <- -2*(null$likelihood-alt$likelihood)
  print(paste0("LR stat = ",LRstat[1]))
  # pvalue using chi square approximation
  pValChi2 <- 1-pchisq(LRstat[1],1) # uses distributional theory

  # Perform bootstrapping
  for (iboot in 2:(nBootSims+1)) {
    # simulate under Null
    bootdata <- simulate_ar1(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT)
    # fit under null and alt
    nullBoot <- fit_ar1_grid(bootdata,hypothesis="null")
    altBoot <- fit_ar1_grid(bootdata,hypothesis="alt")
    # nullBoot <- fit_ar1_opt(bootdata,null$rhoEst,hypothesis="null")
    # altBoot <- fit_ar1_opt(bootdata,null$rhoEst,hypothesis="alt")
    # statisicic
    LRstat[iboot] <- -2*(nullBoot$likelihood-altBoot$likelihood)
  } # end bootstrap

  # now we can calculate the p-value based on the bootstrapping
  pVal_boot <- sum(LRstat >= LRstat[1])/(nBootSims+1)
  print(paste0("pval_boot = ",pVal_boot))

  # plot fits under the null and alternative
  #png("Figure7.png",width=900,height=600,units="px")
  par(mai=c(1,1.5,0,0),oma=c(0,0,1,1))
  plot(dataSet$year,dataSet$riverflow,type="l",xlab="Year",ylab=expression(Discharge~ (f^3/s)),
       cex.lab=2.5,cex.axis=2,lwd=2)
  lines(dataSet$year,rep(null$betaEst,nT),col="black",lty=2,lwd=2)
  lines(dataSet$year,alt$betaEst[1]+alt$betaEst[2]*c(1:nT),col="black",lty=3,lwd=2)
  #dev.off()
  return(list(null=null, alt=alt,pValChi2=pValChi2))
}
