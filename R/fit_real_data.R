#' Run test on real data set
#'
#' Performs a parametric bootstrap to test the significance of a linear trend in time under
#' the assumption of autoregressive order 1 error.
#'
#' @param dataSet Data frame.
#' Two of the fields must be named x (time, equally spaced), and y (response variable)
#' @param nBootSims Numeric scalar. Number of bootstrap samples to perform
#' @param printFig Boolean. Print data an fit in figure window (Default = F)
#'
#' @section: ecodata
#'
#'This function is used in ecodata::geom_lm()
#'
#'@export

fit_real_data <- function(dataSet,nBootSims=499,printFig=F) {

  dataValidation <- check_data_validation(dataSet)
  dataSet <- dataValidation$dataSet
  missingValues <- dataValidation$missingValues

  data <- dataSet
  nT <- nrow(data)
  # fit under the null and alternative
  null <- fit_ar1_opt(data,rho = 0,hypothesis ="null")
  alt <- fit_ar1_opt(data,rho = 0,hypothesis="alt")

  # preallocate likelihood ratio statistic vector
  LRstat <- vector(mode="numeric",length=nBootSims+1)
  # LR stat for data
  LRstat[1] <- -2*(null$likelihood-alt$likelihood)
  #print(paste0("LR stat = ",LRstat[1]))
  # pvalue using chi square approximation
  pValChi2 <- 1-pchisq(LRstat[1],1) # uses distributional theory

  # Perform bootstrapping
  for (iboot in 2:(nBootSims+1)) {
    # simulate under Null
    bootdata <- simulate_ar1(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT,missingValues = missingValues)

    dataValidation <- check_data_validation(bootdata)
    bootdata <- dataValidation$dataSet

    # fit under null and alt
    nullBoot <- fit_ar1_opt(bootdata,null$rhoEst,hypothesis="null")
    altBoot <- fit_ar1_opt(bootdata,null$rhoEst,hypothesis="alt")

    # statisicic
    LRstat[iboot] <- -2*(nullBoot$likelihood-altBoot$likelihood)
  } # end bootstrap

  # now we can calculate the p-value based on the bootstrapping
  pVal_boot <- sum(LRstat >= LRstat[1])/(nBootSims+1)


  if(printFig) {
    print(paste0("pval_boot = ",pVal_boot))
    par(mai=c(1,1.5,0,0),oma=c(0,0,1,1))
    plot(dataSet$x,dataSet$y,type="l",xlab="Year",ylab="Response",
         cex.lab=2.5,cex.axis=2,lwd=2)
    lines(dataSet$x,rep(null$betaEst,nT),col="black",lty=2,lwd=2)
  #  lines(dataSet$x,alt$betaEst[1]+alt$betaEst[2]*c(1:nT),col="black",lty=3,lwd=2)
    lines(dataSet$x,alt$betaEst[1]+alt$betaEst[2]*dataSet$x,col="black",lty=3,lwd=2)
  }


  return(list(null=null, alt=alt,pValue=pVal_boot))
}
