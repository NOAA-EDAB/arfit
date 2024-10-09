#' Fit linear model with trend using optimization routine
#'
#'
#'
#'
#' @export

fit_ar2_opt <- function(data,rho,hypothesis) {

  nT <- nrow(data)

  if(tolower(hypothesis) == "null") {
    xt <- cbind(rep(1,nT))
  } else {
    xt <- cbind(rep(1,nT),data$x)
  }
  yt <- data$y


  outList <- optim(par=rho,fn=likelihood_ar2,gr=NULL,
                   data,hypothesis,method="L-BFGS-B",
                   lower=c(-.99,-.99),upper=c(.99,.99))

  maxRho <- outList$par
  maxLike <- -outList$value

  # MLE for beta conditional on rho
  maxBeta <- est_beta_given_rho2(xt,yt,maxRho[1],maxRho[2])

  # calculate MLE for sigma
  rho1 <- maxRho[1]
  rho2 <- maxRho[2]
  ut <- yt- xt%*%maxBeta
  p1 <- (1-rho2^2)*(ut[1]^2) - 2*rho1*(1+rho2)*ut[1]*ut[2] + ((1-rho2^2)*ut[2]^2)
  p2 <-  sum((ut[3:nT]-rho1*ut[2:(nT-1)] - rho2*ut[1:(nT-2)])^2)
  maxSigma <- sqrt(sum(p1+p2)/nT)


  # calculate the variance of Beta
  varBeta <- NA


  return(list(betaEst=maxBeta,rhoEst=maxRho,sigmaEst=maxSigma,varBeta=varBeta,likelihood=maxLike))

}
