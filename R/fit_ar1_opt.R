#' Fit linear model with trend using optimization routine
#'
#'
#'
#'
#' @export

fit_ar1_opt <- function(data,rho,hypothesis) {

  nT <- nrow(data)

  if(tolower(hypothesis) == "null") {
    xt <- cbind(rep(1,nT))
  } else {
    xt <- cbind(rep(1,nT),data$x)
  }
  yt <- data$y

  #data <- remove_missing_data(data)

  outList <- optim(par=rho,fn=likelihood_ar1,gr=NULL,data,hypothesis,method="L-BFGS-B",lower=-.99,upper=.99)

  maxRho <- outList$par
  maxLike <- -outList$value

  # identify missing values
  indexOfMissingyt <- which(!is.na(yt))
  nTn <- length(indexOfMissingyt)
  # adjust for missing values
  ytn <- yt[!is.na(yt)]
  xtn <- as.matrix(xt[indexOfMissingyt,])

  # MLE for beta conditional on rho
  maxBeta <- est_beta_given_rho(xtn,ytn,maxRho)

  # calculate MLE for sigma
  ut <- yt- xt%*%maxBeta
  maxSigma <- sqrt((((1-maxRho^2)*ut[1]^2) + sum((ut[2:nT]-maxRho*ut[1:(nT-1)])^2,na.rm = T))/nTn)




  # with missing values this need to be reworked
  if(0) {
    # calculate the variance of Beta
    xt1 <- xt[2:nT]
    xtm1 <- xt[1:(nT-1)]
    varBeta <- 1/((1/maxSigma^2)*((1-maxRho^2)*xt[1]^2 + sum((xt1-maxRho*xtm1)^2)))
  } else {
    varBeta <- NA
  }


  return(list(betaEst=maxBeta,rhoEst=maxRho,sigmaEst=maxSigma,varBeta=varBeta,likelihood=maxLike))

}
