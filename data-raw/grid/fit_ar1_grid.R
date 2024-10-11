#' Fit linear model with AR(1) error by GLS
#'
#'
#'
#'
#' @noRd

fit_ar1_grid <- function(data,rhoVec=NULL,hypothesis) {

  nT <- nrow(data)

  if(tolower(hypothesis) == "null") {
    xt <- cbind(rep(1,nT))
  } else {
    xt <- cbind(rep(1,nT),data$x)
  }
  yt <- data$y
  if (is.null(rhoVec)) {
    rhoVec <- seq(-.99,.99,.01)
  }
  minLike <- Inf
  for (rho in rhoVec) {

    #evaluate the likelihood. Returns -ve likelihood
    like <- likelihood_ar1(rho,data,hypothesis)

    if (like < minLike) {
      minLike <- like
      maxLike <- -minLike
      # MLE for beta conditional on rho
      maxBeta <- est_beta_given_rho(xt,yt,rho)
      maxRho <- rho
      ut <- yt- xt%*%maxBeta
      maxSigma <- sqrt((((1-maxRho^2)*ut[1]^2) + sum((ut[2:nT]-maxRho*ut[1:(nT-1)])^2))/nT)
      xt1 <- xt[2:nT]
      xtm1 <- xt[1:(nT-1)]
      varBeta <- 1/((1/maxSigma^2)*((1-maxRho^2)*xt[1]^2 + sum((xt1-maxRho*xtm1)^2)))
    }

  }

  return(list(betaEst=maxBeta,rhoEst=maxRho,sigmaEst=maxSigma,varBeta=varBeta,likelihood=maxLike))
}

