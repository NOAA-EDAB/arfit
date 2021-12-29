#' Fit linear model with AR(1) error by GLS
#'
#'
#'
#'
#' @export

fit_ar1 <- function(data,nT,hypothesis) {

  if(tolower(hypothesis) == "null") {
    xt <- cbind(rep(1,nT))
  } else {
    xt <- cbind(rep(1,nT),data$x)
  }
  yt <- data$y
  rhoVec <- seq(-.99,.99,.01)
  maxLike <- -Inf

  for (rho in rhoVec) {
    Q <- matrix(0,nT,nT)
    # Prais & Winsten transformation
    diag(Q) <- 1
    Q[1,1] <- sqrt(1-rho^2)
    index <- row(Q)-col(Q)
    Q[index==1] <- -rho

    X <- Q%*%xt
    Y <- Q%*%yt
    # mle of beta given rho
    betaEst <- solve(t(X)%*%X)%*% t(X)%*%Y
    #evaluate the likelihood
    like <- likelihood_ar1(betaEst,rho,data,hypothesis)
   # likeKeep[ik] <- like
    if (like > maxLike) {
      maxLike <- like
      maxBeta <- betaEst
      maxRho <- rho
      ut <- yt- xt%*%maxBeta
      maxSigma <- sqrt((((1-maxRho^2)*ut[1]^2) + sum((ut[2:nT]-maxRho*ut[1:(nT-1)])^2))/nT)
      xt1 <- xt[2:nT]
      xtm1 <- xt[1:(nT-1)]
      varBeta <- 1/((1/maxSigma^2)*((1-maxRho^2)*xt[1]^2 + sum((xt1-maxRho*xtm1)^2)))
    }

  }
  #plot(rhoVec,likeKeep)
  return(list(betaEst=maxBeta,rhoEst=maxRho,sigmaEst=maxSigma,varBeta=varBeta,likelihood=maxLike))
}

