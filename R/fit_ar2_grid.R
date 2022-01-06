#' Fit linear model with AR(1) error by GLS
#'
#'
#'
#'
#' @export

fit_ar2_grid <- function(data,rhoVec=NULL,hypothesis) {


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
  for (rho1 in rhoVec) {
    for (rho2 in rhoVec) {
      # check for stationarity, pass if not
      pass <- isStationary_ar2(c(rho1,rho2))
      if (!pass) next #not stationary

      #evaluate the likelihood. Returns -ve likelihood
      like <- likelihood_ar2(rhos=c(rho1,rho2),data,hypothesis)
      # print(c(rho1,rho2))
      # print(like)

      if (like < minLike) {
        minLike <- like
        # MLE for beta conditional on rho
        maxBeta <- est_beta_given_rho2(xt,yt,rho1,rho2)
        maxRho <- c(rho1,rho2)
        ut <- yt- xt%*%maxBeta


        p1 <- (1-rho1^2)*(ut[1]^2) - 2*rho1*(1+rho2)*ut[1]*ut[2] + ((1-rho2^2)*ut[2]^2)
        p2 <-  sum((ut[3:nT]-rho1*ut[2:(nT-1)] - rho2*ut[1:(nT-2)])^2)
        maxSigma <- sqrt(sum(p1+p2)/nT)
        varBeta <- NA

      }
    }

  }

  return(list(betaEst=maxBeta,rhoEst=maxRho,sigmaEst=maxSigma,varBeta=varBeta,likelihood=-minLike))
}

