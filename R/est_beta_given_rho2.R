#' Estimate beta conditional on known AR1 parameter
#'
#' Uses Prais and Winsten transformation which utilizes first observation
#'
#'@param xt Numeric vector. Time series of time points (1, ..., n)
#'@param yt Numeric vector. Time series of observations (length = n)
#'@param rho1 Numeric scalar. AR1 parameter
#'@param rho2 Numeric scalar. AR2 parameter
#'
#'@export


est_beta_given_rho2 <- function(xt,yt,rho1,rho2=0) {

  nT <- nrow(xt)

  # preallocate Q
  Q <- matrix(0,nT,nT)
  # Prais & Winsten transformation
  diag(Q) <- 1
  index <- row(Q)-col(Q)
  Q[index==1] <- -rho1
  Q[index==2] <- -rho2

  a2 <- (1-rho2^2)-(rho1^2)*(1+rho2)/(1-rho2)
  b2 <- (rho1^2)*(1+rho2)/(1-rho2)
  c2 <- 1-rho2^2
  bc <- -rho1*(1-rho2)

  Q[1,1] <- sqrt(a2)
  Q[2,1] <- sqrt(b2)*-1*sign(rho1) # need to check this
  Q[2,2] <- sqrt(c2)

  # design matrix and response
  X <- Q%*%xt
  Y <- Q%*%yt

  # mle of beta given rho
  betaEst <- solve(t(X)%*%X)%*% t(X)%*%Y

  return(betaEst)
}

