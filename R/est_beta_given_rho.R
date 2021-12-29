#' Estimate beta conditional on known AR1 parameter
#'
#' Uses Prais and Winsten transformation which utilizes first observation
#'
#'


est_beta_given_rho <- function(xt,yt,rho,hypothesis) {

  nT <- nrow(xt)

  # preallocate Q
  Q <- matrix(0,nT,nT)
  # Prais & Winsten transformation
  diag(Q) <- 1
  Q[1,1] <- sqrt(1-rho^2)
  index <- row(Q)-col(Q)
  Q[index==1] <- -rho

  # design matrix and response
  X <- Q%*%xt
  Y <- Q%*%yt

  # mle of beta given rho
  betaEst <- solve(t(X)%*%X)%*% t(X)%*%Y

  return(betaEst)
}

