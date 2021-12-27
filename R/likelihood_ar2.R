#' Evaluates the likelihood of AR2
#'
#' Uses full likelihood as explained in Beach and MacKinnon
#' The first two data points are used in the likelihood 
#'
#'@param data
#'@param beta Numeric vector. intercept and slope term
#'@param rho1
#'@param rho2
#'
#'

likelihood_ar2 <- function(data,beta,rho1,rho2) {
  
  T <- nrow(data)
  xt <- data$x
  yt <- data$y
  
  # first two residuals
  A1 <- yt[1] - beta[1] - xt[1]*beta[2]
  A2 <- yt[2] - beta[1] - xt[2]*beta[2]
  # all residuals
  At <- yt - beta[1] - xt*beta[2]
  res <- At[3:T]- rho1*At[2:(T-1)] - rho2*At[1:(T-2)]
  
  # Evaluate the Sum Squares of the Residuals
  SSR1 <- (1-rho1^2)*A1^2
  SSR2 <- -2*rho1*(1+rho2)*A1*A2
  SSR3 <- (1-rho2^2)*A2^2
  SSR4 <- sum(res^2)
  
  part1 <- -(T/2)*log(SSR1+SSR2+SSR3+SSR4)
  part2 <- log(1+rho2) + 0.5*log(1-rho1-rho2) + 0.5*log(1+rho1-rho2)
  
  # likelihood
  like <- part1 + part2
  
  return(like)
  
}