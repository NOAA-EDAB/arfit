#' Evaluates the likelihood of AR2
#'
#' Uses full likelihood as explained in Beach and MacKinnon
#' The first two data points are used in the likelihood
#'
#'@param data dd
#'@param rhos Numeric vector (1x2). 1st and 2nd order Autoregressive parameters
#'@param hypothesis Character string. "null" or "alt" determines whether to evaluate
#'under the null or alternative hypothesis. "null" fixes linear component, beta = 0.
#'
#'@return -ve Log likelihood
#'
#'@export

likelihood_ar2 <- function(rhos,data,hypothesis="null") {

  rho1 <- rhos[1]
  rho2 <- rhos[2]
  #print(rhos)

  pass <- isStationary_ar2(rhos)
  if (!pass) {
    return(1e10)
  }

  hypothesis <- tolower(hypothesis)
  nT <- nrow(data)
  x <- data$x
  yt <- data$y

  if (tolower(hypothesis) == "null"){
    # Design matrix (2:nT)
    xt <- cbind(rep(1,nT))
  } else {
    # Design matrix (2:nT)
    xt <- cbind(rep(1,nT),as.vector(x[1:nT]))
  }

  # MLE for beta conditional on rhos
  beta <- est_beta_given_rho2(xt,yt,rho1,rho2)

  # all residuals
  At <- yt - xt%*%beta
  # first two residuals
  A1 <- At[1]
  A2 <- At[2]
  res <- At[3:nT]- rho1*At[2:(nT-1)] - rho2*At[1:(nT-2)]

  # Evaluate the Sum Squares of the Residuals
  SSR1 <- (1-rho2^2)*A1^2
  SSR2 <- -2*rho1*(1+rho2)*A1*A2
  SSR3 <- (1-rho2^2)*A2^2
  SSR4 <- sum(res^2)

  part1 <- -(nT/2)*log(SSR1+SSR2+SSR3+SSR4)
  part2 <- log(1+rho2) + 0.5*log(1-rho1-rho2) + 0.5*log(1+rho1-rho2)

  # likelihood
  like <- part1 + part2

  return(-like)

}
