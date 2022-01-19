#' Simulate data from an linear process with AR2 error
#'
#'#' model=   y_t = a + b.t + z_t
#' where z_t = rho1. z_t-1 + rho2. z_t-2 + e_t   (e_t ~ N(0,sigma^2))
#' Uses unconditional mean and variance to simulate first two data points
#'
#'@param alpha Numeric scalar. Intercept
#'@param beta Numeric Scalar. Slope
#'@param sigma Numeric scalar. Variance of error process
#'@param rhos Numeric vector (1x2). AR2 parameters
#'@param nT Numeric scalar. Length of time series
#'
#'
#'
#' @export


simulate_ar2 <- function(alpha,beta=0,sigma,rhos,nT){
  # model->   y_t = a + b*t + z_t
  # where z_t = rho* z_t-1 + e_t   (e_t ~ N(0,sigma^2))
  xt <- c(1:nT)
  rho1 <- rhos[1]
  rho2 <- rhos[2]

  # stationarity check
  if(!isStationary_ar2(rhos)) {
    stop(paste0("AR2 process NOT stationary. rho1 = ",rho1," and rho2 = ",rho2))
  }

  # simulate AR - error process
  zt <- vector(mode = "numeric",length=nT)
  yt <- vector(mode = "numeric",length=nT)

  # unconditional mean of process
  meanAR2 <-  0
  # unconditional variance of the process
  varAR2 <- ((1-rho2)*sigma^2)/((1+rho2)*(((1-rho2)^2)-rho1^2))
  zt[1] <- rnorm(1,mean=meanAR2,sd=sqrt(varAR2))
  zt[2] <- rnorm(1,mean=meanAR2,sd=sqrt(varAR2))
  for (it in 3:nT) {
    zt[it] <- rho1*zt[it-1] + rho2*zt[it-2] + rnorm(1,0,sigma)
  }

  # simulate y
  yt <- as.vector(alpha) + xt*as.vector(beta) + zt
  data <- data.frame(x=xt,y=yt)
  return(data)
}
