#' Simulate data from an linear process with AR1 error
#'
#' model=   y_t = a + b*t + z_t
#' where z_t = rho* z_t-1 + e_t   (e_t ~ N(0,sigma^2))
#'
#'@param alpha Numeric scalar. Intercept
#'@param beta Numeric scalar. Slope
#'@param rho Numeric scalar. Auto regressive parameter
#'@param sigma Numeric scalar. Standard deviation of error term
#'@param n Numeric scalar. Length of time series
#'
#'
#' @export

simulate_ar <- function(alpha,beta=0,sigma,rho,nT){
  #
  #
  xt <- c(1:n)
  # simulate AR - error process
  zt <- vector(mode = "numeric",length=n)
  zt[1] <- rnorm(1,mean=0,sd=sqrt((sigma^2)/(1-rho^2)))
  for (it in 2:n) {
    zt[it] <- rho*zt[it-1]+rnorm(1,0,sigma)
  }

  # simulate y
  yt <- rep(alpha,n) + xt*beta + zt
  data <- data.frame(x=xt,y=yt)
  return(data)
}
