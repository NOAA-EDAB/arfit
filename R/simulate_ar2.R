#' simulare AR2
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'

simulate_ar2 <- function(alpha,beta=0,sigma,rhos,nT){
  # model->   y_t = a + b*t + z_t
  # where z_t = rho* z_t-1 + e_t   (e_t ~ N(0,sigma^2))
  xt <- c(1:nT)
  rho1 <- rhos[1]
  rho2 <- rhos[2]
  # stationarity check


  # simulate AR - error process
  zt <- vector(mode = "numeric",length=nT)
  meanAR2 <-  0/(1-rho1 -rho2)
  varAR2 <- ((1-rho1)*sigma^2)/((1+rho2)*(((1-rho2)^2)-rho1^2))
  zt[1] <- rnorm(1,mean=meanAR2,sd=sqrt(varAR2))
  zt[2] <- rnorm(1,mean=meanAR2,sd=sqrt(varAR2))
  for (it in 3:nT) {
    zt[it] <- rho1*zt[it-1] + rho2*zt[it-2] + rnorm(1,0,sigma)
  }

  # simulate y
  yt <- alpha + xt*beta + zt
  data <- data.frame(x=xt,y=yt)
  return(data)
}
