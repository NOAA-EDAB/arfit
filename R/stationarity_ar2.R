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

# checks to see if AR2 parameters satisfy stationarity
stationarity_ar2 <- function(rhos) {
  rho1 <- rhos[1]
  rho2 <- rhos[2]
  if (((rho1 + rho2) < 1) & ((rho1-rho2)> -1) & (rho2>-1)) {
    return(1)
  } else {
    print(paste0("process not stationary with phi = ",rhos))
    return(0)
  }
}
