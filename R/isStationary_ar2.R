#' Determines stationarity of AR2
#'
#' Given parameters of AR2 determines if process is stationary
#'
#' @param rhos Numeric vector. AR2 parameters
#'
#' @return Boolean. Process is stationary or not
#'
#'
#'
#' @export

# checks to see if AR2 parameters satisfy stationarity
isStationary_ar2 <- function(rhos) {

  pass <- FALSE
  rho1 <- rhos[1]
  rho2 <- rhos[2]

  if (((rho1 + rho2) < 1) & ((rho1-rho2)> -1) & (rho2>-1)) {
    pass <- TRUE
  }
  return(pass)

}
