#' Evaluate the likelihood of linear model with errors from an AR order 1 process
#'
#'
#'
#'@param beta Numeric vector. (px1). Defining alpha and beta.
#'@param rho Numeric scalar. Autoregressive parameter
#'@param dataf Data frame. (nx2). Defining xt and yt values
#'@param hypothesis Character string. "null" or "alt" determines whether to evaluate
#'under the null or alternative hypothesis. "null" fixes linear component, beta = 0.
#'
#'
#'
#' @return
#'
#' @export

likelihood_ar1 <- function(beta,rho,dataf,hypothesis="null") {

  hypothesis <- tolower(hypothesis)

  x <- dataf$x
  y <- dataf$y
  if (hypothesis == "null"){
    x1 <- x[1]
  } else {
    x1 <- cbind(1,x[1])
  }
  y1 <- y[1]
  nT <- length(x)

  # response
  yt <- as.vector(y[2:nT] )
  ytm1 <- as.vector(y[1:(nT-1)])

  if (hypothesis == "null") {
    part1 <- .5*log(1-rho^2)
    part2a <- (1-rho^2)*((y1-beta)^2)
    part2b <- sum((yt- rep(beta,nT-1) - rho*ytm1 + rho*rep(beta,nT-1))^2)
    like <- part1 - (nT/2)*log(part2a + part2b)
  } else {
    # Design matrix (2:nT)
    xt <- cbind(rep(1,nT-1),as.vector(x[2:nT]))
    # Design matrix (1:nT-1)
    xtm1 <- cbind(rep(1,nT-1),as.vector(x[1:(nT-1)]))
    # Likelihood evaluation
    part1 <- .5*log(1-rho^2)
    part2a <- (1-rho^2)*((y1-x1%*%beta)^2)
    part2b <- sum((yt- xt%*%beta - rho*ytm1 + rho*xtm1%*%beta)^2)
    like <- part1 - (nT/2)*log(part2a + part2b)
  }

  return(like)

}
