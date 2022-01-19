#' Simulates n data sets and performs bootstrapping
#'
#' Performs n tests
#'
#' For convenience the intercept, beta_0 is set to zero
#'
#'@param outDir Character string. Path to output file
#'@param betaVec Numeric vector. Values for beta_1 (slope/trend parameter)
#'@param rhoVec Numeric vector. Values for autoregressive parameter
#'@param sigmaVec Numeric vector. Values of standard deviation of noise
#'@param nVec Numeric vector. Values for the length of time series to simulate
#'@param nSims Numeric scalar. Number of time series to simulate
#'@param nBootSims Numeric scalar. Number of bootstrap data sets
#'
#'@examples
#'\dontrun{
#'}
#'
#'
#'@export


sim_study_opt <- function(outDir=here::here("out.txt"),
                          betaVec = c(0,.12,.25,.5),
                          rhoVec = c(0, 0.25, 0.5,0.75, 0.95),
                          sigmaVec = c(0.25,0.5,.75),
                          nTVec =  c(10,25,50),
                          nSims = 200,
                          nBootSims = 500) {


  vecHeader <- c("beta","rho","nT","sigma","pValueChi2","pvalueBoot")

  write(vecHeader,file=outDir,ncolumns=length(vecHeader),append=TRUE)

  ic <- 0
  nn <- length(betaVec)*length(rhoVec)*length(sigmaVec)*length(nTVec)

  return(nn)
  for (beta in betaVec) {
    for (rho in rhoVec) {
      for (nT in nTVec){
        for (sigma in sigmaVec) {
          ic <- ic + 1
          message(paste0("Simulation ",ic," of ", nn))

          print(c(beta,rho,nT,sigma))

          sigStats <- sim_single_opt(beta,rho,sigma,nT,nSims,nBootSims)

          vec <- c(beta,rho,nT,sigma,sigStats$pvChi,sigStats$pValue)
          print(vec)
          write(vec,file=outDir,ncolumns=length(vec),append=TRUE)

        } #sigma
      } #nT
    } #rho
  } #beta




}


