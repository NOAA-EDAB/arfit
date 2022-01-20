#' Simulation study using optimization routine, multiple cores
#'
#' Performs a simulation study to assess the performance of the test.
#' Utilizes multiple cores to spread bootstrap samples over multiple cores
#'
#' For convenience the intercept, beta_0 is set to zero
#'
#'@param outDir Character string. Path to output file
#'@param betaVec Numeric vector. Values for beta_1 (slope/trend parameter)
#'@param rhoVec1 Numeric vector. Values for autoregressive parameter order 1
#'@param rhoVec2 Numeric vector. Values for autoregressive parameter order 2
#'@param sigmaVec Numeric vector. Values of standard deviation of noise
#'@param nVec Numeric vector. Values for the length of time series to simulate
#'@param nSims Numeric scalar. Number of time series to simulate
#'@param nBootSims Numeric scalar. Number of bootstrap data sets
#'@param setSeed Numeric scalar. Value of the seed for simulations. (Default = NULL, a random number between 1-e7 is selected)
#'
#'@examples
#'\dontrun{
#'}
#'
#'
#'@export


sim_study_opt_ar2 <- function(outDir=here::here("out.txt"),
                          betaVec = c(0,.12,.25,.5),
                          rhoVec1 = c(0, 0.25, 0.5,0.75, 0.95),
                          rhoVec2 = c(0, 0.25, 0.5,0.75, 0.95),
                          sigmaVec = c(0.25,0.5,.75),
                          nTVec =  c(10,25,50),
                          nSims = 200,
                          nBootSims = 500,
                          setSeed=NULL) {

  if(is.null(setSeed)) {
    setSeed <- sample(1e7,1)
  }

  nC <- parallel::detectCores()
  cl <- parallel::makeCluster(nC-1)
  doParallel::registerDoParallel(cl)
  doRNG::registerDoRNG(seed = setSeed)
  starttime <- Sys.time()

  vecHeader <- c("beta","rho1","rho2","nT","sigma","pValueChi2","pvalueBoot")

  write(vecHeader,file=outDir,ncolumns=length(vecHeader),append=TRUE)

  ic <- 0
  nn <- length(betaVec)*length(rhoVec1)*length(sigmaVec)*length(nTVec)

  for (beta in betaVec) {
    irho <- 0
    for (rho1 in rhoVec1) {
      irho <- irho+1
      rho2 <- rhoVec2[irho]
      isStationary <- isStationary_ar2(c(rho1,rho2))
      if (!isStationary) next

      for (nT in nTVec){
        for (sigma in sigmaVec) {
          tictoc::tic()
          ic <- ic + 1
          message(paste0("Simulation ",ic," of ", nn))

          print(c(beta,rho1,rho2,nT,sigma))

          # bootstrap samples simulated in parallel
          sigStats <- sim_single_opt_ar2(beta,rho1,rho2,sigma,nT,nSims,nBootSims)

          vec <- c(beta,rho1,rho2,nT,sigma,sigStats$pvChi,sigStats$pValue)
          print(vec)
          write(vec,file=outDir,ncolumns=length(vec),append=TRUE)
          tictoc::toc()

        } #sigma
      } #nT
    } #rho
  } #beta
  endtime <- Sys.time()

  message(paste0("Elapsed time = ",endtime-starttime))


  parallel::stopCluster(cl)

}


