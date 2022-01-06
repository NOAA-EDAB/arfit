#' Simulation study using grid search
#'
#' Performs a simulation study to assess the performance of the test.
#' Depending on user inputs this could take a while since
#' Rho is estimated using a grid search from -.99:.99:.01
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
#'@param rhoSearch Numeric Vector. Values to use to for rho in grid search eg (seq(-.99,.99,.01))
#'
#'@examples
#'\dontrun{
#'simulation_study_grid(here::here("output.txt"),betaVec = 0.5,rhovec =0.5, sigmaVec = 0.25,
#'nTVec = 10, nSims = 200, nBootSims = 500)
#'}
#'
#'@export

simulation_study_grid_ar1 <- function(outDir=here::here("out.txt"),
                             betaVec = c(.12,.25,.5),
                             rhoVec = c(0, 0.25, 0.5,0.75, 0.95),
                             sigmaVec = c(0.25,0.5,.75),
                             nTVec =  c(10),
                             nSims = 200,
                             nBootSims = 500,
                             rhoSearch=NULL) {
  #[.004,.051,.147] and phi = [0,0.433,0.8] with sigma = .735?



  # betaVec <- c(.12,.25,.5)# c(0,.5,1)
  # rhoVec <- c(0, 0.25, 0.5,0.75, 0.95)
  # sigmaVec <- c(0.25,0.5,.75)
  # nTVec <-  c(10)#c(10,25,50,100)
  # nSims <- 200
  # nBootSims = 500

  # betaVec <- 0
  # rhoVec <- .25
  # sigmaVec <- .25
  # nTVec <- 25
  #
  # betaVec <- .12 # c(0,.5,1)
  # rhoVec <- .25 #c(0.25, 0.5,0.75, 0.95)
  # sigmaVec <- .25#c(0.25,0.5,.75)
  # nTVec <-  10#c(10,25,50,100)

  vecHeader <- c("beta","rho","nT","sigma","pValueChi2","pvalueBoot")

  write(vecHeader,file=outDir,ncolumns=length(vecHeader),append=TRUE)

  for (beta in betaVec) {
    for (rho in rhoVec) {
      for (nT in nTVec){
        for (sigma in sigmaVec) {
        print(c(beta,rho,nT,sigma))
        icount <- 0
        # allocate memory to a bunch of vectors
        LRstat <- vector(mode="numeric",length=nBootSims) # likelihood ratio statistic
        pVal_boot <- vector(mode="numeric",length=nSims) # pvalue for bootstrap
        pValChi2 <- vector(mode="numeric",length=nSims) # pvalue for chi sq

        for (isim in 1:nSims) {
          print(paste0("sample = ", isim," of ",nSims))
          data <- simulate_ar1(alpha=0,beta=beta,sigma,rho,nT)

          null <- fit_ar1_grid(data,rhoSearch,hypothesis="null")
          alt <- fit_ar1_grid(data,rhoSearch,hypothesis="alt")

          LRstat[1] <- -2*(null$likelihood-alt$likelihood)
          # pvalue using chi square approximation
          pValChi2[isim] <- 1-pchisq(LRstat[1],1) # uses distributional theory

          # bootstrapping
          for (iboot in 2:nBootSims) {

            # simulate under Null
            bootdata <- simulate_ar1(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT)
        #    print(data)
        #    print(bootdata)
            # fit under null and alt
            nullBoot <- fit_ar1_grid(bootdata,rhoSearch,hypothesis="null")
            altBoot <- fit_ar1_grid(bootdata,rhoSearch,hypothesis="alt")
            # statisicic
            LRstat[iboot] <- -2*(nullBoot$likelihood-altBoot$likelihood)
          } # end bootstrap
          # now we can calculate the p-value based on the bootstrapping
          pVal_boot[isim] <- sum(LRstat >= LRstat[1])/nBootSims
          print(paste0("pval_chi = ",pValChi2[isim],". pval_boot = ",pVal_boot[isim]))
          print(sum((pVal_boot<=.05) & (pVal_boot > 0)))
        } # end sims


        # summarises the proportion of times the test rejects the null
        pValue <- sum(pVal_boot <= 0.05)/nSims
        pvChi <- sum(pValChi2 <= 0.05) /nSims

        # writes to file
        vec <- c(beta,rho,nT,sigma,pvChi,pValue)
        print(vec)
        write(vec,file=outDir,ncolumns=length(vec),append=TRUE)
        } # sigma
      }# end nT
    } # end rho
  } # end beta

  return(list(data=data,null=null,alt=alt))

}


