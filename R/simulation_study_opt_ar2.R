#' Simulation study using optimization routine
#'
#' Performs a simulation study to assess the performance of the test.
#'
#' For convenience the intercept, beta_0 is set to zero
#'
#'@param outDir Character string. Path to output file
#'@param betaVec Numeric vector. Values for beta_1 (slope/trend parameter)
#'@param rho1Vec Numeric vector. Values for autoregressive parameter, rho_!
#'@param rho2Vec Numeric vector. Values for autoregressive parameter, rho_2
#'@param sigmaVec Numeric vector. Values of standard deviation of noise
#'@param nVec Numeric vector. Values for the length of time series to simulate
#'@param nSims Numeric scalar. Number of time series to simulate
#'@param nBootSims Numeric scalar. Number of bootstrap data sets
#'
#'@examples
#'\dontrun{
#'simulation_study_opt(here::here("output.txt"),betaVec = 0.5,rhovec =0.5, sigmaVec = 0.25,
#'nTVec = 10, nSims = 200, nBootSims = 500)
#'}
#'
#'@export


simulation_study_opt_ar2 <- function(outDir=here::here("out_ar2.txt"),
                                 betaVec = c(0),
                                 rho1Vec = c(0),
                                 rho2Vec = c(0),
                                 sigmaVec = c(0.25),
                                 nTVec =  c(10),
                                 nSims = 200,
                                 nBootSims = 500,
                                 rhoSearch = NULL) {

  #
  # betaVec <- .12 # c(0,.5,1)
  # rhoVec <- .25 #c(0.25, 0.5,0.75, 0.95)
  # sigmaVec <- .25#c(0.25,0.5,.75)
  # nTVec <-  10#c(10,25,50,100)

  vecHeader <- c("beta","rho1","rho2","nT","sigma","pValueChi2","pvalueBoot")

  write(vecHeader,file=outDir,ncolumns=length(vecHeader),append=TRUE)

  for (beta in betaVec) {
    irho <- 0
    for (rho1 in rho1Vec) { # pairs of rho
      irho <- irho + 1
      rho2 <- rho2Vec[irho]
      for (nT in nTVec){
        for (sigma in sigmaVec) {
          print(c(beta,rho1,rho2,nT,sigma))
          isStationary <- isStationary_ar2(c(rho1,rho2))
          if (!isStationary) next
          icount <- 0
          # allocate memory to a bunch of vectors
          LRstat <- vector(mode="numeric",length=nBootSims) # likelihood ratio statistic
          pVal_boot <- vector(mode="numeric",length=nSims) # pvalue for bootstrap
          pValChi2 <- vector(mode="numeric",length=nSims) # pvalue for chi sq



          for (isim in 1:nSims) {
            print(paste0("sample = ", isim," of ",nSims))
            data <- simulate_ar2(alpha=0,beta=beta,sigma,c(rho1,rho2),nT)

            null <- fit_ar2_opt(data,rho=c(rho1,rho2),hypothesis="null")
            alt <- fit_ar2_opt(data,rho=c(rho1,rho2),hypothesis="alt")

            LRstat[1] <- -2*(null$likelihood-alt$likelihood)
            # pvalue using chi square approximation
            pValChi2[isim] <- 1-pchisq(LRstat[1],1) # uses distributional theory

            # bootstrapping
            for (iboot in 2:nBootSims) {
              # simulate under Null
#              print(iboot)
              bootdata <- simulate_ar2(alpha=null$betaEst,beta=0,null$sigmaEst,null$rhoEst,nT)

              # fit under null and alt
              nullBoot <- fit_ar2_opt(bootdata,rho=,null$rhoEst,hypothesis="null")
              # print("alt")
              # if(iboot == 152) {
              #   print(nullBoot)
              #   return(bootdata)
              # }
              altBoot <- fit_ar2_opt(bootdata,rho=nullBoot$rhoEst,hypothesis="alt")
              # statistic
              LRstat[iboot] <- -2*(nullBoot$likelihood-altBoot$likelihood)
            } # end bootstrap

            # now we can calculate the p-value based on the bootstrapping
            pVal_boot[isim] <- sum(LRstat >= LRstat[1])/nBootSims
            print(paste0("pval_chi = ",pValChi2[isim],". pval_boot = ",pVal_boot[isim]))

          } # end sims


          # summarises the proportion of times the test rejects the null
          pValue <- sum(pVal_boot <= 0.05)/nSims
          pvChi <- sum(pValChi2 <= 0.05) /nSims

          # writes to file
          vec <- c(beta,rho1,rho2,nT,sigma,pvChi,pValue)
          print(vec)
          write(vec,file=outDir,ncolumns=length(vec),append=TRUE)
        } # sigma
      }# end nT
    } # end rho
  } # end beta

  return(vec)
}


