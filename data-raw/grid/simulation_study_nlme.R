#' Simulation study using nlme package for GLS
#'
#' Uses conditional likelihood. This is not equivalent since it omits first observation
#' when fitting AR1
#'
#'@param outDir Character string. Path to output file
#'
#'@noRd


simulation_study_nlme <- function(outDir=here::here("nlme.txt")) {
  set.seed(143)

  # setwd("C:/Users/Andrew.Beet/Documents/MyWork/paperIdeas/trendAR")
  #[.004,.051,.147] and phi = [0,0.433,0.8] with sigma = .735?

  # Parametric Bootstrap to assess the significance of trend in the presence of autocorrelation
  # uses nlme::gls to do fitting by maximum likelihood.
  # model = a + b.x + error
  #                   error ~ AR(1)     - see simulateData function below

  #alpha <- 0
  #betaVec <- c(0.00,.004,.051,.147,.5)
  #rhoVec <- c(0.00,.43,.8)
  #sigma <- .735
  #nTVec <- c(10,20,30)
  nSims <- 1
  nBootSims = 500

  # betaVec <- c(0)
  # rhoVec <- c(0.25, 0.5,0.75, 0.95)
  # sigmaVec <- c(0.25,0.5,.75)
  # nTVec <-  c(25,50)

  betaVec <- 0
  rhoVec <- 0.25
  sigmaVec <- c(.25)
  nTVec <-c(25)

  vecHeader <- c("beta","rho","nT","sigma","pValueChi2","pvalueBoot")
  write(vecHeader,file=outDir,ncolumns=length(vecHeader),append=TRUE)

  for (beta in betaVec) {
    for (sigma in sigmaVec) {
      for (nT in nTVec) {
        for  (rho in rhoVec) {
        print(c(beta,rho,sigma,nT))

        # allocate memory to a bunch of vectors

        pVal_boot <- vector(mode="numeric",length=nSims) # pvalue for bootstrap
        pValChi2 <- vector(mode="numeric",length=nSims) # pvalue for chi sq

        # simulate nSims data sets
        for (isim in 1:nSims) {
          LRstat <- vector(mode="numeric",length=nBootSims) # likelihood ratio statistic
          print(paste0("sample = ", isim," of ",nSims))
          # simulate data using the current values of beta, rho, sigma, n in loop
          # Sean - use your own simulation routine
          data <- simulate_ar1(alpha=0,beta=beta,sigma,rho,nT)
          # fit under the null (beta = 0) and the alternative (estimate beta)
          glsObjNull <- nlme::gls(y ~ 1,data=data, correlation = nlme::corAR1(form= ~x),method="ML")
          glsObjAlt <- nlme::gls(y ~ x,data=data, correlation = nlme::corAR1(form= ~x),method="ML")

          # likelihood ratio statistic
          LRstat[1] <- -2*(glsObjNull$logLik- glsObjAlt$logLik)
          # pvalue using chi quare approcimation
          pValChi2[isim] <- 1-pchisq(LRstat[1],1) # uses distributional theory

          # these are estimates of a, sigma and rho under the null. We need these to do the bootstrap
          alphaEst <- glsObjNull$coefficients[1]
          sigmaEst <- glsObjNull$sigma
          rhoEst <- coef(glsObjNull$modelStruct$corStruct,unconstrained=FALSE)

          print(c(alphaEst,sigmaEst,rhoEst))
       #   return()
          # bootstrapping. simulate and fit nBootSims samples
          for (iboot in 2:nBootSims) {

            # simulate under Null
            # you already have your own simulation method, use that
            # uses the parameter estimates above - fitting under the null
            bootdata <- simulate_ar1(alpha=alphaEst,beta=0,sigmaEst,rhoEst,nT)
            # fit under null and alt
            # i think you call y = series and x = time. You'll need to change that
            # i also simulate with alpha (the intercept) = 0, i dont think you do.
            #In that case you will remove the "-1"
            # from the calls to gls
            glsObjBootNull <- nlme::gls(y ~ 1,data=bootdata, correlation = nlme::corAR1(form= ~x),method="ML")
            glsObjBootAlt <- nlme::gls(y ~ x,data=bootdata, correlation = nlme::corAR1(form= ~x),method="ML")

            # calculate and store the likelihood ratio statistic for each bootstrapped sample
            LRstat[iboot] <- -2*(glsObjBootNull$logLik- glsObjBootAlt$logLik)
          } # end bootstrap loop


          # now we can calculate the p-value based on the bootstrapping
          pVal_boot[isim] <- sum(LRstat >= LRstat[1])/nBootSims
          print(paste0("pval_chi = ",pValChi2[isim],". pval_boot = ",pVal_boot[isim]))
          print(sum((pVal_boot<=.05) & (pVal_boot > 0)))


          if ((pValChi2[isim] < .05) & (pVal_boot[isim] > 0.05) ) {
            plot(data$x,data$y,type="l",xlab="time",ylab="Indicator")#,ylim=c(-.3,.3)

            return(data)
          }


        } # end sims loop

        # summarises the proportion of times the test rejects the null
        pValue <- sum(pVal_boot <= 0.05)/nSims
        pvChi <- sum(pValChi2 <= 0.05) /nSims
        # writes to file
        vec <- c(beta,rho,nT,sigma,pvChi,pValue)
        print(vec)
        write(vec,file=outDir,ncolumns=length(vec),append=TRUE)
        }
       }
    }
  }

return(list(data=data,glsObjNull=glsObjNull,glsObjAlt=glsObjAlt))


}

