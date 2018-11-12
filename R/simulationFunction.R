#' simulations with two gaussian distributed predictors
#'
#' @param B simulations runs
#' @param n sample size
#' @param beta1 first coefficient
#' @param beta2 second coefficient
#' @param censrate censoring rate for exponential model
#' @return list
#' @examples
#' \dontrun{
#' require("survival")
#' require("lava")
#' # setup
#' simmat  <- expand.grid(n=c(100, 200, 350), beta1=c(0,0.5,1),beta2=c(0,1))
#' # simulations w/o censoring
#' simres <- t(apply(simmat, 1, function(d) unlist(simNonNestedNormal(B=10, n=d[1], beta1=d[2], beta2=d[3], censrate=0))))
#' simres <- cbind(simmat,simres)
#' }
#'
#' @export
#'

simNonNestedNormal <- function(B, n, beta1, beta2, censrate=0) {

   resultsNonN <- resultsN <- resultsNonNlog <- resultsNonNc <- resultsNc <- resultsNonNlogc <- list()
   censprob <- rep(NA,B)

   ### simulation set-up
   m <- lvm()
   # exponentially distributed, no censoring
   distribution(m,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=1)
   distribution(m,"censtime")  <- coxWeibull.lvm(scale=censrate,shape=1)
   m <- eventTime(m,time~min(eventtime=1,censtime=0),"status")
   #  predictors
   distribution(m,"V1") <- gaussian.lvm()
   distribution(m,"V2") <- gaussian.lvm()
   # coefficients
   regression(m,from="V1",to="eventtime") <- beta1
   regression(m,from="V2",to="eventtime") <- beta2

   ### simulation runs
   for (runs in 1:B) {

      # seed
      set.seed(123+runs)
      # simulate data
      simdat       <- sim(m,n)
      simdat$expV1 <- exp(simdat$V1)
      # fit competing models
      m1 <- coxph(Surv(time, status) ~ V1, simdat,x=T)
      m2 <- coxph(Surv(time, status) ~ V2, simdat,x=T)
      m3 <- coxph(Surv(time, status) ~ V1    + V2, simdat,x=T)
      m4 <- coxph(Surv(time, status) ~ expV1 + V2, simdat,x=T)

      # tests
      resultsN[[runs]]       <- unlist(plrtest(m1, m3, nested=T)[1:7])
      resultsNonN[[runs]]    <- unlist(plrtest(m1, m2, nested=F)[1:7])
      resultsNonNlog[[runs]] <- unlist(plrtest(m4, m3, nested=F)[1:7])
      # cindex
      resultsNc[[runs]]       <- .conctest(m1, m3)["pvalue"]
      resultsNonNc[[runs]]    <- .conctest(m1, m2)["pvalue"]
      resultsNonNlogc[[runs]] <- .conctest(m4, m3)["pvalue"]

      # censoring
      censprob[runs]         <- sum(1 - simdat$status)/n
   }

   ### test results
   # LRT
   resNonN    <- data.frame(do.call("rbind",resultsNonN))
   resNonNlog <- data.frame(do.call("rbind",resultsNonNlog))
   resN       <- data.frame(do.call("rbind",resultsN))

   # average censoring
   resCens <- mean(censprob, na.rm=T)

   ### power
   # c-index
   cN       <- .powerfunc(unlist(resultsNc))
   cNonN    <- .powerfunc(unlist(resultsNonNc))
   cNonNlog <- .powerfunc(unlist(resultsNonNlogc))

   # nested
   H0r <- .powerfunc(resN$pLRTAB)
   H0  <- .powerfunc(resN$pLRT)
   resnest <- c("H0 robust"=H0r, "H0 classical"=H0, "cindex"=cN)

   # non-nested
   H1    <- .powerfunc(resNonN$pOmega1)
   H1b   <- .powerfunc(resNonN$pOmega2)
   H0    <- .powerfunc(resNonN$pLRTAB)
   seqH0 <- .powerfunc(pmax(resNonN$pOmega1, resNonN$pLRTAB))
   resnonnest <- c("H1 Fine"=H1,"H1 Vuong"=H1b,"H0"=H0,"H0seq"= seqH0,"cindex"=cNonN)

   H1    <- .powerfunc(resNonNlog$pOmega1)
   H1b   <- .powerfunc(resNonNlog$pOmega2)
   H0    <- .powerfunc(resNonNlog$pLRTAB)
   seqH0 <- .powerfunc(pmax(resNonNlog$pOmega1, resNonNlog$pLRTAB))
   resnonnestlog <- c("H1 Fine"=H1,"H1 Vuong"=H1b,"H0"=H0,"H0seq"= seqH0,"cindex"=cNonNlog)

   return(list(B1vsB1B2=resnest, B1vsB2=resnonnest, logB1B2vsB1B2=resnonnestlog, censProp=resCens))
}

#' simulations with two binomial distributed predictors
#'
#' @param B simulations runs
#' @param n sample size
#' @param beta1 first coefficient
#' @param beta2 second coefficient
#' @param censrate censoring rate for exponential model
#' @return list
#' @examples
#' \dontrun{
#' require("lava")
#' # setup
#' simmat  <- expand.grid(n=c(100, 200, 350), beta1=c(0,0.5,1),beta2=c(0,1))
#' # simulations w/o censoring
#' simres <- t(apply(simmat, 1, function(d) unlist(simNonNestedBinomial(B=10, n=d[1], beta1=d[2], beta2=d[3], censrate=0))))
#' simres <- cbind(simmat,simres)
#' }
#'
#'
#' @export


simNonNestedBinomial <- function(B, n, beta1, beta2, censrate=0) {

  resultsNonN <- resultsN <- resultsNonNc <- resultsNc <- list()
  censprob <- rep(NA,B)

  ### simulation set-up
  m <- lvm()
  # exponentially distributed, no censoring
  distribution(m,"eventtime") <- coxWeibull.lvm(scale=1/100,shape=1)
  distribution(m,"censtime")  <- coxWeibull.lvm(scale=censrate,shape=1)
  m <- eventTime(m,time~min(eventtime=1,censtime=0),"status")
  #  predictors
  distribution(m,"V1") <- binomial.lvm(p=0.5)
  distribution(m,"V2") <- binomial.lvm(p=0.5)
  # coefficients
  regression(m,from="V1",to="eventtime") <- beta1
  regression(m,from="V2",to="eventtime") <- beta2

  ### simulation runs
  for (runs in 1:B) {

    # seed
    set.seed(123+runs)
    # simulate data
    simdat       <- sim(m,n)
    # fit competing models
    m1 <- coxph(Surv(time, status) ~ V1, simdat,x=T)
    m2 <- coxph(Surv(time, status) ~ V2, simdat,x=T)
    m3 <- coxph(Surv(time, status) ~ V1 + V2, simdat,x=T)

    # tests
    resultsN[[runs]]       <- unlist(plrtest(m1, m3, nested=T)[1:7])
    resultsNonN[[runs]]    <- unlist(plrtest(m1, m2, nested=F)[1:7])
    # cindex
    resultsNc[[runs]]       <- .conctest(m1, m3)["pvalue"]
    resultsNonNc[[runs]]    <- .conctest(m1, m2)["pvalue"]

    # censoring
    censprob[runs]         <- sum(1 - simdat$status)/n
  }

  ### test results
  # LRT
  resNonN    <- data.frame(do.call("rbind",resultsNonN))
  resN       <- data.frame(do.call("rbind",resultsN))

  # average censoring
  resCens <- mean(censprob, na.rm=T)

  ### power
  # c-index
  cN       <- .powerfunc(unlist(resultsNc))
  cNonN    <- .powerfunc(unlist(resultsNonNc))

  # nested
  H0r <- .powerfunc(resN$pLRTAB)
  H0  <- .powerfunc(resN$pLRT)
  resnest <- c("H0 robust"=H0r, "H0 classical"=H0, "cindex"=cN)

  # non-nested
  H1    <- .powerfunc(resNonN$pOmega1)
  H1b   <- .powerfunc(resNonN$pOmega2)
  H0    <- .powerfunc(resNonN$pLRTAB)
  seqH0 <- .powerfunc(pmax(resNonN$pOmega1, resNonN$pLRTAB))
  resnonnest <- c("H1 Fine"=H1,"H1 Vuong"=H1b,"H0"=H0,"H0seq"= seqH0, "cindex"=cNonN)

  return(list(B1vsB1B2=resnest, B1vsB2=resnonnest, censProp=resCens))
}

### power function
.powerfunc <- function(x) sum(as.numeric(x[!is.na(x)] < 0.05))/length(x[!is.na(x)])

### test on difference in c-index (see concordance vignette survival package)
.conctest <- function(m1, m2){
  ctest <- survival::concordance(m1, m2)
  contr <- c(-1, 1)
  dtest <- contr %*% coef(ctest)
  dvar  <- contr %*% vcov(ctest) %*% contr
  zstat <- dtest/sqrt(dvar)
  pval  <- 2 * (1 - pnorm(abs(zstat)))
  res   <- c(cindex1=coef(ctest)[1], cindex2=coef(ctest)[2], cindexDiff=dtest, sd=sqrt(dvar), z=zstat, pvalue=pval)
  return(res)
}
