#' Partial Likelihood Ratio Test for non-nested \code{coxph} models
#'
#' @author Thomas Hielscher
#' @param object1 \code{coxph} model with \code{x=T}
#' @param object2 \code{coxph} model with \code{x=T}, based on the identical data and order of observations as for \code{object1}
#' @param nested  specify if models are nested, default ist \code{FALSE}
#' @param adjusted  specify if test statistic for non-nested models should be adjusted for different complexities of the models using BIC-type adjustment (p/2)log n - (q/2)log n, default ist \code{FALSE}
#' @return a object of class \code{finetest}
#' @references Fine J.P., Comparing nonnested Cox models, Biometrika (2002), 89, 3, 635-647.
#' @references Vuong Q.H., Likelihood Ratio Tests for Model Selection and Non-Nested Hypotheses, Econometrika (1989), 57, 2, 307-333.
#' @references Merkle E.C. and You D., Testing Nonnested Structural Equation Models, Psychological Methods (2016), 21, 2, 151-163.
#' @references Edgar Merkle and Dongjun You (2018). nonnest2: Tests of Non-Nested Models. R package version 0.5-1.
#' @references P. Duchesne, P. Lafaye de Micheaux, Computing the distribution of quadratic forms: Further comparisons between the Liu-Tang-Zhang approximation and exact methods, Computational Statistics and Data Analysis, Volume 54, (2010), 858-862
#' @references Schwarz G. (1978). Estimating the Dimension of a Model. Annals of Statistics, 6, 461.464.
#' @examples
#' \dontrun{
#' ### example data set from Fine paper, section 5
#' require("survival")
#' pbc  <- subset(pbc, !is.na(trt))
#' mod1 <- coxph(Surv(time, status==2) ~ age, data=pbc, x=T)
#' mod2 <- coxph(Surv(time, status==2) ~ age + albumin + bili + edema + protime, data=pbc,  x=T)
#' mod3 <- coxph(Surv(time, status==2) ~ age + log(albumin) + log(bili) + edema + log(protime), data=pbc, x=T)
#' plrtest(mod3, mod2, nested=F) # non-nested models
#' plrtest(mod3, mod1, nested=T) # nested models
#' }
#'
#' @export
#'

plrtest <- function (object1, object2, nested = FALSE, adjusted=FALSE) {

  if (is.null(object1$x) | is.null(object2$x)) stop("coxph object without x=T option fitted")

  # basic check on data
  if (any(object1$y[,1]!=object2$y[,1]) | any(object1$y[,2]!=object2$y[,2])) stop("models not fitted on the same data or data not in the same order")

  if (nested) {
    if (logLik(object2)[1] > logLik(object1)[1]) {
      tmp     <- object1
      object1 <- object2
      object2 <- tmp
    }
  }

  llA <- llcont(object1)
  llB <- llcont(object2)
  lr  <- sum(llA - llB)

  n   <- length(llA)
  z1  <- object1$x
  z2  <- object2$x
  p1  <- ncol(z1)
  p2  <- ncol(z2)

  # information matrix
  I1    <- chol2inv(chol(n * vcov(object1)))
  I2    <- chol2inv(chol(n * vcov(object2)))
  zero1 <- matrix(0,p1,p2)
  zero2 <- matrix(0,p2,p1)
  # score matrix
  S1  <- matrix(crossprod(sandwich::estfun(object1))/n, nrow(I1), nrow(I1))
  S2  <- matrix(crossprod(sandwich::estfun(object2))/n, nrow(I2), nrow(I2))
  S12 <- crossprod(sandwich::estfun(object1), sandwich::estfun(object2))/n
  S21 <- t(S12)
  ### composite matrices
  # Sigma12
  Sigma12 <- cbind(rbind(S1, S21), rbind(S12, S2))
  # K12
  K12    <- cbind(rbind(I1, zero2), rbind(zero1, I2))
  K12inv <- chol2inv(chol(K12))

  # initialize p-values
  pLRTA <- pLRTB <- pLRTAB <- pLRT <- pOmega1 <- pOmega2 <- NA

  ### theorem 1, Fine
  Itilde    <- diag(c(rep(1,p1),rep(-1,p2)), ncol=p1+p2, nrow=p1+p2)
  A12       <- Itilde %*% Sigma12 %*% K12inv # -W from Vuong paper
  eigenPHI  <- sort(Re(eigen(A12, only.values = TRUE)$values))

  ### theorem 3, Fine
  J        <- (n/(n-1))* cbind(rbind(var(z1), -cov(z2,z1)), rbind(-cov(z1,z2), var(z2)))
  B12      <- J %*% K12inv  %*% Sigma12 %*% K12inv
  eigenPSI <- sort(Re(eigen(B12,only.values = TRUE)$values))

  # test statistic for difference in hazards using linear predictor
  lp1      <- drop(z1 %*% coef(object1))
  lp2      <- drop(z2 %*% coef(object2))
  varlp    <- (n - 1)/n  * var(lp1 - lp2)
  pOmega1  <- pmax(0,CompQuadForm::imhof(n*varlp, eigenPSI)$Qq)

  ### test statistic for difference in hazards using Vuong approach based on LLi
  varll    <- (n - 1)/n * var(llA - llB)
  pOmega2  <- pmax(0,CompQuadForm::imhof(n*varll, eigenPHI^2)$Qq)

  if (nested) {
     teststat  <- 2*lr
     pLRTAB    <- pmax(0,CompQuadForm::imhof(teststat, eigenPHI)$Qq)
     pLRT      <- anova(object1, object2)[2,4]
  }

  ### theorem 2, Fine
  if (!nested) {
     if (adjusted) {
       lr <- lr - (p1/2*log(n) - p2/2*log(n))
     }
     teststat   <- (1/sqrt(n)) * lr/sqrt(varll)
     pLRTA      <- pnorm(teststat, lower.tail = FALSE)
     pLRTB      <- pnorm(teststat)
     pLRTAB     <- 2 * min(pLRTA, pLRTB) # two-sided
  }

  rval <- list(pOmega1 = pOmega1, pOmega2 = pOmega2, LRTstat = teststat, pLRT=pLRT, pLRTA= pLRTA, pLRTB = pLRTB, pLRTAB=pLRTAB ,nested = nested)
  class(rval) <- "finetest"
  return(rval)
}

