#' Getting partial log-likelihood of a \code{coxph} object for individual cases
#'
#' @author Thomas Hielscher
#' @param x \code{coxph} model fitted with \code{x=T}
#' @details ties are handled according to Efron
#' @return a vector of individual likelihoods (sorted by decreasing event times)
#' @examples
#' \dontrun{
#' ### example data set from Fine paper, section 5
#' require("survival")
#' mod <- coxph(Surv(time, status==2) ~ age + log(albumin) + log(bili) + edema + log(protime), data=subset(pbc, !is.na(trt)),  x=T)
#' # individual LLs sum up to total model LL
#' sum(llcont(mod));logLik(mod)[1]
#' }
#' @importFrom data.table data.table
#' @export

llcont <- function(x) UseMethod("llcont")

#' Getting partial log-likelihood of a \code{coxph} object for individual cases
#'
#' @author Thomas Hielscher
#' @param x \code{coxph} model fitted with \code{x=T}
#' @details ties are handled according to Efron
#' @return a vector of individual likelihoods (sorted by increasing event times)
#' @examples
#' \dontrun{
#' ### example data set from Fine paper, section 5
#' require("survival")
#' mod <- coxph(Surv(time, status==2) ~ age + log(albumin) + log(bili) + edema + log(protime), data=subset(pbc, !is.na(trt)),  x=T)
#' # individual LLs sum up to total model LL
#' sum(llcont(mod));logLik(mod)[1]
#' }
#' @importFrom data.table data.table
#' @export

llcont.coxph <- function(x) {

  if (is.null(x$x)) stop("coxph object without x=T option fitted")
  if (any(x$y[,"status"]>1)) stop("competing risk not supported")

  tmpdat          <- data.frame(time=x$y[,"time"], status=x$y[,"status"], elp=1)
  if(ncol(x$x)>0) tmpdat$elp <- exp(drop(x$x %*% coef(x)))
  tmpdat          <- tmpdat[order(-tmpdat$time),]
  tmpdat$cumelp   <- cumsum(tmpdat$elp)
  # ties handling (Efron)
  tmpdat          <- data.table(tmpdat)
  tmpdat          <- tmpdat[, cumelp := max(cumelp), by = time]
  tmpdat          <- tmpdat[, corr   := sum(elp[status==1]), by = time]
  tmpdat          <- tmpdat[, weight := pmax(0,cumsum(status)-1)/pmax(1,sum(status)), by = time]
  tmpdat$cumelp   <- tmpdat$cumelp - tmpdat$weight * tmpdat$corr
  # individual log-likelihoods
  tmpdat$lli      <- log(tmpdat$elp/tmpdat$cumelp) * tmpdat$status
  # order by time, needed for l3 term computation in plrtest
  tmpdat          <- tmpdat[order(tmpdat$time),]
  return(tmpdat$lli)
}


