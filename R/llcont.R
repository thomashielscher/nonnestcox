#' Getting partial log-likelihood of a \code{coxph} object for individual cases
#'
#' @author Thomas Hielscher
#' @param model \code{coxph} model with \code{x=T}
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

llcont <- function(x, ...) UseMethod("llcont")

#' @export

llcont.coxph <- function(object) {

  if (is.null(object$x)) stop("coxph object without x=T option fitted")
  if (any(object$y[,"status"]>1)) stop("competing risk not supported")

  tmpdat          <- data.frame(time=object$y[,"time"], status=object$y[,"status"], elp=exp(drop(object$x %*% coef(object))))
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
  return(tmpdat$lli)
}


