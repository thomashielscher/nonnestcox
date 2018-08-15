#' Getting partial log-likelihood of a \code{coxph} object for individual cases
#'
#' @author Thomas Hielscher
#' @param model \code{coxph} model with \code{x=T}
#' @details ties are handled according to Efron
#' @return a vector of individual likelihoods
#' @examples
#' \dontrun{
#' ### example data set from Fine paper, section 5
#' require("survival")
#' mod <- coxph(Surv(time, status==2) ~ age + log(albumin) + log(bili) + edema + log(protime), data=subset(pbc, !is.na(trt)),  x=T)
#' # individual LLs sum up to total model LL
#' sum(llcont(mod));logLik(mod)[1]
#' }
#' @importFrom plyr ddply ldply llply mutate .
#' @export

llcont <- function(x, ...) UseMethod("llcont")

#' @export

llcont.coxph <- function(object) {

  if (is.null(object$x)) stop("coxph object without x=T option fitted")
  # save original ordering of data before calculations
  tmpdat         <- data.frame(time=object$y[,1], status=object$y[,2], elp=exp(drop(object$x %*% coef(object))), ordering=1:nrow(tmpdat))
  tmpdat         <- tmpdat[order(tmpdat$time),]
  tmpdat$cumelp  <- rev(cumsum(rev(tmpdat$elp)))

  for (i in 2:nrow(tmpdat)) if(tmpdat$time[i]==tmpdat$time[i-1]) tmpdat$cumelp[i] <- tmpdat$cumelp[i-1]

  tmpdat         <- ddply(tmpdat, .(time), mutate, corr=sum(elp[status==1]), weight=pmax(0,cumsum(status)-1)/pmax(1,sum(status)))
  tmpdat$cumelp  <- tmpdat$cumelp - tmpdat$weight * tmpdat$corr
  tmpdat$lli     <- ifelse(tmpdat$status==1, log(tmpdat$elp/tmpdat$cumelp), 0)
  # restore original ordering
  tmpdat         <- tmpdat[order(tmpdat$ordering),]
  return(tmpdat$lli)
}


