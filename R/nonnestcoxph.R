################################################################
## print method for Fine test
################################################################
#' @method print finetest
#' @export

print.finetest <- function (x, ...)
{
  cat("\nVariance test \n")
  cat("  H0: Model 1 and Model 2 are indistinguishable", "\n")
  cat("  H1: Model 1 and Model 2 are distinguishable", "\n")
  cat("Fine: p = ", format.pval(x$pOmega1, digits = 3L), "\n", sep = "")
  cat("Vuong: p = ", format.pval(x$pOmega2, digits = 3L), "\n\n", sep = "")

  if (x$nested) {
    cat("Robust likelihood ratio test of distinguishable models \n")
    cat("  H0: Both models fit equally well \n")
    cat("  H1: Full model fits better than reduced model \n")
    cat("    LR = ", formatC(x$LRTstat, digits = 3L, format = "f"),
        ",   ", "p = ", format.pval(x$pLRTAB, digits = 3L),
        "\n", sep = "")
    cat("Classical likelihood ratio test\n")
    cat("  H0: Both models fit equally well \n")
    cat("  H1: Full model fits better than reduced model \n")
    cat("    LR = ", formatC(x$LRTstat, digits = 3L, format = "f"),
        ",   ", "p = ", format.pval(x$pLRT, digits = 3L),
        "\n", sep = "")

  }
  else {
    cat("Non-nested likelihood ratio test \n")
    cat("  H0: Model fits are equally close to true Model \n")
    cat("  H1A: Model 1 fits better than Model 2 \n")
    cat("    z = ", formatC(x$LRTstat, digits = 3L, format = "f"),
        ",   ", "p = ", format.pval(x$pLRTA, digits = 3L),
        "\n", sep = "")
    cat("  H1B: Model 2 fits better than Model 1 \n")
    cat("    z = ", formatC(x$LRTstat, digits = 3L, format = "f"),
        ",   ", "p = ", format.pval(x$pLRTB, digits = 4L),
        "\n", sep = "")
  }
}
