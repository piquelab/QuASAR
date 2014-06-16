#' @title log-likelihood of beta-binomial model for heterozygotes
#'
#' @description
#' \code(logLikBetaBinomialRhoEps) returns the log likelihood of reference & alternate read count data
#' given rho, & dispersion.
#'
#' @param rho rho value \in [0, 1]
#' @param eps estimate of error
#' @param D dispersion estimate
#' @param R reference read count
#' @param A alternate read count
#' @return log likelihood
logLikBetaBinomial2 <- function(logit.p,D,R,A){
  p <- plogis(logit.p)
  aux <- (lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) + lgamma(D));
  -sum(aux)
}
