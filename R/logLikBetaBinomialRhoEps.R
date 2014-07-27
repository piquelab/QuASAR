#' @title log-likelihood of beta-binomial model
#'
#' @description
#' returns the log likelihood of reference & alternate read count data
#' given rho, dispersion, and error
#'
#' @param rho rho value in [0, 1]
#' @param eps estimate of error
#' @param D dispersion estimate
#' @param R reference read count
#' @param A alternate read count
#' @return log likelihood
#' @export
logLikBetaBinomialRhoEps <- function(rho,eps,D,R,A){
  p <- (rho*(1-eps)+(1-rho)*eps)
  aux <- (lgamma(R+p*D) + lgamma(A + (1-p)*D) - lgamma(R+A+D) - lgamma(p*D) - lgamma((1-p)*D) + lgamma(D)) ##+ lgamma(R+A+1) - lgamma(A+1) - lgamma(R+1)
  aux
}
