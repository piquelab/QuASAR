#' @title gradient of log-likelihood of beta-binomial model
#'
#' @description
#' \code(gLogLikBetaBinomial) returns the gradient of the log likelihood of reference & alternate 
#' read count data given logit(\rho) and adispersion estimate.
#'
#' @param logit.p logit(\rho) 
#' @param D dispersion estimate
#' @param R reference read count
#' @param A alternate read count
#' @return gradient evaulated at logit(\rho)
#'
gLogLikBetaBinomial <- function(logit.p,D,R,A){
  p <- plogis(logit.p)
 # aux <- D * (digamma(R+p*D) - digamma(A + (1-p)*D) - digamma(p*D) + digamma((1-p)*D)) 
  aux <- (p * (1-p))*D*(digamma(R+p*D) - digamma(A + (1-p)*D) - digamma(p*D) + digamma((1-p)*D)) 
  as.matrix(-aux)
}
