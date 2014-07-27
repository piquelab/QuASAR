#' @title critical value for likelihood ratio test with misspecififed genotypes
#'
#' @description
#' returns the critical value of llrt comparing with 
#' H0: the het is mis-specified or rho=o.5 & H1: the het has imbalance with estimated rho
#'
#' @param ref reference read count
#' @param alt alternate read count
#' @param eps estimate of error
#' @param rho rho estimate
#' @return llrt critical value
#' @export
lrtEpsRhoBinom <- function(ref,alt,eps,rho){
  log.gt <- cbind(g0   = ref * log(1 - eps) + alt * log(eps),
                  g1t0 = (ref + alt) * log(0.5),
                  g1t1 = ref * log((1 - eps) * rho + eps * (1 - rho)) + alt * log((1 - eps) * (1 - rho) + eps * rho),
                  g2   = ref * log(eps) + alt * log(1 - eps)
                  )
  lrt <- log.gt[,"g1t1"]-pmax(log.gt[,"g0"],log.gt[,"g1t0"],log.gt[,"g2"])
}
