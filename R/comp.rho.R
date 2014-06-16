#' @title rho calculation
#'
#' @description
#' given reference and alternate read counts, & an error estimate, return the 
#' estimated alternate alle frequency
#'
#' @param ref integer. reference allele count 
#' @param alt integer. alternate allele count 
#' @param eps float. error estimate in [0,1]
#' @return rho float. elstimate of alternate allele frequency
comp.rho <- function(ref,alt,eps){
  rho <- (ref * (1 - eps) - eps * alt) / ((ref + alt) * (1 - 2 * eps))
  rho[rho<0] <- 0
  rho[rho>1] <- 1
  rho
}	
