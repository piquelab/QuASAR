#' Title fitQuasarMpra
#'
#' @description
#' For an MPRA experiment assess allele specific expression. 
#' 
#' @param ref reference read count
#' @param alt alternate read count
#' @param prop reference DNA proportion (def.=0.5)
#' @param eps base-calling error rate (def.=0.001)
#' @param nbreaks number of breaks of coverage to estimate dispersion (def.=10)
#'
#' @return data.frame with the following fields: 
#' - bin total coverage bin used
#' - betas.beta.binom  logit transfomation of the allelic skew
#' - betas_se standard error of beta 
#' - betas_z zscore for beta-plogis(propr) 
#' - pval3 p.value 
#' - padj_quasar BH adjusted p.value 
#' 
#' @export
#'
#' @examples
#' 
fitQuasarMpra <- function(ref,alt,prop=0.5,eps=0.001,nbreaks=10){
  tot <- ref + alt;
  cov_breaks <- unique(c(0,quantile(tot,(1:nbreaks)/nbreaks)))
  bin <- cut(tot,cov_breaks)
  M <- exp((0:500)/50)
  ## aux ~ loglikelihood of the beta bionmial model across D values
  ##               with null \rho value
  ## Mmax ~ disperison which maximizes the llk
  Mvec <- sapply(levels(bin),function(mybin){
    aux <- sapply(M,function(M){
      sum(logLikBetaBinomialRhoEps(prop[bin==mybin],eps,M,ref[bin==mybin ],alt[bin==mybin]))
    })
    Mmax <- M[which.max(aux)]
    Mmax
  })
  #Mvec
  #cat("#Disp: ", round(Mvec, 2), "\n")
  aux2 <- t(sapply(1:length(ref),function(ii){
    auxLogis <- optim(0,
                      fn=logLikBetaBinomial2,
                      gr=gLogLikBetaBinomial,
                      D=Mvec[bin[ii]],
                      R=ref[ii],
                      A=alt[ii],
                      method="L-BFGS-B",
                      hessian=TRUE,
                      lower=-5,
                      upper=5)
    c(auxLogis$par,1/(auxLogis$hessian)^.5)
  }))
  ## eps ~ jointly inferred error rate for this sample
  rho3 <- plogis(aux2[,1])
  betas.beta.binom <- aux2[,1]
  lrt3 <-
    logLikBetaBinomialRhoEps(rho3,eps,Mvec[bin],ref,alt) -
    logLikBetaBinomialRhoEps(prop,eps,Mvec[bin],ref,alt)
  pval3 <- (1-pchisq(2*lrt3,df=1))
  ##
  ##betas.se <- abs(betas.beta.binom/qnorm(pval3/2))
  betas_se <- aux2[, 2]
  ##  betas_w <- 1/betas_se^2
  betas_z <- (betas.beta.binom - qlogis(prop))/betas_se
  ##
  padj_quasar <- p.adjust(pval3,method="BH")
  data.frame(bin,betas.beta.binom,betas_se,betas_z,pval3,padj_quasar)
}
