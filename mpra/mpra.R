#processing
#DNA library 
#summed replicate counts
#require ref+alt counts >100 
#require 5 counts each ref and alt
#DNA_prop = ref/(ref+alt)
#RNA library
#require being in DNA library
#require 5 counts each ref and alt

## QuASAR ############################

##It should be in an RData object, /wsu/home/fh/fh85/fh8591/piquelab/cindy/bitStarr/Tewhey/RNA/hepg2/summed.RNA.RData, LCL ind 1: ~/piquelab/cindy/bitStarr/Tewhey/summed.RNA.base.RData , and LCL ind 2: ~/piquelab/cindy/bitStarr/Tewhey/summed.RNA.base.LCL2.RData

##Loading LCL.1
##load("~/piquelab/cindy/bitStarr/Tewhey/summed.RNA.base.RData")


library(QuASAR)
library(qqman)


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

#individual 1 forward
dd <- unique(mpra_Filt[,-c(3:12)])
aux <- fitQuasarMpra(dd$R,dd$A,dd$DNA_prop)
zzz <- cbind(dd,aux)

#individual 1 reverse
dd <- unique(RC_Filt[,-c(2:11)])
aux <- fitQuasarMpra(dd$R,dd$A,dd$DNA_prop)
zzz_RC <- cbind(dd,aux)

#individual 2 forward
dd <- unique(mpra_2_Filt[,-c(2:7,10:11)])
aux <- fitQuasarMpra(dd$R,dd$A,dd$DNA_prop)
zzz_2 <- cbind(dd,aux)

# individual 2 reverse
dd <- unique(RC_2_Filt[,-c(2:7,10:11)])
aux <- fitQuasarMpra(dd$R,dd$A,dd$DNA_prop)
zzz_RC_2 <- cbind(dd,aux)


#combining individuals using fixed effects method
LCL <- merge(zzz[,-c(6:8,12)], zzz_2[,-c(4:6,10)], by=c("rsID","alt_hap"), all=T)
RC_LCL <- merge(zzz_RC[,-c(6:8,12)], zzz_RC_2[,-c(4:6,10)], by=c("rsID","alt_hap"), all=T)
LCL1 <- LCL[!is.na(LCL$betas_w.x), ]
LCL1 <- LCL1[!is.na(LCL1$betas_w.y), ]
LCL1 <- unique(LCL1)
RC_LCL1 <- RC_LCL[!is.na(RC_LCL$betas_w.x), ]
RC_LCL1 <- RC_LCL1[!is.na(RC_LCL1$betas_w.y), ]
RC_LCL1 <- unique(RC_LCL1)

#no DNA NA, removed
LCL <- transform(LCL1, DNA_prop=ifelse(DNA_prop.x=="NA", DNA_prop.y, DNA_prop.x), betas_T=((betas_w.x*betas.beta.binom.x) + (betas_w.y*betas.beta.binom.y))/(betas_w.x+betas_w.y), 
                 betas_se_comb = sqrt(1/(betas_w.x+betas_w.y)))
LCL <- transform(LCL, betas_T_low = betas_T-(1.96*betas_se_comb), 
                 betas_T_high = betas_T+(1.96*betas_se_comb),
                 betas_z_comb = (betas_T - qlogis(DNA_prop)) /betas_se_comb)
RC_LCL <- transform(RC_LCL1, DNA_prop=ifelse(DNA_prop.x=="NA", DNA_prop.y, DNA_prop.x), betas_T=((betas_w.x*betas.beta.binom.x) + (betas_w.y*betas.beta.binom.y))/(betas_w.x+betas_w.y), 
                    betas_se_comb = sqrt(1/(betas_w.x+betas_w.y)))
RC_LCL <- transform(RC_LCL, betas_T_low = betas_T-(1.96*betas_se_comb), 
                    betas_T_high = betas_T+(1.96*betas_se_comb),
                    betas_z_comb = (betas_T - qlogis(DNA_prop)) /betas_se_comb)

LCL <- transform(LCL, p_comb = 2 * pnorm(-abs(betas_z_comb)))
RC_LCL <- transform(RC_LCL, p_comb = 2 * pnorm(-abs(betas_z_comb)))

#genomic inflation test (grabbed this from hepg2 script)
require(IDPmisc) 
#NaRV.omit in IDPmisc package
#paper ref: http://www.nature.com/ejhg/journal/v19/n7/full/ejhg201139a.html
#http://stats.stackexchange.com/questions/110755/how-calculate-inflation-observed-and-expected-p-values-from-uniform-distribution
#test for genomic inflation
chisq <- qchisq(1-RC_hepg2_MPRA$pval3,1)
chisq <- NaRV.omit(chisq)
#ks <- ks.test(RC_hepg2_MPRA$pval3, "punif", 0, 1) another test optional
# For z-scores as association, just square them
# chisq <- data$z^2
#For chi-squared values, keep as is
#chisq <- data$chisq
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 0.7474663
chisq <- qchisq(1-RC_hepg2_MPRA$C.Skew.logP,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 1.546896
chisq <- qchisq(1-RC_hepg2_MPRA$pval_binom,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 27.43236
chisq <- qchisq(1-RC_hepg2_MPRA$pval_ttest,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 2.725156
chisq <- qchisq(1-RC_hepg2_MPRA$pval_fisher,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda
#[1] 10.00775

#test for genomic inflation reverse
chisq <- qchisq(1-hepg2_MPRA$pval3,1)
chisq <- NaRV.omit(chisq)
# For z-scores as association, just square them
# chisq <- data$z^2
#For chi-squared values, keep as is
#chisq <- data$chisq
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 0.7636204
chisq <- qchisq(1-hepg2_MPRA$C.Skew.logP,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 1.655326
chisq <- qchisq(1-hepg2_MPRA$pval_binom,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 27.06694
chisq <- qchisq(1-hepg2_MPRA$pval_ttest,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda 
#[1] 2.850532
chisq <- qchisq(1-hepg2_MPRA$pval_fisher,1)
chisq <- NaRV.omit(chisq)
lambda = median(chisq)/qchisq(0.5,1)
lambda
#[1] 9.945027
