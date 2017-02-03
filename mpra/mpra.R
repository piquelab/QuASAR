## QuASAR ############################

library(QuASAR)
library(qqman)

## Import and preprocess the data by first running preprocess.R
## source(preprocess.R)

#Hepg2  Forward
dd <- unique(mpra_Filt[,-c(2:11)])
aux <- fitQuasarMpra(dd$R,dd$A,dd$DNA_prop)
zzz <- cbind(dd,aux)

#Hepg2  Reverse
dd <- unique(RC_Filt[,-c(2:11)])
aux <- fitQuasarMpra(dd$R,dd$A,dd$DNA_prop)
zzz_RC <- cbind(dd,aux)

## Number of significant hits below 0.1
sum(zzz$padj_quasar<0.1)
sum(zzz_RC$padj_quasar<0.1)

## QQ-plot 
qq(zzz$pval3)

qq(zzz_RC$pval3)


