## QuASAR ############################

library(QuASAR)
library(qqman)

## Import and preprocess the data by first running preprocess.R
## Here we already have it saved as a text file as example
HepG2 <- read.table("http://genome.grid.wayne.edu/quasar/samplempra/HepG2.mpra.txt",sep='\t',as.is=T,header=T)

## Fitting the QuASAR model
HepG2.res <- fitQuasarMpra(HepG2$R,HepG2$A,HepG2$DNA_prop)

## Number of significant hits 10% FDR
sum(HepG2.res$padj_quasar<0.1)

## QQ-plot 
qq(HepG2.res$pval3)



