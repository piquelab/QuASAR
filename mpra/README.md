# QuASAR-MPRA: Accurate allele-specific analysis for massively parallel  reporter assays
We have further developed our method for allele specific analysis QuASAR (quantitative allele-specific analysis of reads) ([Harvey et al, 2015]) to analyze allele specific signals in barcoded read counts data from MPRAs. Using this approach, we can take into account the uncertainty on the original plasmid proportions, over-dispersion, and sequencing errors. Here, we demonstrate how to use QuASAR-MPRA to analyze the MPRA data by [Tewhey et al,2016].

The current software is still in development and we will kindly appreciate any comments and bug reports. We assume that you already have installed the QuASAR library. The [mpra.R] script contains the instructions to run the test, and [preprocessing.R] contains the steps we did to preprocess the HepG2 data for the forward strand in [Tewhey et al,2016] in this example. 

```R
library(QuASAR)

## Loading the sample data:
HepG2 <- read.table("http://genome.grid.wayne.edu/quasar/samplempra/HepG2.mpra.txt",sep='\t',as.is=T,header=T)

## Fitting the QuASAR model:
HepG2.res <- fitQuasarMpra(HepG2$R,HepG2$A,HepG2$DNA_prop)

## Number of significant hits 10% FDR:
sum(HepG2.res$padj_quasar<0.1)

## QQ-plot: 
library(qqman)
qq(HepG2.res$pval3)

```

To fit the model we use the `fitQuasarMpra` function that needs three vectors: `ref` = number of RNA reads for the reference allele, `alt` = number of RNA reads for the alternate allele, `prop` the reference DNA proportion in the plasmid library. 

<!-- links -->
[Harvey et al, 2015]:http://bioinformatics.oxfordjournals.org/content/31/8/1235
[mpra.R]:mpra.R
[Tewhey et al,2016]:https://www.ncbi.nlm.nih.gov/pubmed/27259153

