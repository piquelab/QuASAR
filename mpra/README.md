# QuASAR-MPRA: Accurate allele-specific analysis for massively parallel  reporter assays
We have further developed our method for allele specific analysis QuASAR (quantitative allele-specific analysis of reads) ([Harvey et al, 2015]) to analyze allele specific signals in barcoded read counts data from MPRAs. Using this approach, we can take into account the uncertainty on the original plasmid proportions, over-dispersion, and sequencing errors. Here, we demonstrate how to use QuASAR-MPRA to analyze the MPRA data by [Tewhey et al,2016].

The current software is still in development and we will kindly appreciate any comments and bug reports. We assume that you already have installed the QuASAR library. 

The [mpra.R] script contains the instructions to run the test. As an example we provide the HepG2 data for the forward strand in [Tewhey et al,2016].[preprocessing.R] contains the steps to prepare this input file.

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

To fit the model we use the `fitQuasarMpra` function that needs three input vectors: 
- `ref` = number of RNA reads for the reference allele 
- `alt` = number of RNA reads for the alternate allele 
- `prop` = reference DNA proportion in the plasmid library

The returned data frame `HepG2.res` has the following fileds:
- `bin` total coverage bin used
- `betas.beta.binom`  logit transfomation of the RNA allelic skew
- `betas_se` standard error for the beta parameter estimate
- `betas_z` zscore for `betas.beta.binom` -plogis(`propr`) 
- `pval3` p.value 
- `padj_quasar` BH adjusted p.value 

```R
> head(HepG2.res)
                  bin betas.beta.binom  betas_se    betas_z       pval3 padj_quasar
1   (2.7e+03,3.3e+03]      -0.12164712 0.1409832 -0.2849424 0.774655868   1.0000000
2   (2.1e+03,2.7e+03]      -0.22659942 0.1630645 -1.0884486 0.273956397   1.0000000
3  (3.3e+03,3.96e+03]      -0.05672033 0.1282004 -0.6543586 0.513012947   1.0000000
4 (7.56e+03,2.22e+06]      -0.46219112 0.2338646 -0.1999347 0.838163133   1.0000000
5   (2.7e+03,3.3e+03]      -0.04063415 0.1403660  0.8889998 0.376112030   1.0000000
6             (0,859]      -1.17993132 0.5370852 -2.8797701 0.001298818   0.1951728
```

<!-- links -->
[Harvey et al, 2015]:http://bioinformatics.oxfordjournals.org/content/31/8/1235
[mpra.R]:mpra.R
[process.R]:process.R
[Tewhey et al,2016]:https://www.ncbi.nlm.nih.gov/pubmed/27259153

